"""EXPERIMENTAL: insertions in the region upstream of a gene.

Not part of the default MISFIT run. MISFIT proper reads coding sequence, so a
gene whose promoter carries an IS element is reported — correctly — as an
intact ORF. This module looks at the other half of that picture.

The reference here is the cohort, not an external genome: intergenic sequence
diverges much faster than coding sequence, but across a few hundred isolates of
one species the modal upstream architecture is the wild type. For each
(panel, gene) the module picks a conserved anchor upstream of the start codon,
measures the anchor-to-ATG distance in every isolate, and treats the mode as
wild type. Isolates carrying substantially more sequence than the mode have an
insertion, which is then localised by aligning against a wild-type isolate.

The module reports the size of the insertion and how far upstream of the start
codon it sits. It deliberately does not decide whether that is disruptive --
distance matters and the cutoff is a judgement call, so the numbers are left
for downstream analysis to weigh.

usage:
    python -m misfit.scripts.upstream <misfit.tsv> <assembly_dir> <out.tsv>
        [--species-map map.tsv] [--genes ompC,ompF] [--window 3000]
        [--min-insertion 100]
"""

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

from Bio import Align

from .misc import load_fasta, reverse_complement

ANCHOR_K = 40
# Offsets upstream of the start codon to try as anchor positions, far first: an
# anchor beyond the insertion point is what makes the extra sequence visible.
ANCHOR_OFFSETS = (2500, 2000, 1600, 1300, 1000, 800, 600, 450, 350)
MIN_ANCHOR_SUPPORT = 0.35
# The window must clear the furthest anchor offset by more than the biggest
# insertion we expect: an IS pushes the anchor further from the start codon, and
# an anchor shoved past the window edge reads as "not found" instead of a hit.
DEFAULT_WINDOW = 8000
# Exact k-mer anchoring is fast but a single SNP destroys an anchor, and some
# upstream regions are too variable to anchor at all. Those isolates fall back
# to aligning against a wild-type representative over this shorter window.
ALIGN_WINDOW = 3000
# Reference-free repeat scan. Where the cohort has no agreed wild-type upstream
# architecture, comparing to a mode is meaningless -- but a mobile element still
# betrays itself by occurring many times in the isolate's own assembly.
TILE = 80
TILE_STEP = 150
REPEAT_SCAN_BP = 2400

_ALIGNER = Align.PairwiseAligner(
    mode="global", match_score=2, mismatch_score=-3,
    open_gap_score=-8, extend_gap_score=-0.5,   # cheap extension: favour one long gap
)


def upstream_window(contig_seq, start, end, strand, window):
    """Sequence immediately 5' of the CDS, oriented so it ends at the start codon."""
    if strand == "+":
        lo = max(0, start - 1 - window)
        return contig_seq[lo:start - 1], (start - 1) >= window
    hi = min(len(contig_seq), end + window)
    return reverse_complement(contig_seq[end:hi]), (len(contig_seq) - end) >= window


def pick_anchors(windows, max_anchors=5):
    """Choose several conserved, single-copy anchors at different offsets.

    One anchor is too brittle: a panel spanning several species, or an isolate
    whose insertion landed on the anchor itself, loses it and the isolate drops
    out of the analysis entirely. Holding a ladder of anchors lets each isolate
    use the furthest one it actually has.
    """
    n = len(windows)
    found = []
    for off in ANCHOR_OFFSETS:
        counts = Counter()
        for w in windows.values():
            if len(w) >= off:
                kmer = w[-off:-off + ANCHOR_K]
                if len(kmer) == ANCHOR_K and "N" not in kmer:
                    counts[kmer] += 1
        if not counts:
            continue
        kmer, support = counts.most_common(1)[0]
        if support < MIN_ANCHOR_SUPPORT * n:
            continue
        # An anchor that repeats within a window cannot place anything. Tolerate a
        # few such isolates rather than discarding an otherwise good anchor.
        repeated = sum(1 for w in windows.values() if w.count(kmer) > 1)
        if repeated > 0.05 * n:
            continue
        found.append((kmer, off, support))
        if len(found) >= max_anchors:
            break
    return found


def wt_representative(windows, span):
    """Pick a wild-type window without needing an anchor.

    Most isolates are wild type, so the commonest sequence immediately 5' of
    the start codon is the native promoter. An isolate carrying an insertion
    close in has a different proximal sequence and is excluded by construction.
    """
    proximal = Counter()
    for w in windows.values():
        if len(w) >= 200:
            proximal[w[-200:]] += 1
    if not proximal:
        return None, None
    seq, n = proximal.most_common(1)[0]
    if n < 3:
        return None, None
    for asm, w in sorted(windows.items()):
        if len(w) >= 200 and w[-200:] == seq:
            return asm, w[-span:]
    return None, None


def locate_insertion(query_win, wt_win):
    """Align an outlier's upstream region to a wild-type one.

    Returns (insertion_length, distance_upstream_of_start_codon, inserted_seq)
    for the largest gap, measuring distance from the 3' end of the insertion to
    the start codon -- i.e. how far upstream of the ATG the element sits.
    """
    aln = _ALIGNER.align(wt_win, query_win)[0]
    ref_aln, qry_aln = aln[0], aln[1]

    best = (0, None, "")
    q_consumed = 0
    i = 0
    n = len(ref_aln)
    while i < n:
        if ref_aln[i] == "-":
            start_q = q_consumed
            seq = ""
            while i < n and ref_aln[i] == "-":
                seq += qry_aln[i]
                q_consumed += 1
                i += 1
            if len(seq) > best[0]:
                best = (len(seq), start_q + len(seq), seq)
        else:
            if qry_aln[i] != "-":
                q_consumed += 1
            i += 1
    if best[1] is None:
        return 0, None, ""
    # distance from the end of the insertion to the start codon
    return best[0], len(query_win) - best[1], best[2]


def repeat_copies(genome, genome_rc, probe):
    """How many times a sequence occurs in the assembly. >1 marks a mobile element."""
    if len(probe) < 60:
        return None
    core = probe[len(probe) // 2 - 40: len(probe) // 2 + 40]
    return genome.count(core) + genome_rc.count(core)


def scan_upstream_repeats(window, genome, genome_rc, min_copies=2):
    """Find multi-copy sequence upstream, without reference to any wild type.

    Tiles the region 5' of the start codon and counts each tile in the isolate's
    own assembly. A tile present several times is repeated sequence -- an IS
    element, most often. Reports the multi-copy tile nearest the start codon and
    how far upstream it sits. Repeated does not by itself mean mobile: rRNA
    operons and REP elements also recur, so the copy count is reported rather
    than interpreted.
    """
    nearest, best_copies, hits = None, 0, []
    for off in range(TILE, min(REPEAT_SCAN_BP, len(window)) + 1, TILE_STEP):
        tile = window[-off:-off + TILE]
        if len(tile) < TILE or "N" in tile:
            continue
        n = genome.count(tile) + genome_rc.count(tile)
        if n >= min_copies:
            hits.append(off)
            best_copies = max(best_copies, n)
            if nearest is None:
                nearest = off
    span = (max(hits) - min(hits) + TILE) if hits else 0
    return nearest, best_copies, span


def collect_windows(rows, assembly_dir, fasta_for, genes, window):
    """Extract upstream windows, loading each assembly exactly once."""
    by_assembly = defaultdict(list)
    for r in rows:
        if genes and r["Gene"] not in genes:
            continue
        if r["Mutation_Type"] in ("Not found", "Error") or r["Contig"] in ("--", ""):
            continue
        by_assembly[r["Assembly"]].append(r)

    windows, genomes, genomes_rc, truncated = {}, {}, {}, set()
    for asm, items in sorted(by_assembly.items()):
        fname = fasta_for.get(asm)
        path = os.path.join(assembly_dir, fname) if fname else None
        if not path or not os.path.exists(path):
            continue
        seqs = dict(load_fasta(path))
        g = "".join(seqs.values())
        genomes[asm] = g
        genomes_rc[asm] = reverse_complement(g)
        for r in items:
            contig = seqs.get(r["Contig"])
            if contig is None:
                continue
            w, full = upstream_window(contig, int(r["Contig_Start"]),
                                      int(r["Contig_End"]), r["Strand"], window)
            windows[(asm, r["Gene"])] = w
            if not full:
                truncated.add((asm, r["Gene"]))
    return windows, genomes, genomes_rc, truncated


def main():
    ap = argparse.ArgumentParser(
        description="EXPERIMENTAL: detect insertions upstream of genes called by MISFIT. "
                    "Run separately; not part of the standard MISFIT pipeline.")
    ap.add_argument("misfit_tsv", help="MISFIT output (must include Contig/Contig_Start/Strand)")
    ap.add_argument("assembly_dir")
    ap.add_argument("out")
    ap.add_argument("--species-map", default=None,
                    help="TSV with fasta_file and species_wgs, to group isolates by panel")
    ap.add_argument("--genes", default=None, help="comma-separated gene subset")
    ap.add_argument("--window", type=int, default=DEFAULT_WINDOW)
    ap.add_argument("--min-insertion", type=int, default=100,
                    help="report insertions at least this many bp longer than the cohort mode")
    args = ap.parse_args()

    print("[misfit-upstream] EXPERIMENTAL module. Reports insertion size and distance "
          "upstream of the start codon; it does not judge whether an insertion is "
          "disruptive.", file=sys.stderr)

    rows = list(csv.DictReader(open(args.misfit_tsv), delimiter="\t"))
    if "Contig_Start" not in (rows[0] if rows else {}):
        sys.exit("error: input lacks Contig_Start; re-run MISFIT to get position columns")

    fasta_for, panel_of = {}, {}
    if args.species_map:
        sys.path.insert(0, os.path.dirname(os.path.abspath(args.species_map)))
        for r in csv.DictReader(open(args.species_map), delimiter="\t"):
            stem = r["fasta_file"].split(".")[0]
            fasta_for[stem] = r["fasta_file"]
            panel_of[stem] = r.get("species_wgs", "-")
    else:
        for f in os.listdir(args.assembly_dir):
            if f.endswith((".fasta", ".fa", ".fna")):
                fasta_for.setdefault(f.split(".")[0], f)

    # rows may already carry a panel column from the frequency annotation
    if rows and "panel" in rows[0]:
        for r in rows:
            panel_of[r["Assembly"]] = r["panel"]

    genes = set(args.genes.split(",")) if args.genes else None
    windows, genomes, genomes_rc, truncated = collect_windows(
        rows, args.assembly_dir, fasta_for, genes, args.window)
    print(f"[misfit-upstream] extracted {len(windows)} upstream windows", file=sys.stderr)

    groups = defaultdict(dict)
    for (asm, gene), w in windows.items():
        groups[(panel_of.get(asm, "-"), gene)][asm] = w

    out_rows = []
    for (panel, gene), wins in sorted(groups.items()):
        if len(wins) < 5:
            for asm in wins:
                out_rows.append(dict(assembly=asm, panel=panel, gene=gene,
                                     status="too_few_isolates", spacing_bp="",
                                     cohort_mode_bp="", insertion_bp="",
                                     distance_upstream_bp="", repeat_copies="",
                                     note=f"only {len(wins)} isolate(s) in this panel/gene"))
            continue

        anchors = pick_anchors(wins)
        per_anchor = []

        # spacing and wild-type mode are anchor-specific, so keep them separate
        for kmer, off, _support in anchors:
            sp = {}
            for asm, w in wins.items():
                i = w.find(kmer)
                if i != -1:
                    sp[asm] = len(w) - i
            if not sp:
                continue
            m = Counter(sp.values()).most_common(1)[0][0]
            wt = min((a for a, v in sp.items() if v == m), default=None)
            per_anchor.append({"offset": off, "spacing": sp, "mode": m,
                               "wt_win": wins[wt][-m:] if wt else None})
        per_anchor.sort(key=lambda a: -a["offset"])

        def measure(asm):
            """Furthest anchor this isolate actually carries."""
            for a in per_anchor:
                if asm in a["spacing"]:
                    return a
            return None

        mode = per_anchor[0]["mode"] if per_anchor else ""
        fb_asm, fb_win = wt_representative(wins, ALIGN_WINDOW)

        # How much do these isolates agree on their upstream architecture? Where
        # they do not, the cohort mode means nothing and the reference-free
        # repeat scan is the only usable signal -- so pay for it there.
        prox = Counter(w[-200:] for w in wins.values() if len(w) >= 200)
        agreement = (prox.most_common(1)[0][1] / len(wins)) if prox else 0.0
        scan_group = agreement < 0.5

        for asm, w in wins.items():
            base = dict(assembly=asm, panel=panel, gene=gene, cohort_mode_bp=mode,
                        anchor_offset_bp="", method="", insertion_bp="",
                        distance_upstream_bp="", repeat_copies="", note="")
            unplaced = measure(asm) is None
            if scan_group or unplaced:
                near, ncopies, nspan = scan_upstream_repeats(
                    w, genomes.get(asm, ""), genomes_rc.get(asm, ""))
                base["repeat_nearest_bp"] = near if near is not None else ""
                base["repeat_max_copies"] = ncopies or ""
                base["repeat_span_bp"] = nspan or ""
                base["cohort_agreement"] = f"{agreement:.2f}"
            else:
                base["repeat_nearest_bp"] = base["repeat_max_copies"] = ""
                base["repeat_span_bp"] = ""
                base["cohort_agreement"] = f"{agreement:.2f}"
            a = measure(asm)
            if a is None:
                short = (asm, gene) in truncated and len(w) < (mode or 0) + args.min_insertion
                if short or fb_win is None or asm == fb_asm:
                    base.update(
                        status="not_assessable" if short else
                                ("wild_type_reference" if asm == fb_asm else "anchor_not_found"),
                        spacing_bp="", method="",
                        note=(f"only {len(w)} bp of contig upstream"
                              if short else
                              "used as the wild-type reference for this group"
                              if asm == fb_asm else
                              "no anchor and no wild-type representative available"))
                    out_rows.append(base)
                    continue
                length, dist, seq = locate_insertion(w[-ALIGN_WINDOW:], fb_win)
                copies = repeat_copies(genomes.get(asm, ""), genomes_rc.get(asm, ""), seq) if seq else None
                base.update(
                    status="insertion" if length >= args.min_insertion else "no_insertion",
                    spacing_bp="", method="alignment_vs_wild_type",
                    insertion_bp=length if length >= args.min_insertion else "",
                    distance_upstream_bp=dist if length >= args.min_insertion else "",
                    repeat_copies=copies if copies is not None else "",
                    note=("multi-copy in this assembly, consistent with a mobile element"
                          if copies and copies > 1 else
                          "single copy in this assembly" if copies == 1 else ""))
                out_rows.append(base)
                continue
            base["method"] = "anchor_spacing"

            s = a["spacing"][asm]
            mode, wt_win = a["mode"], a["wt_win"]
            base["spacing_bp"] = s
            base["cohort_mode_bp"] = mode
            base["anchor_offset_bp"] = a["offset"]
            if s - mode < args.min_insertion:
                base["status"] = "no_insertion" if s >= mode - args.min_insertion else "shorter_than_mode"
                if base["status"] == "shorter_than_mode":
                    base["insertion_bp"] = s - mode
                out_rows.append(base)
                continue

            length, dist, seq = locate_insertion(w[-s:], wt_win)
            copies = repeat_copies(genomes.get(asm, ""), genomes_rc.get(asm, ""), seq) if seq else None
            base.update(status="insertion", insertion_bp=length or (s - mode),
                        distance_upstream_bp=dist if dist is not None else "",
                        repeat_copies=copies if copies is not None else "",
                        note=("multi-copy in this assembly, consistent with a mobile element"
                              if copies and copies > 1 else
                              "single copy in this assembly" if copies == 1 else ""))
            out_rows.append(base)

    cols = ["assembly", "panel", "gene", "status", "spacing_bp", "cohort_mode_bp",
            "anchor_offset_bp", "method", "insertion_bp", "distance_upstream_bp",
            "repeat_copies", "repeat_nearest_bp", "repeat_max_copies",
            "repeat_span_bp", "cohort_agreement", "note"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    tally = Counter(r["status"] for r in out_rows)
    print(f"[misfit-upstream] {dict(tally)}", file=sys.stderr)
    print(f"[misfit-upstream] wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
