#!/usr/bin/env python3
"""Add alleles, genes and panels to the MISFIT reference set.

Every header in the reference panels follows one format:

    >gene|protein_accession|nuccore_accession|locus_tag

and every sequence has a row in PROVENANCE.tsv. That discipline is what makes a
call auditable, so this tool will not invent an identifier. It resolves the
header from a real record by one of these routes, in decreasing authority:

  --from-assembly GCF_...   Look the gene up in that RefSeq assembly's
                            annotation. All four header fields come straight
                            from the record, and the supplied sequence is
                            checked against it. Use this when you can.

  --from-nuccore CP_/NC_..  The same, from a single nucleotide record.

  --protein-accession WP_.. Verify the accession exists at NCBI and confirm it
                            translates to the supplied sequence.

  (sequence only)           BLAST the translated protein against RefSeq and
                            print the candidates. This route never writes: it
                            exits so you can re-run with the accession you
                            accept, so nothing enters a panel by homology alone.

  --unpublished             Record an allele that has no public accession.

Any route that resolves an accession also requires it to be >=99% identical to
the supplied sequence, so a mistyped accession cannot be attached to a sequence
it does not describe.

Before anything is written the candidate is screened the same way the panels
were built: intact reading frame, sane length, and close enough to the existing
alleles to be the same gene rather than a paralog.
"""
import argparse
import csv
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Seq import Seq

from .misfit import available_panels, reference_root, resolve_panel

STARTS = ("ATG", "GTG", "TTG")
STOPS = ("TAA", "TAG", "TGA")
ORTHOLOG_MIN_ID = 80.0
LEN_TOLERANCE = 0.20
UA = {"User-Agent": "misfit-manage-refs/2.0"}

_AA = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-2,
                            open_gap_score=-3, extend_gap_score=-0.5)


def _get(url, data=None, timeout=180, tries=4):
    """Fetch a URL, retrying transient failures.

    NCBI truncates large responses often enough that a single attempt is not
    reliable: fetching every CDS from a bacterial genome can come back as an
    IncompleteRead partway through.
    """
    last = None
    for attempt in range(tries):
        try:
            return urllib.request.urlopen(
                urllib.request.Request(url, data=data, headers=UA),
                timeout=timeout).read().decode()
        except Exception as exc:
            last = exc
            if attempt < tries - 1:
                wait = 3 * (attempt + 1)
                print(f"  NCBI request failed ({type(exc).__name__}); "
                      f"retrying in {wait}s...", flush=True)
                time.sleep(wait)
    raise SystemExit(f"NCBI request failed after {tries} attempts: "
                     f"{type(last).__name__}: {last}")


def translate(nt):
    return str(Seq(nt[: len(nt) // 3 * 3]).translate(table="Bacterial")).rstrip("*")


def identity(a, b):
    aln = _AA.align(a, b)[0]
    m = sum(1 for x, y in zip(aln[0], aln[1]) if x == y and x != "-")
    n = sum(1 for x, y in zip(aln[0], aln[1]) if x != "-" and y != "-")
    return 100.0 * m / n if n else 0.0


def infer_gene(description):
    """Guess the gene name from a FASTA header.

    Understands '[gene=ompC]', MISFIT's own 'ompC|WP_...|...' form, and a plain
    '>ompC ...' first token. Returns None when nothing looks like a gene name.
    """
    m = re.search(r"\[gene=([^\]]+)\]", description)
    if m:
        return m.group(1).strip()
    first = description.split()[0] if description.split() else ""
    if "|" in first:
        return first.split("|")[0] or None
    if first and re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]{1,15}", first):
        return first
    return None


def read_sequence(path_or_seq):
    """Accept a FASTA file or a bare nucleotide string.

    Returns (sequence, gene_from_header_or_None).
    """
    p = Path(path_or_seq)
    if p.exists():
        recs = list(SeqIO.parse(p, "fasta"))
        if not recs:
            raise SystemExit(f"no FASTA records in {p}")
        if len(recs) > 1:
            raise SystemExit(f"{p} holds {len(recs)} records; supply one allele at a time")
        return str(recs[0].seq).upper().replace("\n", ""), infer_gene(recs[0].description)
    s = re.sub(r"\s+", "", path_or_seq).upper()
    if s and set(s) <= set("ACGTN"):
        return s, None
    raise SystemExit("--sequence is neither a readable FASTA file nor a nucleotide string")


def check_orf(nt):
    """Return a list of problems; empty means the reading frame is sound."""
    bad = []
    if len(nt) % 3:
        bad.append(f"length {len(nt)} is not a multiple of 3")
    if nt[:3] not in STARTS:
        bad.append(f"starts with {nt[:3]}, not a start codon")
    if nt[-3:] not in STOPS:
        bad.append(f"ends with {nt[-3:]}, not a stop codon")
    aa = str(Seq(nt[: len(nt) // 3 * 3]).translate(table="Bacterial"))
    internal = aa.count("*") - (1 if aa.endswith("*") else 0)
    if internal:
        bad.append(f"{internal} internal stop codon(s)")
    return bad


def resolve_from_assembly(accession, gene, seq, cache="ncbi_cache"):
    """Pull the gene's CDS out of a RefSeq assembly and build the header from it."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from urllib.request import urlopen
    import io as _io, zipfile as _zip
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
           f"{accession}/download?include_annotation_type=CDS_FASTA")
    blob = urlopen(urllib.request.Request(url, headers=UA), timeout=300).read()
    with _zip.ZipFile(_io.BytesIO(blob)) as z:
        name = next(n for n in z.namelist() if n.endswith("cds_from_genomic.fna"))
        raw = z.read(name).decode()
    best = None
    header, chunk = None, []
    for line in raw.splitlines() + [">"]:
        if line.startswith(">"):
            if header:
                tags = dict(re.findall(r"\[(\w+)=([^\]]*)\]", header))
                s = "".join(chunk).upper()
                if tags.get("gene", "").lower() == gene.lower() and s:
                    nuc = header.split("|")[1].split("_cds_")[0] if "|" in header else ""
                    cand = {"protein_id": tags.get("protein_id", "NA"), "nuccore": nuc,
                            "locus_tag": tags.get("locus_tag", "NA"),
                            "product": tags.get("protein", ""), "seq": s}
                    score = identity(translate(seq), translate(s))
                    if best is None or score > best[0]:
                        best = (score, cand)
            header, chunk = line[1:], []
        else:
            chunk.append(line.strip())
    if best is None:
        raise SystemExit(f"{accession} has no CDS annotated as gene={gene!r}")
    return best[1], best[0]


def resolve_from_nuccore(accession, gene, seq):
    """Resolve the header from a nucleotide accession plus the gene name.

    A CDS has no accession of its own, so this pulls the record's annotated CDS
    set and picks the one named for the gene -- giving protein_id, locus_tag and
    the nuccore accession together.
    """
    u = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
         + urllib.parse.urlencode({"db": "nuccore", "id": accession,
                                   "rettype": "fasta_cds_na", "retmode": "text"}))
    txt = _get(u, timeout=300)
    best, header, chunk = None, None, []
    for line in txt.splitlines() + [">"]:
        if line.startswith(">"):
            if header:
                tags = dict(re.findall(r"\[(\w+)=([^\]]*)\]", header))
                s = "".join(chunk).upper()
                if tags.get("gene", "").lower() == gene.lower() and s:
                    nuc = header.split("|")[1].split("_cds_")[0] if "|" in header else accession
                    cand = {"protein_id": tags.get("protein_id", "NA"), "nuccore": nuc,
                            "locus_tag": tags.get("locus_tag", "NA"),
                            "product": tags.get("protein", ""), "seq": s}
                    score = identity(translate(seq), translate(s))
                    if best is None or score > best[0]:
                        best = (score, cand)
            header, chunk = line[1:], []
        else:
            chunk.append(line.strip())
    if best is None:
        raise SystemExit(f"{accession} has no CDS annotated as gene={gene!r}")
    return best[1], best[0]


def verify_protein_accession(acc, seq):
    """Confirm the accession exists and matches the translated sequence."""
    u = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
         + urllib.parse.urlencode({"db": "protein", "id": acc,
                                   "rettype": "fasta", "retmode": "text"}))
    fa = _get(u, timeout=60)
    lines = fa.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise SystemExit(f"NCBI returned no protein record for {acc!r}")
    prot = "".join(l.strip() for l in lines[1:])
    return lines[0][1:], identity(translate(seq), prot)


BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def blast_identify(seq, timeout_s=1800, entrez_query="Enterobacterales[Organism]",
                   rid=None):
    """blastp the translation against RefSeq protein and return the top hits.

    Anything that is not a genuine "no significant similarity" result raises
    with an actionable message, so a slow queue or a network problem is never
    reported back as an unrecognised allele.
    """
    prot = translate(seq)
    if rid:
        print(f"  collecting existing BLAST job {rid}...", flush=True)
    else:
        fields = {"CMD": "Put", "PROGRAM": "blastp", "DATABASE": "refseq_protein",
                  "QUERY": prot, "HITLIST_SIZE": "5"}
        if entrez_query:
            fields["ENTREZ_QUERY"] = entrez_query
        resp = _get(BLAST_URL, data=urllib.parse.urlencode(fields).encode())
        rid = next((l.split("= ")[1].strip() for l in resp.splitlines()
                    if l.strip().startswith("RID =")), None)
        if not rid:
            raise SystemExit("  NCBI did not accept the BLAST submission. Try again, or "
                             "pass --protein-accession / --from-nuccore / --from-assembly.")
        rtoe = next((int(m.group(1)) for m in
                     (re.match(r"\s*RTOE\s*=\s*(\d+)", l) for l in resp.splitlines()) if m), 30)
        print(f"  BLAST submitted (RID {rid}); NCBI estimates ~{rtoe}s, waiting up to "
              f"{timeout_s // 60} min. The public queue is often much slower than the "
              f"estimate.", flush=True)
        time.sleep(min(rtoe, 20))

    link = f"{BLAST_URL}?CMD=Get&RID={rid}"
    t0, ready = time.time(), False
    while time.time() - t0 < timeout_s:
        info = _get(BLAST_URL + "?" + urllib.parse.urlencode(
            {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}))
        if "Status=READY" in info:
            ready = True
            break
        if "Status=FAILED" in info or "Status=UNKNOWN" in info:
            raise SystemExit(f"  BLAST job {rid} failed at NCBI. Try again, or pass "
                             f"--protein-accession / --from-nuccore.")
        print(f"    ...still running ({int(time.time() - t0)}s elapsed)", flush=True)
        time.sleep(15)
    if not ready:
        raise SystemExit(
            f"  BLAST job {rid} was still running after {timeout_s}s. The job is not "
            f"lost - results will appear at\n    {link}\n"
            f"  Collect it later with  --blast-rid {rid}  (no re-submission), give a "
            f"longer --blast-timeout, or pass --protein-accession.")

    # FORMAT_TYPE must be a value NCBI actually implements. "Tabular" is not one:
    # it is silently ignored and the response is an empty page, which is
    # indistinguishable from a search that found nothing.
    out = _get(BLAST_URL + "?" + urllib.parse.urlencode(
        {"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid,
         "ALIGNMENTS": "5", "DESCRIPTIONS": "5"}))
    try:
        root = ET.fromstring(out)
    except ET.ParseError:
        raise SystemExit(f"  could not parse the BLAST result for {rid}. See {link}")

    qlen = int(root.findtext(".//BlastOutput_query-len") or len(prot) or 1)
    hits = []
    for hit in list(root.iter("Hit"))[:5]:
        hsp = hit.find(".//Hsp")
        if hsp is None:
            continue
        ident = int(hsp.findtext("Hsp_identity") or 0)
        alen = int(hsp.findtext("Hsp_align-len") or 0)
        acc = re.search(r"\|([^|]+)\|", hit.findtext("Hit_id") or "")
        hits.append({"accession": acc.group(1) if acc else
                                  (hit.findtext("Hit_accession") or "?"),
                     "title": (hit.findtext("Hit_def") or "").strip(),
                     "identity": 100.0 * ident / alen if alen else 0.0,
                     "coverage": 100.0 * alen / qlen,
                     "aln_len": alen})
    return hits


def screen_against_panel(gene_file, seq, gene):
    """Is this the same gene as the alleles already in the file, and is it sound?

    Returns (problems, notes). Problems block the write; notes are advisory.
    """
    problems, notes = [], []
    problems.extend(check_orf(seq))
    aa = translate(seq)
    existing = list(SeqIO.parse(gene_file, "fasta"))
    if not existing:
        return problems, notes

    canonical = str(existing[0].seq).upper()
    can_aa = translate(canonical)
    best = max(identity(translate(str(r.seq).upper()), aa) for r in existing)
    if best < ORTHOLOG_MIN_ID:
        problems.append(f"only {best:.1f}% identical to the closest allele already in "
                        f"ref_{gene}.fasta (needs >={ORTHOLOG_MIN_ID:.0f}%); "
                        f"this looks like a different gene or a paralog")
    lo, hi = len(can_aa) * (1 - LEN_TOLERANCE), len(can_aa) * (1 + LEN_TOLERANCE)
    if not lo <= len(aa) <= hi:
        problems.append(f"protein is {len(aa)} aa against a canonical {len(can_aa)} aa "
                        f"(outside +/-{LEN_TOLERANCE:.0%}); a short allele can win hit "
                        f"selection and make intact genes look truncated")
    for r in existing:
        if translate(str(r.seq).upper()) == aa:
            problems.append(f"identical protein to {r.id}, already in the panel")
            break
    if best < 95.0:
        notes.append(f"closest existing allele is {best:.1f}% identical")
    return problems, notes


def append_allele(gene_file, header, seq):
    text = gene_file.read_text()
    with open(gene_file, "a") as fh:
        if text and not text.endswith("\n"):
            fh.write("\n")
        fh.write(f">{header}\n{seq}\n")


def append_provenance(root, row):
    path = root / "PROVENANCE.tsv"
    cols = list(csv.DictReader(open(path), delimiter="\t").fieldnames) if path.exists() else None
    if cols is None:
        raise SystemExit(f"{path} not found; cannot record provenance")
    with open(path, "a", newline="") as fh:
        csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                       extrasaction="ignore").writerow(row)


def cmd_add_allele(args, root):
    panel = resolve_panel(args.panel, root) or args.panel
    gene_file = root / panel / f"ref_{args.gene}.fasta"
    if not gene_file.exists():
        have = sorted(p.stem[4:] for p in (root / panel).glob("ref_*.fasta")) \
            if (root / panel).is_dir() else []
        raise SystemExit(f"no gene {args.gene!r} in panel {panel!r}."
                         + (f" Available: {', '.join(have)}" if have else
                            f" Panels: {', '.join(n for n, _ in available_panels(root))}"))
    seq, header_gene = read_sequence(args.sequence)
    if header_gene and header_gene.lower() != args.gene.lower():
        print(f"  note: FASTA header says {header_gene!r}, --gene says {args.gene!r}; "
              f"using --gene")
    print(f"panel {panel} / gene {args.gene}: {len(seq)} nt, {len(translate(seq))} aa\n")

    src, ident_note = identify(
        args, seq, args.gene,
        rerun_hint=(f"misfit-db add-allele --panel {panel} --gene {args.gene} "
                    f"--sequence {args.sequence} --protein-accession {{acc}}"))

    header = (f"{args.gene}|{src['protein_id'] or 'NA'}|"
              f"{src['nuccore'] or 'NA'}|{src['locus_tag'] or 'NA'}")

    # ---- screen before writing
    problems, notes = screen_against_panel(gene_file, seq, args.gene)
    for n in notes:
        print(f"  note: {n}")
    if problems:
        print("\n  REJECTED:")
        for p in problems:
            print(f"    - {p}")
        if not args.force:
            raise SystemExit("\n  Nothing written. Pass --force to override.")
        print("\n  --force given; writing anyway.")

    print(f"\n  header: >{header}")
    if args.dry_run:
        print("  --dry-run: nothing written.")
        return

    append_allele(gene_file, header, seq)
    append_provenance(root, {
        "panel": panel, "gene": args.gene, "ref_id": header,
        "source": "user_added", "assembly_accession": args.from_assembly or "-",
        "nuccore_accession": src.get("nuccore", "-"),
        "protein_accession": src.get("protein_id", "-"),
        "locus_tag": src.get("locus_tag", "-"),
        "organism": args.organism or "-", "product": src.get("product", "-"),
        "location": "-", "evidence": ident_note,
        "length_nt": len(seq), "length_aa": len(translate(seq)),
        "multiple_of_3": len(seq) % 3 == 0,
    })
    print(f"  added to {gene_file}")
    print(f"  provenance row appended to {root/'PROVENANCE.tsv'}")


def identify(args, seq, gene, rerun_hint=None):
    """Resolve the four header fields by whichever route the user supplied.

    Routes are tried in decreasing order of authority: a genome annotation, then
    a protein accession the user vouches for, then - only as a last resort - a
    BLAST search, which is advisory and never writes anything on its own.
    """
    if args.from_assembly:
        rec, score = resolve_from_assembly(args.from_assembly, gene, seq)
        print(f"  {args.from_assembly}: gene={gene} -> {rec['protein_id']} "
              f"locus_tag={rec['locus_tag']} ({score:.1f}% identical to your sequence)")
        if score < 99.0 and not args.force:
            raise SystemExit(f"  only {score:.1f}% identical to that record. "
                             f"Check the assembly/gene, or pass --force.")
        if rec.get("seq") and rec["seq"] != seq:
            print("  note: your sequence differs from the annotated CDS; yours is stored.")
        return rec, f"from assembly {args.from_assembly} annotation"

    if args.from_nuccore:
        rec, score = resolve_from_nuccore(args.from_nuccore, gene, seq)
        print(f"  {args.from_nuccore}: gene={gene} -> {rec['protein_id']} "
              f"locus_tag={rec['locus_tag']} ({score:.1f}% identical to your sequence)")
        if score < 99.0 and not args.force:
            raise SystemExit(f"  only {score:.1f}% identical to that record. "
                             f"Check the accession/gene, or pass --force.")
        if rec.get("seq") and rec["seq"] != seq:
            print("  note: your sequence differs from the annotated CDS; yours is stored.")
        return rec, f"from nuccore {args.from_nuccore} annotation"

    if getattr(args, "unpublished", False):
        print("  --unpublished: recording without a public accession.")
        return ({"protein_id": "NA", "nuccore": args.nuccore or "NA",
                 "locus_tag": args.locus_tag or "NA",
                 "product": args.product or "user-supplied allele"},
                "user-supplied, not from a public record")

    if args.protein_accession:
        title, score = verify_protein_accession(args.protein_accession, seq)
        print(f"  {args.protein_accession}: {title[:66]}")
        print(f"  translation identity to your sequence: {score:.1f}%")
        if score < 99.0 and not args.force:
            raise SystemExit("  that accession does not match your sequence; "
                             "check it or pass --force.")
        return ({"protein_id": args.protein_accession, "nuccore": args.nuccore or "NA",
                 "locus_tag": args.locus_tag or "NA",
                 "product": title.split(" ", 1)[-1]},
                "protein accession verified at NCBI")

    print("  no accession given - asking NCBI what this protein is...")
    hits = blast_identify(seq, timeout_s=getattr(args, "blast_timeout", 1800),
                          rid=getattr(args, "blast_rid", None))
    if not hits:
        raise SystemExit(
            "  BLAST found no similar protein in RefSeq. If the allele really is "
            "novel, record it with --unpublished; otherwise check the sequence.")
    print(f"\n  {'accession':<16} {'ident':>7} {'cov':>7}  product")
    for h in hits:
        print(f"  {h['accession']:<16} {h['identity']:6.1f}% {h['coverage']:6.1f}%  "
              f"{h['title'][:54]}")
    hint = rerun_hint or "--protein-accession {acc}"
    # to stdout, so it stays below the table it refers to
    print("\n  Identification is advisory, so nothing was written. Re-run with the "
          "accession you accept, e.g.:\n    " + hint.format(acc=hits[0]["accession"]),
          flush=True)
    raise SystemExit(1)


def cmd_add_gene(args, root):
    """Create a new gene file in an existing panel."""
    panel = resolve_panel(args.panel, root) or args.panel
    pdir = root / panel
    if not pdir.is_dir():
        raise SystemExit(f"no panel {panel!r}. Panels: "
                         f"{', '.join(n for n, _ in available_panels(root))}")
    gene_file = pdir / f"ref_{args.gene}.fasta"
    if gene_file.exists():
        raise SystemExit(f"{gene_file} already exists - use add-allele to extend it")

    seq, header_gene = read_sequence(args.sequence)
    if header_gene and header_gene.lower() != args.gene.lower():
        print(f"  note: FASTA header says {header_gene!r}, --gene says {args.gene!r}; "
              f"using --gene")
    print(f"panel {panel} / NEW gene {args.gene}: {len(seq)} nt, {len(translate(seq))} aa\n")

    problems = check_orf(seq)
    # A new gene must not duplicate one already in the panel under another name:
    # two files holding the same locus is what put calls on the wrong paralog.
    aa = translate(seq)
    for other in sorted(pdir.glob("ref_*.fasta")):
        for r in SeqIO.parse(other, "fasta"):
            pid = identity(translate(str(r.seq).upper()), aa)
            if pid >= 90.0:
                problems.append(f"{pid:.1f}% identical to {r.id} in {other.name} - "
                                f"the panel would hold this locus twice")
                break
        else:
            continue
        break

    src, ident_note = identify(args, seq, args.gene)
    header = (f"{args.gene}|{src['protein_id'] or 'NA'}|"
              f"{src['nuccore'] or 'NA'}|{src['locus_tag'] or 'NA'}")
    if problems:
        print("\n  REJECTED:")
        for p in problems:
            print(f"    - {p}")
        if not args.force:
            raise SystemExit("\n  Nothing written. Pass --force to override.")
        print("\n  --force given; writing anyway.")
    print(f"\n  header: >{header}")
    if args.dry_run:
        print(f"  --dry-run: would create {gene_file}")
        return
    gene_file.write_text(f">{header}\n{seq}\n")
    append_provenance(root, {
        "panel": panel, "gene": args.gene, "ref_id": header, "source": "user_added",
        "assembly_accession": args.from_assembly or "-",
        "nuccore_accession": src.get("nuccore", "-"),
        "protein_accession": src.get("protein_id", "-"),
        "locus_tag": src.get("locus_tag", "-"), "organism": args.organism or "-",
        "product": src.get("product", "-"), "location": "-", "evidence": ident_note,
        "length_nt": len(seq), "length_aa": len(translate(seq)),
        "multiple_of_3": len(seq) % 3 == 0})
    print(f"  created {gene_file}")
    print(f"  provenance row appended")


def cmd_add_panel(args, root):
    """Create a new organism panel, optionally seeded from a RefSeq assembly."""
    name = args.panel.strip().lower().replace(" ", "_")
    pdir = root / name
    if pdir.exists():
        raise SystemExit(f"panel {name!r} already exists at {pdir}")
    if args.dry_run:
        print(f"  --dry-run: would create {pdir}")
        if args.from_assembly:
            print(f"  seeded from {args.from_assembly} with genes: {args.genes or '(none given)'}")
        return
    pdir.mkdir(parents=True)
    aliases = [a.strip() for a in (args.aliases or "").split(",") if a.strip()]
    if aliases:
        (pdir / "aliases.txt").write_text(
            "# species names or regexes routed to this panel, one per line\n"
            + "\n".join(aliases) + "\n")
        print(f"  routing: {', '.join(aliases)} -> {name}")
    print(f"  created {pdir}")

    if args.from_assembly and args.genes:
        wanted = [g.strip() for g in args.genes.split(",") if g.strip()]
        added, missing = [], []
        for gene in wanted:
            try:
                rec, _ = resolve_from_assembly(args.from_assembly, gene, "ATG" * 10)
            except SystemExit:
                missing.append(gene)
                continue
            header = (f"{gene}|{rec['protein_id'] or 'NA'}|"
                      f"{rec['nuccore'] or 'NA'}|{rec['locus_tag'] or 'NA'}")
            (pdir / f"ref_{gene}.fasta").write_text(f">{header}\n{rec['seq']}\n")
            append_provenance(root, {
                "panel": name, "gene": gene, "ref_id": header, "source": "user_added",
                "assembly_accession": args.from_assembly,
                "nuccore_accession": rec["nuccore"], "protein_accession": rec["protein_id"],
                "locus_tag": rec["locus_tag"], "organism": args.organism or "-",
                "product": rec.get("product", "-"), "location": "-",
                "evidence": f"seeded from assembly {args.from_assembly} annotation",
                "length_nt": len(rec["seq"]), "length_aa": len(translate(rec["seq"])),
                "multiple_of_3": len(rec["seq"]) % 3 == 0})
            added.append(gene)
        print(f"  seeded {len(added)} gene(s) from {args.from_assembly}")
        if missing:
            print(f"  not annotated in that assembly: {', '.join(missing)}")
    else:
        print("  empty panel; add genes with `misfit-db add-gene`")


def main():
    ap = argparse.ArgumentParser(
        description="Add alleles, genes and panels to the MISFIT reference set.")
    ap.add_argument("--ref-root", default=None, help="Reference root (default: installed panels)")
    sub = ap.add_subparsers(dest="command", required=True)

    a = sub.add_parser("add-allele", help="Add one allele to an existing gene in a panel")
    a.add_argument("--panel", required=True, help="Panel or species, e.g. e_coli / 'E coli'")
    a.add_argument("--gene", required=True)
    a.add_argument("--sequence", required=True, help="FASTA file or a nucleotide string")
    a.add_argument("--from-assembly", default=None,
                   help="RefSeq assembly accession, e.g. GCF_000005845.2 (best route)")
    a.add_argument("--from-nuccore", default=None,
                   help="Nucleotide accession, e.g. CP049085.2 or NC_000913.3. With --gene "
                        "this resolves protein_id, locus_tag and nuccore together.")
    a.add_argument("--unpublished", action="store_true",
                   help="The allele is not in any public record (e.g. from your own "
                        "assembly). Records it with NA accessions, flagged in provenance.")
    a.add_argument("--product", default=None, help="Product description for --unpublished")
    a.add_argument("--protein-accession", default=None, help="e.g. WP_000865538.1")
    a.add_argument("--nuccore", default=None)
    a.add_argument("--locus-tag", default=None)
    a.add_argument("--organism", default=None)
    a.add_argument("--blast-rid", default=None,
                    help="Collect an earlier BLAST job by its RID instead of "
                         "submitting a new one")
    a.add_argument("--blast-timeout", type=int, default=1800,
                    help="Seconds to wait for a BLAST identification (default: 1800). "
                         "The public queue routinely takes 10-15 min.")
    a.add_argument("--dry-run", action="store_true", help="Show the header, write nothing")
    a.add_argument("--force", action="store_true", help="Write despite screening problems")

    g = sub.add_parser("add-gene", help="Add a new gene to an existing panel")
    g.add_argument("--panel", required=True)
    g.add_argument("--gene", required=True)
    g.add_argument("--sequence", required=True, help="FASTA file or a nucleotide string")
    g.add_argument("--from-assembly", default=None)
    g.add_argument("--from-nuccore", default=None)
    g.add_argument("--protein-accession", default=None)
    g.add_argument("--unpublished", action="store_true")
    g.add_argument("--nuccore", default=None)
    g.add_argument("--locus-tag", default=None)
    g.add_argument("--product", default=None)
    g.add_argument("--organism", default=None)
    g.add_argument("--blast-rid", default=None,
                    help="Collect an earlier BLAST job by its RID instead of "
                         "submitting a new one")
    g.add_argument("--blast-timeout", type=int, default=1800,
                    help="Seconds to wait for a BLAST identification (default: 1800). "
                         "The public queue routinely takes 10-15 min.")
    g.add_argument("--dry-run", action="store_true")
    g.add_argument("--force", action="store_true")

    n = sub.add_parser("add-panel", help="Create a new organism panel")
    n.add_argument("--panel", required=True, help="Directory name, e.g. serratia")
    n.add_argument("--aliases", default=None,
                   help="Comma-separated species names/regexes routed to this panel, "
                        "e.g. '^serratia'. Written to aliases.txt so --species works.")
    n.add_argument("--from-assembly", default=None,
                   help="Seed the panel from this RefSeq assembly's annotation")
    n.add_argument("--genes", default=None, help="Comma-separated genes to seed")
    n.add_argument("--organism", default=None)
    n.add_argument("--dry-run", action="store_true")

    args = ap.parse_args()
    root = Path(args.ref_root) if args.ref_root else reference_root()
    if args.command == "add-allele":
        cmd_add_allele(args, root)
    elif args.command == "add-gene":
        cmd_add_gene(args, root)
    elif args.command == "add-panel":
        cmd_add_panel(args, root)


if __name__ == "__main__":
    main()
