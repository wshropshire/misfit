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

  add-ortholog / extend-panel / add-panel --genes-like
                            For genes RefSeq leaves unnamed: take the gene by
                            name if annotated, else by reciprocal-best-hit
                            orthology to an allele in another panel. Evidence
                            (identity, coverage, RBH genome) goes to PROVENANCE.

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


def _download_cds_fasta(accession, cache="ncbi_cache"):
    """The assembly's annotated CDS FASTA (cds_from_genomic.fna), cached on disk.

    Ortholog searches read the same reference genomes over and over (every new
    gene is checked back against its seed panel's genome), so the download is
    kept under `cache/<accession>.cds.fna`.
    """
    import io as _io, zipfile as _zip
    cdir = Path(cache)
    cpath = cdir / f"{accession}.cds.fna"
    if cpath.exists() and cpath.stat().st_size > 0:
        return cpath.read_text()
    url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
           f"{accession}/download?include_annotation_type=CDS_FASTA")
    last = None
    for attempt in range(4):
        try:
            blob = urllib.request.urlopen(urllib.request.Request(url, headers=UA),
                                          timeout=300).read()
            with _zip.ZipFile(_io.BytesIO(blob)) as z:
                name = next(n for n in z.namelist() if n.endswith("cds_from_genomic.fna"))
                raw = z.read(name).decode()
            break
        except StopIteration:
            raise SystemExit(f"{accession}: the download has no CDS annotation "
                             f"(is it a RefSeq GCF_ assembly?)")
        except Exception as exc:
            last = exc
            if attempt < 3:
                print(f"  NCBI download failed ({type(exc).__name__}); retrying...", flush=True)
                time.sleep(3 * (attempt + 1))
    else:
        raise SystemExit(f"could not download {accession}: {type(last).__name__}: {last}")
    cdir.mkdir(parents=True, exist_ok=True)
    cpath.write_text(raw)
    return raw


def assembly_cds(accession, cache="ncbi_cache"):
    """Every annotated CDS of a RefSeq assembly as a list of dicts:
    gene, protein_id, nuccore, locus_tag, product, location, seq (nt), aa."""
    out, header, chunk = [], None, []
    for line in _download_cds_fasta(accession, cache).splitlines() + [">"]:
        if line.startswith(">"):
            if header:
                tags = dict(re.findall(r"\[(\w+)=([^\]]*)\]", header))
                s = "".join(chunk).upper()
                if s and "pseudo" not in tags:
                    nuc = header.split("|")[1].split("_cds_")[0] if "|" in header else ""
                    out.append({"gene": tags.get("gene", ""),
                                "protein_id": tags.get("protein_id", "NA"), "nuccore": nuc,
                                "locus_tag": tags.get("locus_tag", "NA"),
                                "product": tags.get("protein", ""),
                                "location": tags.get("location", "-"),
                                "seq": s, "aa": translate(s)})
            header, chunk = line[1:], []
        else:
            chunk.append(line.strip())
    return out


def resolve_from_assembly(accession, gene, seq, cache="ncbi_cache"):
    """Pull the gene's CDS out of a RefSeq assembly and build the header from it."""
    raw = _download_cds_fasta(accession, cache)
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
    if getattr(args, "genes_like", None):
        return _add_panel_like(args, root, name, pdir)
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


# ---------------------------------------------------------------- orthologs
#
# RefSeq does not name every gene (cirA/fepA in Enterobacter, ramR in some
# Klebsiella, most of the panel in Proteus/Providencia), and the panels were
# built by refusing to guess: those genes sit in UNRESOLVED.tsv. The route below
# fills them without guessing either. The gene is found in the target genome by
# homology to an allele already in another panel, and only accepted when the
# relationship is one-to-one: the target protein's best match back in the seed
# allele's own genome must be the seed allele itself (reciprocal best hit). That
# is what separates cirA from the dozen other TonB-dependent receptors, or ompF
# from ompC/phoE, which a plain best-hit would confuse.

ORTHOLOG_DEFAULT_MIN_ID = 60.0
ORTHOLOG_DEFAULT_MIN_COV = 80.0


def _kmers(aa, k=3):
    return {aa[i:i + k] for i in range(len(aa) - k + 1)}


# Ortholog search uses a conventional protein alignment: BLOSUM62, gap open 11 /
# extend 1, free end gaps. Identity is BLAST-style - matches over every aligned
# column between the first and last aligned residue, internal gaps included. (With
# cheap gaps and identity counted over gap-free columns only, two unrelated
# proteins of equal length score ~75%, which makes any identity threshold useless.)
_ORTH = None


def _orth_aligner():
    global _ORTH
    if _ORTH is None:
        from Bio.Align import substitution_matrices
        _ORTH = Align.PairwiseAligner(mode="global",
                                      substitution_matrix=substitution_matrices.load("BLOSUM62"),
                                      open_gap_score=-11, extend_gap_score=-1)
        _ORTH.end_gap_score = 0
    return _ORTH


def aa_hit(q, t):
    """(identity %, coverage % of q) from a BLOSUM62 global alignment with free end gaps."""
    q, t = q.replace("*", "X"), t.replace("*", "X")
    aln = _orth_aligner().align(q, t)[0]
    cols = list(zip(aln[0], aln[1]))
    idx = [i for i, (x, y) in enumerate(cols) if x != "-" and y != "-"]
    if not idx:
        return 0.0, 0.0
    span = cols[idx[0]: idx[-1] + 1]
    m = sum(1 for x, y in span if x == y and x != "-")
    return 100.0 * m / len(span), 100.0 * len(idx) / len(q)


def best_hits(query_aa, proteome, top=10, prefilter=40, len_tol=0.30):
    """Proteome entries ranked by identity to query_aa (k-mer prefilter, then align).

    Only proteins within +/-len_tol of the query length are considered. A global
    alignment against a much longer protein (a 3,400 aa adhesin against a 660 aa
    receptor) can thread the query through scattered matching residues with cheap
    gaps and report a spuriously high identity, so length is screened first.
    """
    qk = _kmers(query_aa)
    if not qk:
        return []
    lo, hi = len(query_aa) * (1 - len_tol), len(query_aa) * (1 + len_tol)
    pool = [c for c in proteome if lo <= len(c["aa"]) <= hi]
    scored = sorted(pool, key=lambda c: -len(qk & _kmers(c["aa"])))[:prefilter]
    hits = []
    for c in scored:
        pid, cov = aa_hit(query_aa, c["aa"])
        hits.append((pid, cov, c))
    hits.sort(key=lambda h: (-h[0] * min(h[1], 100.0), -h[0]))
    return hits[:top]


def _provenance_rows(root):
    path = root / "PROVENANCE.tsv"
    return list(csv.DictReader(open(path), delimiter="\t")) if path.exists() else []


def seed_allele(root, seed_panel, gene):
    """The canonical (first) allele of `gene` in `seed_panel`, with the assembly it
    came from (from PROVENANCE), so the reciprocal check can use that genome."""
    panel = resolve_panel(seed_panel, root) or seed_panel
    f = root / panel / f"ref_{gene}.fasta"
    if not f.exists():
        return None
    rec = next(SeqIO.parse(f, "fasta"))
    asm = next((r["assembly_accession"] for r in _provenance_rows(root)
                if r["panel"] == panel and r["ref_id"] == rec.id), "-")
    return {"panel": panel, "ref_id": rec.id, "seq": str(rec.seq).upper(),
            "aa": translate(str(rec.seq).upper()), "assembly": asm}


def find_neighbor(gene, target_cds, anchor, product_regex, window=2):
    """Locate `gene` by synteny: a CDS within `window` genes of the anchor gene whose
    RefSeq product matches `product_regex` (e.g. ramR: the TetR-family repressor
    divergently transcribed next to ramA; ampR: the LysR regulator next to ampC).
    The anchor is found by its own allele (>=70% identity, >=80% coverage)."""
    hits = best_hits(anchor["aa"], target_cds, top=1)
    if not hits or hits[0][0] < 70 or hits[0][1] < 80:
        return None, f"anchor {anchor['ref_id'].split('|')[0]} not found in this genome"
    a = hits[0][2]
    idx = next(i for i, c in enumerate(target_cds) if c is a)
    rx = re.compile(product_regex, re.I)
    near = []
    for j in range(max(0, idx - window), min(len(target_cds), idx + window + 1)):
        c = target_cds[j]
        if j != idx and c["nuccore"] == a["nuccore"] and rx.search(c["product"] or ""):
            near.append((abs(j - idx), c))
    if not near:
        return None, (f"no CDS matching /{product_regex}/ within {window} genes of "
                      f"{anchor['ref_id'].split('|')[0]} ({a['locus_tag']})")
    near.sort(key=lambda x: x[0])
    c = near[0][1]
    return c, (f"neighbour of {anchor['ref_id'].split('|')[0]} ({a['locus_tag']}, "
               f"{near[0][0]} gene(s) away) with product '{c['product']}'")


def find_ortholog(gene, target_cds, seed, seed_proteome=None,
                  min_id=ORTHOLOG_DEFAULT_MIN_ID, min_cov=ORTHOLOG_DEFAULT_MIN_COV,
                  anchor=None, product_regex=None):
    """Locate `gene` in a target genome.

    1. RefSeq already names it (gene=...)                  -> evidence refseq_gene_name
    2. otherwise the best homolog of the seed allele, kept only if identity and
       coverage pass and it is the reciprocal best hit of the seed allele in the
       seed genome                                          -> evidence ortholog
    Returns (cds dict, evidence string) or (None, reason).
    """
    named = [c for c in target_cds if c["gene"].lower() == gene.lower()]
    if named:
        if seed:
            named.sort(key=lambda c: -aa_hit(seed["aa"], c["aa"])[0])
        return named[0], "refseq_gene_name"
    if anchor and product_regex:
        cand, ev = find_neighbor(gene, target_cds, anchor, product_regex)
        if cand is not None:
            if seed:     # synteny found it; report how it compares with the seed allele
                pid, cov = aa_hit(seed["aa"], cand["aa"])
                ev += f"; {pid:.1f}% aa identity to {seed['panel']} {gene}"
            return cand, ev
        if not seed:
            return None, ev
    if not seed:
        return None, f"no CDS annotated gene={gene} and no seed allele to search with"
    hits = best_hits(seed["aa"], target_cds)
    if not hits:
        return None, "no homolog found"
    pid, cov, cand = hits[0]
    if pid < min_id or cov < min_cov:
        return None, (f"best homolog {cand['protein_id']} ({cand['product'][:40]}) only "
                      f"{pid:.1f}% identity / {cov:.0f}% coverage to {seed['ref_id']}")
    lo, hi = len(seed["aa"]) * (1 - LEN_TOLERANCE), len(seed["aa"]) * (1 + LEN_TOLERANCE)
    if not lo <= len(cand["aa"]) <= hi:
        return None, (f"best homolog {cand['protein_id']} is {len(cand['aa'])} aa vs "
                      f"{len(seed['aa'])} aa seed (outside +/-{LEN_TOLERANCE:.0%})")
    rbh = "reciprocal check skipped (seed genome unknown)"
    if seed_proteome:
        back = best_hits(cand["aa"], seed_proteome, top=1)
        if not back:
            return None, "reciprocal search returned nothing"
        bpid, _, bc = back[0]
        seed_pid = aa_hit(cand["aa"], seed["aa"])[0]
        if bc["aa"] != seed["aa"] and bpid > seed_pid + 0.5:
            return None, (f"not a reciprocal best hit: {cand['protein_id']} matches "
                          f"{bc['gene'] or bc['locus_tag']} ({bc['product'][:40]}, {bpid:.1f}%) "
                          f"in {seed['assembly']} better than the seed ({seed_pid:.1f}%)")
        rbh = f"reciprocal best hit in {seed['assembly']}"
    return cand, (f"ortholog of {seed['panel']} {seed['ref_id'].split('|')[0]} "
                  f"({seed['ref_id'].split('|')[1]}) by {pid:.1f}% aa identity, "
                  f"{cov:.0f}% coverage; {rbh}")


def append_unresolved(root, panel, gene, assembly, organism, reason):
    path = root / "UNRESOLVED.tsv"
    cols = list(csv.DictReader(open(path), delimiter="\t").fieldnames) if path.exists() else \
        ["panel", "gene", "assembly_accession", "organism", "searched_gene_names",
         "searched_products"]
    rows = [r for r in csv.DictReader(open(path), delimiter="\t")] if path.exists() else []
    if any(r["panel"] == panel and r["gene"] == gene and r["assembly_accession"] == assembly
           for r in rows):
        return
    new = not path.exists()
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow({"panel": panel, "gene": gene, "assembly_accession": assembly,
                    "organism": organism, "searched_gene_names": gene,
                    "searched_products": f"ortholog search: {reason}"})


def reference_assembly_for_taxon(taxon):
    """RefSeq reference genome for a taxon name; falls back to the first complete
    RefSeq genome when the taxon has no designated reference."""
    base = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/"
            + urllib.parse.quote(taxon) + "/dataset_report?")
    import json
    for extra in ({"filters.reference_only": "true"},
                  {"filters.assembly_level": "complete_genome"}):
        q = {"filters.assembly_source": "refseq", "page_size": "5", **extra}
        rep = json.loads(_get(base + urllib.parse.urlencode(q), timeout=120) or "{}")
        for r in rep.get("reports", []):
            acc = r.get("accession", "")
            if acc.startswith("GCF_"):
                org = r.get("organism", {}).get("organism_name", taxon)
                strain = (r.get("organism", {}).get("infraspecific_names", {}) or {}).get("strain", "")
                kind = "reference" if "filters.reference_only" in extra else "complete RefSeq genome"
                return acc, f"{org} {strain}".strip(), kind
    raise SystemExit(f"no RefSeq genome found for taxon {taxon!r}")


# Table 11 start codons. RefSeq annotates real CDS that start with ATT/ATC/ATA/CTG
# (e.g. Providencia ampC/ampD); for an annotated CDS these are accepted.
TABLE11_STARTS = ("ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA")


def check_orf_annotated(nt):
    """check_orf for a RefSeq-annotated CDS: any table-11 start codon is valid."""
    probs = check_orf(nt)
    if nt[:3] in TABLE11_STARTS:
        probs = [p for p in probs if "not a start codon" not in p]
    return probs


def _collision(pdir, aa, gene, min_id=90.0):
    """Is this protein already in the panel under another gene name?

    Uses the BLOSUM62/BLAST-style identity and requires the two proteins to be of
    comparable length; the gap-free identity() scores a short protein (arnE, exbD,
    pmrD) as ~90% 'identical' to unrelated longer ones."""
    for other in sorted(pdir.glob("ref_*.fasta")):
        if other.stem[4:].lower() == gene.lower():
            continue
        for r in SeqIO.parse(other, "fasta"):
            oa = translate(str(r.seq).upper())
            if min(len(oa), len(aa)) / max(len(oa), len(aa)) < 0.8:
                continue
            pid, cov = aa_hit(aa, oa)
            if pid >= min_id and cov >= 80:
                return f"{pid:.1f}% identical to {r.id} in {other.name} - " \
                       f"the panel would hold this locus twice"
    return None


def _write_gene(root, panel, gene, cand, assembly, organism, evidence, dry_run, force,
                allele_min_id=ORTHOLOG_MIN_ID):
    """Create ref_<gene>.fasta, or append to it as a further allele after the usual screen."""
    pdir = root / panel
    gene_file = pdir / f"ref_{gene}.fasta"
    header = (f"{gene}|{cand['protein_id'] or 'NA'}|{cand['nuccore'] or 'NA'}|"
              f"{cand['locus_tag'] or 'NA'}")
    seq = cand["seq"]
    problems = check_orf_annotated(seq)
    if seq[:3] not in ("ATG", "GTG", "TTG"):
        evidence += f"; RefSeq CDS starts with {seq[:3]} (table-11 alternative start)"
    if gene_file.exists():
        existing = [translate(str(r.seq).upper()) for r in SeqIO.parse(gene_file, "fasta")]
        if cand["aa"] in existing:
            print(f"    {panel}/{gene}: identical protein already present - skipped")
            return False
        best = max(aa_hit(cand["aa"], e)[0] for e in existing)
        if best < allele_min_id:
            problems.append(f"only {best:.1f}% identical to the closest allele already in "
                            f"ref_{gene}.fasta (needs >={allele_min_id:.0f}%)")
        can = existing[0]
        if not len(can) * (1 - LEN_TOLERANCE) <= len(cand["aa"]) <= len(can) * (1 + LEN_TOLERANCE):
            problems.append(f"protein is {len(cand['aa'])} aa against a canonical {len(can)} aa")
    else:
        col = _collision(pdir, cand["aa"], gene)
        if col:
            problems.append(col)
    if problems and not force:
        print(f"    {panel}/{gene}: REJECTED - " + "; ".join(problems))
        return False
    print(f"    {panel}/{gene}: >{header}  [{evidence}]")
    if dry_run:
        return True
    if gene_file.exists():
        append_allele(gene_file, header, seq)
    else:
        gene_file.write_text(f">{header}\n{seq}\n")
    append_provenance(root, {
        "panel": panel, "gene": gene, "ref_id": header, "source": "misfit_db_ortholog",
        "assembly_accession": assembly, "nuccore_accession": cand["nuccore"] or "-",
        "protein_accession": cand["protein_id"] or "-", "locus_tag": cand["locus_tag"] or "-",
        "organism": organism or "-", "product": cand["product"] or "-",
        "location": cand.get("location", "-"), "evidence": evidence,
        "length_nt": len(seq), "length_aa": len(cand["aa"]),
        "multiple_of_3": len(seq) % 3 == 0})
    return True


def _add_panel_like(args, root, name, pdir):
    """add-panel --genes-like: the new panel gets every gene the template panel(s)
    carry, each taken from the new organism's genome by name or by reciprocal-best-hit
    orthology to the template allele. Genes not found go to UNRESOLVED.tsv."""
    templates = [t.strip() for t in args.genes_like.split(",") if t.strip()]
    if args.reference_taxon:
        asm, organism, kind = reference_assembly_for_taxon(args.reference_taxon)
        print(f"  {args.reference_taxon}: {kind} {asm} ({organism})")
    elif args.from_assembly:
        asm, organism = args.from_assembly, args.organism or "-"
    else:
        raise SystemExit("--genes-like needs --reference-taxon or --from-assembly")
    genes = ([g.strip() for g in args.genes.split(",")] if args.genes else
             sorted({p.stem[4:] for t in templates
                     for p in (root / (resolve_panel(t, root) or t)).glob("ref_*.fasta")}))
    if not args.dry_run:
        pdir.mkdir(parents=True)
        aliases = [a.strip() for a in (args.aliases or "").split(",") if a.strip()]
        if aliases:
            (pdir / "aliases.txt").write_text(
                "# species names or regexes routed to this panel, one per line\n"
                + "\n".join(aliases) + "\n")
            print(f"  routing: {', '.join(aliases)} -> {name}")
    else:
        print(f"  --dry-run: would create {pdir}")
    target = assembly_cds(asm, args.cache)
    added, missing = [], []
    for gene in genes:
        seed, seed_prot = _seed_for(root, gene, templates, args.cache)
        cand, ev = find_ortholog(gene, target, seed, seed_prot, args.min_identity,
                                 args.min_coverage)
        if cand is None:
            print(f"    {name}/{gene}: not found - {ev}")
            missing.append(gene)
            if not args.dry_run:
                append_unresolved(root, name, gene, asm, organism, ev)
            continue
        if args.dry_run:
            print(f"    {name}/{gene}: {cand['protein_id']} [{ev}]")
            added.append(gene)
            continue
        if _write_gene(root, name, gene, cand, asm, organism, ev, False, args.force):
            added.append(gene)
    print(f"\n  {name}: {len(added)} genes added, {len(missing)} not found"
          + (f" ({', '.join(missing)})" if missing else ""))


def _anchor(root, panel, args):
    """The panel's allele of the --neighbor-of gene (for the synteny route)."""
    g = getattr(args, "neighbor_of", None)
    if not g:
        return None
    a = seed_allele(root, panel, g)
    if a is None:
        raise SystemExit(f"--neighbor-of {g}: panel {panel} has no ref_{g}.fasta")
    return a


def _seed_for(root, gene, seed_panels, cache):
    for sp in seed_panels:
        s = seed_allele(root, sp, gene)
        if s:
            prot = assembly_cds(s["assembly"], cache) if s["assembly"].startswith("GCF_") else None
            return s, prot
    return None, None


def panel_genomes(root, panel):
    """RefSeq genomes the panel draws alleles from, most-used (the primary) first."""
    from collections import Counter
    cnt, org = Counter(), {}
    for r in _provenance_rows(root):
        a = r.get("assembly_accession", "-")
        if r["panel"] == panel and a.startswith("GCF_"):
            cnt[a] += 1
            org.setdefault(a, r.get("organism", "-"))
    return [(a, org[a]) for a, _ in cnt.most_common()]


def cmd_add_ortholog(args, root):
    """Add one gene (or, if the panel already has it, one more allele) found by orthology."""
    panel = resolve_panel(args.panel, root) or args.panel
    if not (root / panel).is_dir():
        raise SystemExit(f"no panel {panel!r}")
    if args.from_provenance:
        return _add_ortholog_from_provenance(args, root, panel)
    if args.reference_taxon:
        asm, organism, kind = reference_assembly_for_taxon(args.reference_taxon)
        print(f"  {args.reference_taxon}: {kind} {asm} ({organism})")
    else:
        asm, organism = args.from_assembly, args.organism or "-"
    if not asm:
        raise SystemExit("give --from-assembly GCF_... or --reference-taxon NAME")
    target = assembly_cds(asm, args.cache)
    seeds = [s.strip() for s in (args.seed or "").split(",") if s.strip()]
    ok = 0
    for gene in [g.strip() for g in args.gene.split(",") if g.strip()]:
        seed, seed_prot = _seed_for(root, gene, seeds or [panel], args.cache)
        cand, ev = find_ortholog(gene, target, seed, seed_prot, args.min_identity,
                                 args.min_coverage, _anchor(root, panel, args), args.product_regex)
        if cand is None:
            print(f"    {panel}/{gene}: not added - {ev}")
            if not args.dry_run:
                append_unresolved(root, panel, gene, asm, organism, ev)
            continue
        ok += _write_gene(root, panel, gene, cand, asm, organism, ev, args.dry_run, args.force)
    print(f"  {ok} written{' (dry run)' if args.dry_run else ''}")


def _add_ortholog_from_provenance(args, root, panel):
    """Try the panel's own genomes in order (primary first) until one carries the gene.

    A gene can be missing from the panel's primary genome yet present in the genus
    (e.g. cirA: absent from the E. hormaechei reference, annotated in E. cloacae,
    E. asburiae, E. roggenkampii ...). The first genome that yields it supplies the
    canonical allele; `extend-panel --from-provenance` then adds the others.
    """
    genomes = panel_genomes(root, panel)
    if not genomes:
        raise SystemExit(f"panel {panel!r} has no RefSeq genomes in PROVENANCE.tsv")
    seeds = [s.strip() for s in (args.seed or "").split(",") if s.strip()]
    for gene in [g.strip() for g in args.gene.split(",") if g.strip()]:
        if (root / panel / f"ref_{gene}.fasta").exists():
            print(f"    {panel}/{gene}: already in the panel - use extend-panel for more alleles")
            continue
        seed, seed_prot = _seed_for(root, gene, seeds or [panel], args.cache)
        anchor = _anchor(root, panel, args)
        reasons = []
        for asm, organism in genomes:
            cand, ev = find_ortholog(gene, assembly_cds(asm, args.cache), seed, seed_prot,
                                     args.min_identity, args.min_coverage,
                                     anchor, args.product_regex)
            if cand is None:
                reasons.append(f"{asm} ({organism}): {ev}")
                continue
            if _write_gene(root, panel, gene, cand, asm, organism, ev, args.dry_run, args.force):
                if reasons:
                    print(f"      (not in {len(reasons)} earlier genome(s): "
                          f"{', '.join(r.split(':')[0] for r in reasons)})")
                break
        else:
            print(f"    {panel}/{gene}: not found in any of {len(genomes)} panel genomes")
            if not args.dry_run:
                append_unresolved(root, panel, gene, genomes[0][0], genomes[0][1],
                                  f"absent from all {len(genomes)} panel genomes; "
                                  + reasons[0].split(": ", 1)[-1])


def cmd_extend_panel(args, root):
    """Add one representative allele per species to every gene (or --genes) of a panel.

    Genomes come from --reference-taxon (each taxon's RefSeq reference) and/or
    --from-provenance (every genome the panel already draws alleles from, so a gene
    added later gets the same species representatives as the genes built with it).
    """
    panel = resolve_panel(args.panel, root) or args.panel
    pdir = root / panel
    genes = ([g.strip() for g in args.genes.split(",")] if args.genes else
             sorted(p.stem[4:] for p in pdir.glob("ref_*.fasta")))
    targets = []
    if args.from_provenance:
        seen = {}
        for r in _provenance_rows(root):
            a = r.get("assembly_accession", "-")
            if r["panel"] == panel and a.startswith("GCF_") and a not in seen:
                seen[a] = r.get("organism", "-")
        targets += [(f"{org} (panel genome)", a, org, "panel genome") for a, org in seen.items()]
    for taxon in [t.strip() for t in (args.reference_taxon or "").split(",") if t.strip()]:
        asm, organism, kind = reference_assembly_for_taxon(taxon)
        targets.append((taxon, asm, organism, kind))
    if not targets:
        raise SystemExit("give --reference-taxon and/or --from-provenance")
    for taxon, asm, organism, kind in targets:
        print(f"\n  {taxon}: {kind} {asm} ({organism})")
        target = assembly_cds(asm, args.cache)
        for gene in genes:
            seed, seed_prot = _seed_for(root, gene, [panel], args.cache)
            cand, ev = find_ortholog(gene, target, seed, seed_prot,
                                     args.min_identity, args.min_coverage)
            if cand is None:
                print(f"    {panel}/{gene}: no {taxon} allele - {ev}")
                continue
            ev = "species representative; " + ev
            _write_gene(root, panel, gene, cand, asm, organism, ev, args.dry_run, args.force,
                        allele_min_id=args.min_identity)


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
    n.add_argument("--reference-taxon", default=None,
                   help="Use this taxon's RefSeq reference genome (e.g. 'Serratia marcescens')")
    n.add_argument("--genes-like", default=None,
                   help="Comma-separated template panel(s), e.g. e_coli,enterobacter: the new "
                        "panel gets every gene they carry, found by name or by reciprocal-best-hit "
                        "orthology to the template allele")
    n.add_argument("--min-identity", type=float, default=50.0,
                   help="Ortholog %% aa identity floor for --genes-like (default 50; distant genera)")
    n.add_argument("--min-coverage", type=float, default=ORTHOLOG_DEFAULT_MIN_COV)
    n.add_argument("--cache", default="ncbi_cache")
    n.add_argument("--force", action="store_true")
    n.add_argument("--dry-run", action="store_true")

    o = sub.add_parser("add-ortholog",
                       help="Add gene(s) to a panel from a genome, by RefSeq name or by "
                            "reciprocal-best-hit orthology to another panel's allele")
    o.add_argument("--panel", required=True)
    o.add_argument("--gene", required=True, help="Gene name, or comma-separated list")
    o.add_argument("--from-assembly", default=None, help="Target RefSeq assembly (GCF_...)")
    o.add_argument("--reference-taxon", default=None,
                   help="Alternative to --from-assembly: the taxon's RefSeq reference genome")
    o.add_argument("--from-provenance", action="store_true",
                   help="Try the panel's own genomes (primary first) until one carries the gene")
    o.add_argument("--neighbor-of", default=None,
                   help="Synteny route: the gene sits within 2 genes of this panel gene, e.g. "
                        "ramR next to ramA, ampR next to ampC (use with --product-regex)")
    o.add_argument("--product-regex", default=None,
                   help="RefSeq product the neighbour must match, e.g. 'TetR|RamR'")
    o.add_argument("--seed", default=None,
                   help="Panel(s) holding the seed allele, e.g. e_coli (default: the panel itself)")
    o.add_argument("--min-identity", type=float, default=ORTHOLOG_DEFAULT_MIN_ID)
    o.add_argument("--min-coverage", type=float, default=ORTHOLOG_DEFAULT_MIN_COV)
    o.add_argument("--organism", default=None)
    o.add_argument("--cache", default="ncbi_cache")
    o.add_argument("--dry-run", action="store_true")
    o.add_argument("--force", action="store_true")

    x = sub.add_parser("extend-panel",
                       help="Add a species-representative allele to every gene of a panel")
    x.add_argument("--panel", required=True)
    x.add_argument("--reference-taxon", default=None,
                   help="Comma-separated species, e.g. 'Serratia nevei,Serratia sarumanii'")
    x.add_argument("--from-provenance", action="store_true",
                   help="Use every RefSeq genome the panel already draws alleles from")
    x.add_argument("--genes", default=None, help="Only these genes (default: all in the panel)")
    x.add_argument("--min-identity", type=float, default=ORTHOLOG_MIN_ID,
                   help="Identity to the panel's existing allele required for a species "
                        "representative (default 80; ~70 for diverse genera)")
    x.add_argument("--min-coverage", type=float, default=ORTHOLOG_DEFAULT_MIN_COV)
    x.add_argument("--cache", default="ncbi_cache")
    x.add_argument("--dry-run", action="store_true")
    x.add_argument("--force", action="store_true")

    args = ap.parse_args()
    root = Path(args.ref_root) if args.ref_root else reference_root()
    if args.command == "add-allele":
        cmd_add_allele(args, root)
    elif args.command == "add-gene":
        cmd_add_gene(args, root)
    elif args.command == "add-panel":
        cmd_add_panel(args, root)
    elif args.command == "add-ortholog":
        cmd_add_ortholog(args, root)
    elif args.command == "extend-panel":
        cmd_extend_panel(args, root)


if __name__ == "__main__":
    main()
