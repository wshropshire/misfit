# Mutation Identification in Sequences and Frameshift/ Indel Tracking (MISFIT)
# misfit.py <assemblies> -o <tsv output>

import argparse
import atexit
import logging
import os
import re
import tempfile
from functools import lru_cache
from pathlib import Path

from Bio import SeqIO

from .kleborate_aligner_core import align_query_to_ref, cull_redundant_hits
from .variants import analyze
from .version import __version__

# ---- Species -> reference panel ----
# Only the matched panel's directory is ever read, so an isolate is compared
# against its own organism's alleles and nothing else.
#
# Patterns are matched against the species name, so a full binomial copied out
# of a metadata sheet resolves without anyone editing this table. Klebsiella is
# split deliberately: the pneumoniae complex, the oxytoca complex and
# K. aerogenes each have their own porin repertoire.
PANEL_PATTERNS = [
    (r"^(escherichia|shigella)\b",                                  "e_coli"),
    (r"^klebsiella (pneumoniae|quasipneumoniae|variicola)\b",       "k_pneumoniae"),
    (r"^klebsiella (oxytoca|michiganensis|grimontii|pasteurii|"
     r"huaxiensis|quasivariicola)\b",                               "k_oxytoca"),
    (r"^klebsiella aerogenes\b",                                    "k_aerogenes"),
    (r"^enterobacter\b",                                            "enterobacter"),
    (r"^citrobacter\b",                                             "citrobacter"),
]

SPECIES_ALIASES = {
    "e coli": "e_coli", "e. coli": "e_coli", "escherichia coli": "e_coli",
    "k pneumoniae": "k_pneumoniae", "k. pneumoniae": "k_pneumoniae",
    "k oxytoca": "k_oxytoca", "k. oxytoca": "k_oxytoca",
    "k aerogenes": "k_aerogenes", "k. aerogenes": "k_aerogenes",
}

# Subdirectories of reference/ that hold retired sequences, not usable panels.
NON_PANEL_DIRS = {"quarantine"}

# ---- Build reference list from a directory ----
FASTA_SUFFIXES = {".fa", ".fasta", ".fna"}

def normalize_species(s: str) -> str:
    if s is None:
        return "e_coli"
    key = re.sub(r"\s+", " ", s.strip().lower())
    return SPECIES_ALIASES.get(key, key.replace(" ", "_"))


def available_panels(ref_root: Path):
    """Every usable reference panel, as (name, gene count)."""
    out = []
    if not ref_root.is_dir():
        return out
    for d in sorted(ref_root.iterdir()):
        if not d.is_dir() or d.name in NON_PANEL_DIRS:
            continue
        n = len([p for p in d.glob("*") if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES])
        if n:
            out.append((d.name, n))
    return out


def panel_aliases(ref_root: Path):
    """Extra species -> panel routing contributed by the panels themselves.

    A panel added with `misfit-db add-panel` cannot appear in PANEL_PATTERNS
    without editing this file, so each panel directory may carry an
    `aliases.txt` of species names or regexes, one per line, that route to it.
    """
    out = []
    for name, _ in available_panels(ref_root):
        f = ref_root / name / "aliases.txt"
        if not f.exists():
            continue
        for line in f.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                out.append((line.lower(), name))
    return out


def resolve_panel(species: str, ref_root: Path):
    """Map a species name to a panel directory name, or None."""
    slug = normalize_species(species)
    names = {n for n, _ in available_panels(ref_root)}
    if slug in names:
        return slug
    key = re.sub(r"\s+", " ", (species or "").strip().lower())
    for pattern, panel in PANEL_PATTERNS:
        if re.search(pattern, key) and panel in names:
            return panel
    for pattern, panel in panel_aliases(ref_root):
        try:
            if re.search(pattern, key) and panel in names:
                return panel
        except re.error:
            if pattern in key and panel in names:
                return panel
    return None

# Column names accepted in a species map, in preference order.
_MAP_KEY_COLS = ("fasta_file", "assembly", "file", "filename", "sample_id")
_MAP_SPECIES_COLS = ("species_wgs", "species", "organism")


def load_species_map(path):
    """Read a TSV/CSV mapping each assembly to its species.

    Keyed on both the file name and the name with extensions stripped, so a
    sheet listing 'ISOLATE.fasta' still matches an argument of 'ISOLATE.fna'.
    """
    import csv
    delim = "\t" if str(path).endswith((".tsv", ".tab")) else None
    with open(path, newline="") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        if delim is None:
            delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        rows = list(csv.DictReader(fh, delimiter=delim))
    if not rows:
        raise SystemExit(f"species map is empty: {path}")
    cols = rows[0].keys()
    key_col = next((c for c in _MAP_KEY_COLS if c in cols), None)
    sp_col = next((c for c in _MAP_SPECIES_COLS if c in cols), None)
    if not key_col or not sp_col:
        raise SystemExit(
            f"species map {path} needs an assembly column ({'/'.join(_MAP_KEY_COLS)}) "
            f"and a species column ({'/'.join(_MAP_SPECIES_COLS)}); found: {list(cols)}")
    out = {}
    for r in rows:
        name, species = (r[key_col] or "").strip(), (r[sp_col] or "").strip()
        if not name or not species:
            continue
        out[name] = species
        out[name.split(".")[0]] = species
    return out


def reference_root() -> Path:
    """Where the packaged reference panels live.

    Derived from this module rather than the `misfit` package, which has no
    __init__.py and so is a namespace package whose __file__ is None.
    """
    return Path(__file__).resolve().parent.parent / "reference"


_GENE_COLS = ("gene", "gene_name", "target")
_GENE_PANEL_COLS = ("panel", "species", "species_wgs", "organism", "reference")


def load_gene_list(path, ref_root=None):
    """Read genes of interest from a file.

    One gene per line:

        ftsI
        ompC

    or a table with a gene column and an optional panel/species column, so a
    different subset can apply to each organism:

        gene     species
        ftsI     Escherichia coli
        ompK36   Klebsiella pneumoniae

    Returns {panel_or_'*': {gene, ...}}. Genes listed without a panel go under
    '*' and apply to every panel.
    """
    import csv as _csv
    lines = [ln for ln in open(path).read().splitlines()
             if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        raise SystemExit(f"gene list is empty: {path}")

    delim = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
    wanted = {}
    if delim:
        rows = list(_csv.DictReader(lines, delimiter=delim))
        cols = list(rows[0].keys()) if rows else []
        gcol = next((c for c in _GENE_COLS if c in cols), None)
        pcol = next((c for c in _GENE_PANEL_COLS if c in cols), None)
        if not gcol:
            raise SystemExit(f"gene list {path} needs a gene column "
                             f"({'/'.join(_GENE_COLS)}); found: {cols}")
        for r in rows:
            gene = (r.get(gcol) or "").strip()
            if not gene:
                continue
            key = "*"
            if pcol and (r.get(pcol) or "").strip():
                raw = r[pcol].strip()
                key = (resolve_panel(raw, ref_root) if ref_root else None) or raw
            wanted.setdefault(key, set()).add(gene)
    else:
        wanted["*"] = {ln.strip() for ln in lines}
    return wanted


def select_genes(ref_genes, wanted, panel=None):
    """Narrow a panel's genes to those asked for, matched case-insensitively.

    Returns (selected, unmatched). A gene asked for under a different panel is
    simply not selected here, and is not reported as unmatched.
    """
    if not wanted:
        return ref_genes, []
    names = set(wanted.get("*", set()))
    if panel and panel in wanted:
        names |= wanted[panel]
    if not names:
        return {}, []
    by_label = {}
    for key in ref_genes:
        label = key[4:] if key.startswith("ref_") else key
        by_label[label.lower()] = key
    selected, unmatched = {}, []
    for name in sorted(names):
        key = by_label.get(name.lower())
        if key:
            selected[key] = ref_genes[key]
        else:
            unmatched.append(name)
    return selected, unmatched


def apply_gene_filter(ref_genes, wanted, panel, ref_dir):
    """Filter a panel's genes, reporting anything that matched nothing."""
    if not wanted:
        return ref_genes
    selected, unmatched = select_genes(ref_genes, wanted, panel)
    if unmatched:
        every = sorted((k[4:] if k.startswith("ref_") else k) for k in ref_genes)
        logging.warning(f"{panel or ref_dir}: no such gene(s): {', '.join(unmatched)}")
        logging.warning(f"  available: {', '.join(every)}")
    if not selected:
        raise SystemExit(
            f"none of the requested genes exist in panel {panel or ref_dir}. "
            f"Run with --list-organisms to see the panels, or drop --genes/--gene-list.")
    return selected


def discover_ref_genes(ref_dir: Path):
    """
    Return a dict {gene_name: fasta_basename} discovered under ref_dir.
    """
    if not ref_dir.exists() or not ref_dir.is_dir():
        raise FileNotFoundError(f"Reference directory not found: {ref_dir}")

    ref_genes = {}
    for p in sorted(ref_dir.glob("*")):
        if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES:
            gene_name = p.stem              # filename without extension
            ref_genes[gene_name] = p.name   # basename only (not full path)
    if not ref_genes:
        raise RuntimeError(f"No FASTA files found in: {ref_dir}")
    return ref_genes


@lru_cache(maxsize=None)
def load_reference(ref_path: str):
    """Parse a gene's reference file once per process.

    Returns (canonical_id, canonical_nt, expected_aa_len). Panels do not change
    during a run, and this file was previously re-read twice for every gene of
    every assembly.

    The canonical allele is the first record: always the panel's model reference
    genome. Calls are numbered against whichever allele matched, and this one
    supplies the equivalent canonical numbering reported in Notes.

    expected_aa_len is the longest translation in the file, used as the yardstick
    for deciding a protein is truncated.
    """
    canonical_id, canonical, max_aa = None, "", 0
    for record in SeqIO.parse(ref_path, "fasta"):
        if canonical_id is None:
            canonical_id, canonical = record.id, str(record.seq).upper()
        try:
            aa = str(record.seq.translate(table="Bacterial", to_stop=False)).rstrip("*")
            max_aa = max(max_aa, len(aa))
        except Exception:
            continue
    return canonical_id, canonical, max_aa


@lru_cache(maxsize=None)
def panel_query(ref_dir: str, ref_files: tuple):
    """One FASTA holding every allele in the panel, plus an id -> gene index.

    minimap2 indexes its target, so calling it once per gene re-indexed the whole
    assembly for every gene in the panel -- 54 times over, and ~48x slower than
    aligning the panel in a single pass. The combined query returns exactly the
    same alignments; they are simply partitioned by gene afterwards.
    """
    tmp = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
    gene_of = {}
    for gene, fname in ref_files:
        label = gene[4:] if gene.startswith("ref_") else gene
        for record in SeqIO.parse(os.path.join(ref_dir, fname), "fasta"):
            gene_of[record.id] = label
            tmp.write(f">{record.id}\n{str(record.seq).upper()}\n")
    tmp.close()
    atexit.register(lambda path=tmp.name: os.path.exists(path) and os.unlink(path))
    return tmp.name, gene_of


MUTATION_LABELS = {
    "frameshift": "Frameshift",
    "premature_stop": "Premature stop",
    "translation_error": "Translation error",
    "inframe_indel": "In-frame indel",
    "missense": "Missense",
    "inframe_indel_missense": "In-frame indel and Missense",
    "truncation": "Truncation",
    "intact": "WT",
}

# The same event under the notation other tools use, so results can be checked
# against their output. Every mapping below was verified by running that tool on
# this cohort, not inferred: AMRFinderPlus agrees on all 246 Escherichia ftsI
# calls, and Kleborate reports OmpK36GD / OmpK36TD for the isolates MISFIT calls
# p.Gly134_Asp135dup / p.Asp135_Thr136dup.
# Keyed by panel as well as gene: AMRFinderPlus ftsI point mutations are
# E. coli determinants and Kleborate's Omp designations are Klebsiella ones, so
# an identical substitution in another genus must not inherit the label.
EXTERNAL_NOTATION = {
    ("e_coli", "ftsI", "p.Tyr334_Asn337dup"): "AMRFinderPlus ftsI_N337NYRIN",
    ("e_coli", "ftsI", "p.Ile336_Asn337insLysTyrArgIle"): "AMRFinderPlus ftsI_I336IKYRI",
    ("e_coli", "ftsI", "p.Ala498Thr"): "AMRFinderPlus ftsI_A498T",
    ("k_pneumoniae", "ompK36", "p.Gly134_Asp135dup"): "Kleborate OmpK36GD",
    ("k_pneumoniae", "ompK36", "p.Asp135_Thr136dup"): "Kleborate OmpK36TD",
    ("k_oxytoca", "ompK36", "p.Gly134_Asp135dup"): "Kleborate OmpK36GD",
    ("k_oxytoca", "ompK36", "p.Asp135_Thr136dup"): "Kleborate OmpK36TD",
    ("k_aerogenes", "ompK36", "p.Gly134_Asp135dup"): "Kleborate OmpK36GD",
    ("k_aerogenes", "ompK36", "p.Asp135_Thr136dup"): "Kleborate OmpK36TD",
}


def external_notation(panel, gene, variants):
    """Cross-reference strings for any variant another tool names differently."""
    seen = []
    for v in variants:
        alias = EXTERNAL_NOTATION.get((panel, gene, v.strip()))
        if alias and alias not in seen:
            seen.append(alias)
    return seen


PRESET_HELP = """\
alignment presets (--preset)
  A gene is located by aligning the panel's alleles to the assembly. The preset
  sets how far the isolate's copy may have drifted from the nearest allele and
  still be found. Too strict and the gene reads as "Not found", which looks like
  absence; too loose and alignment gets noisier at the gene's edges.

  asm10   DEFAULT. Tolerates roughly 10% divergence. Correct whenever the panel
          holds an allele close to the isolate -- which is the normal case, since
          panels carry per-species and per-phylogroup alleles.
            misfit isolate.fasta -o out.tsv --species "Klebsiella pneumoniae"

  asm20   Tolerates roughly 20%. Use for an organism with no dedicated allele in
          its panel: a congener routed to the genus panel (say Citrobacter
          sedlakii against the C. freundii-based panel), or any species listed in
          reference/UNRESOLVED.tsv. On divergent input this recovers genes asm10
          misses outright, and it does not make paralog mis-assignment more
          likely -- the extra cost is occasional alignment noise at gene ends.
            misfit isolate.fasta -o out.tsv --species "Citrobacter sedlakii" --preset asm20

  asm5    Tolerates roughly 5%. Only for confirming a near-exact match to a known
          allele. Too strict for general use: ordinary within-species variation
          is enough to lose a gene entirely.
            misfit isolate.fasta -o out.tsv --species "E coli" --preset asm5

  Unsure? Run the default, then re-run the assemblies whose genes came back
  "Not found" with --preset asm20 and see whether they are recovered.
"""

OUTPUT_HEADER = ("Assembly\tGene\tGene_ID\tContig\tContig_Start\tContig_End\tStrand\t"
                 "Nuc_Len\tNuc_ID\tNuc_COV\tAA_Len\tAA_ID\tAA_Cov\t"
                 "Mutation_Type\tMutation_Description\tHGVS_c\tNotes\n")


def detect_mutations(assembly_file, output_file, ref_dir, ref_genes, preset="asm10"):
    output_lines = []
    asm_id = os.path.basename(assembly_file).split(".")[0]

    # minimap2 runs with the reference gene panel as its query and the assembly
    # as its target, so hit.query_* describes the reference allele and hit.ref_*
    # describes the sample's copy of the gene. The whole panel goes in one pass;
    # hits are then split by which gene each allele belongs to.
    combined, gene_of = panel_query(ref_dir, tuple(sorted(ref_genes.items())))
    by_gene = {}
    for hit in align_query_to_ref(combined, assembly_file, preset=preset):
        by_gene.setdefault(gene_of.get(hit.query_name), []).append(hit)

    for gene, ref_file in ref_genes.items():
        ref_path = os.path.join(ref_dir, ref_file)
        label = gene[4:] if gene.startswith("ref_") else gene

        hits = by_gene.get(label, [])
        if not hits:
            output_lines.append(f"{asm_id}\t{label}\t--\t--\t--\t--\t--\t--\t--\t--\t--\t"
                                f"--\t--\tNot found\t--\t--\t--")
            continue

        best_hit = cull_redundant_hits(hits)
        ref_nucl = best_hit.query_seq.upper()      # reference allele, aligned span
        sample_nucl = best_hit.ref_seq.upper()     # sample copy of the gene

        canonical_id, canonical, expected_ref_aa_len = load_reference(ref_path)

        notes = []
        offset = best_hit.query_start
        # An alignment that starts mid-codon translates the reference out of
        # frame, which invents stop codons in the reference and makes the sample
        # look like it has read through one. Trim both sides to the next codon
        # boundary so the comparison is in frame.
        shift = (-offset) % 3
        if shift:
            ref_nucl = ref_nucl[shift:]
            sample_nucl = sample_nucl[shift:]
            offset += shift
            notes.append(f"alignment began mid-codon at reference nt {best_hit.query_start}; "
                         f"trimmed {shift} base(s) to restore frame")
        if best_hit.query_cov < 100.0:
            notes.append(f"only {best_hit.query_cov:.1f}% of the reference CDS aligned; "
                         "positions outside that span were not examined")

        # The allele minimap2 actually aligned. Re-deriving it by string
        # similarity was both redundant and unreliable: difflib's autojunk
        # heuristic treats every base as junk on DNA over 200 bp and returns 0,
        # which is what produced "NA" here.
        gene_id = best_hit.query_name
        result = analyze(ref_nucl, sample_nucl,
                         expected_ref_aa_len=expected_ref_aa_len,
                         ref_nt_offset=offset,
                         canonical_nt=canonical)
        notes.extend(result.notes)
        # Calls are made against the allele that fits the isolate best, so their
        # numbering follows that allele. Where the canonical allele numbers a
        # variant differently, give that form too -- it is the one that matches
        # the literature and AMRFinderPlus.
        differing = [c for v, c in zip(result.variants, result.canonical) if v != c]
        if differing:
            notes.append("canonical: " + "; ".join(differing))
        notes.extend(external_notation(os.path.basename(str(ref_dir).rstrip("/")),
                                       label, result.variants))

        print(f"→ Final best hit for {label}: {gene_id}")

        aa_len = result.query_aa_len
        aa_cov = (aa_len / expected_ref_aa_len * 100) if expected_ref_aa_len else 0
        msg = MUTATION_LABELS.get(result.kind, "WT")

        line = (
            f"{asm_id}\t{label}\t{gene_id}\t"
            f"{best_hit.ref_name}\t{best_hit.ref_start + 1}\t{best_hit.ref_end}\t{best_hit.strand}\t"
            f"{len(sample_nucl)}\t"
            f"{best_hit.percent_identity:.1f}%\t{best_hit.query_cov:.1f}%\t"
            f"{aa_len}\t{result.aa_identity:.1f}%\t{aa_cov:.1f}%\t{msg}\t"
            f"{result.description}\t{result.hgvs_c}\t"
            f"{'; '.join(notes) if notes else '--'}"
        )
        output_lines.append(line)

    write_mode = 'a' if os.path.exists(output_file) else 'w'
    with open(output_file, write_mode) as f:
        if write_mode == 'w':
            f.write(OUTPUT_HEADER)
        f.write("\n".join(output_lines) + "\n")

def run_with_species_map(args, ref_root, wanted=None):
    """Process a mixed-species set, routing each assembly to its own panel."""
    smap = load_species_map(args.species_map)
    groups, unmapped, no_panel = {}, [], {}
    for fasta in args.assemblies:
        base = os.path.basename(fasta)
        species = smap.get(base) or smap.get(base.split(".")[0])
        if not species:
            unmapped.append(base)
            continue
        panel = resolve_panel(species, ref_root)
        if not panel:
            no_panel.setdefault(species, []).append(base)
            continue
        groups.setdefault(panel, []).append(fasta)

    for panel, files in sorted(groups.items()):
        logging.info(f"{panel}: {len(files)} assembly(ies)")
    if no_panel:
        for species, files in sorted(no_panel.items()):
            logging.warning(f"no reference panel for {species!r}: skipping {len(files)} assembly(ies)")
    if unmapped:
        logging.warning(f"{len(unmapped)} assembly(ies) absent from the species map, skipped: "
                        f"{', '.join(unmapped[:5])}{' ...' if len(unmapped) > 5 else ''}")
    if not groups:
        raise SystemExit("no assemblies could be routed to a reference panel")

    ran = False
    for panel, files in sorted(groups.items()):
        ref_dir = ref_root / panel
        ref_genes = discover_ref_genes(ref_dir)
        if wanted:
            selected, unmatched = select_genes(ref_genes, wanted, panel)
            if not selected:
                logging.warning(f"{panel}: none of the requested genes are in this panel; "
                                f"skipping {len(files)} assembly(ies)")
                continue
            if unmatched:
                logging.warning(f"{panel}: no such gene(s): {', '.join(unmatched)}")
            ref_genes = selected
            logging.info(f"{panel}: {len(ref_genes)} gene(s) selected")
        ran = True
        for fasta in files:
            logging.info(f"Processing {os.path.basename(fasta)} against {panel}")
            detect_mutations(fasta, args.output, str(ref_dir), ref_genes, args.preset)
    if not ran:
        raise SystemExit("no panel contained any of the requested genes")
    logging.info("Mutation detection completed successfully.")


def main():
    parser = argparse.ArgumentParser(
        description="Reference-anchored variant calling in chromosomal resistance "
                    "genes of Enterobacterales.",
        epilog=PRESET_HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("assemblies", nargs='*', help="Input assembly FASTA file(s)")
    parser.add_argument("-o", "--output", help="Output TSV summary file")
    parser.add_argument("--list-organisms", action="store_true",
                        help="List the reference panels available and exit")
    parser.add_argument("--species-map", default=None,
                        help="TSV/CSV mapping each assembly to its species, for a "
                             "directory holding more than one organism. Needs an "
                             "assembly column (fasta_file/assembly/sample_id) and a "
                             "species column (species_wgs/species/organism). Takes "
                             "precedence over --species.")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}", help="Show program's version number and exit")
    parser.add_argument("--species", default="E coli",
                        help="Organism whose reference panel to use; only that panel is read. "
                             "Accepts a full species name ('Enterobacter hormaechei') or a "
                             "panel name ('k_pneumoniae'). See --list-organisms.")
    parser.add_argument("--preset", default="asm10", choices=["asm5", "asm10", "asm20"],
                        help="How much sequence divergence to tolerate when locating a gene "
                             "(default: asm10). See the examples below for when to change it.")
    parser.add_argument("--genes", default=None,
                        help="Comma-separated genes of interest, e.g. 'ftsI,ompC,ompF'. "
                             "Only these are examined instead of the whole panel.")
    parser.add_argument("--gene-list", default=None,
                        help="File of genes of interest: one gene per line, or a table with "
                             "a gene column and an optional panel/species column so a "
                             "different subset applies to each organism.")
    parser.add_argument("--ref-dir",default=None,help="Override reference directory (takes precedence over --species). If set, all FASTA files in this directory are used.")

    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.INFO
    )

    # Determine reference directory
    ref_root = reference_root()

    if args.list_organisms:
        panels = available_panels(ref_root)
        print(f"Reference panels available in {ref_root}\n")
        print(f"  {'panel':16s} {'genes':>6s}  recognised species")
        for name, n in panels:
            species = [p for p, t in PANEL_PATTERNS if t == name]
            hint = {
                "e_coli": "Escherichia, Shigella",
                "k_pneumoniae": "K. pneumoniae / quasipneumoniae / variicola",
                "k_oxytoca": "K. oxytoca complex (oxytoca, michiganensis, grimontii, ...)",
                "k_aerogenes": "K. aerogenes",
                "enterobacter": "Enterobacter (any species)",
                "citrobacter": "Citrobacter (any species)",
            }.get(name, "match by panel name")
            print(f"  {name:16s} {n:6d}  {hint}")
        return

    if not args.assemblies or not args.output:
        parser.error("assemblies and -o/--output are required (unless --list-organisms)")

    wanted = {}
    if args.gene_list:
        wanted = load_gene_list(args.gene_list, ref_root)
    if args.genes:
        wanted.setdefault("*", set()).update(
            g.strip() for g in args.genes.split(",") if g.strip())

    if args.species_map:
        run_with_species_map(args, ref_root, wanted)
        return

    if args.ref_dir:
        ref_dir = Path(args.ref_dir)
    else:
        panel = resolve_panel(args.species, ref_root)
        if panel is None:
            names = ", ".join(n for n, _ in available_panels(ref_root))
            raise SystemExit(
                f"No reference panel for species {args.species!r}.\n"
                f"Available panels: {names}\n"
                f"Run with --list-organisms to see which species each one covers.")
        ref_dir = ref_root / panel

    # Auto-discover references (replaces hard-coded REF_GENES)
    REF_GENES = discover_ref_genes(ref_dir)
    if wanted:
        panel_name = ref_dir.name
        REF_GENES = apply_gene_filter(REF_GENES, wanted, panel_name, ref_dir)
        print(f"[MISFIT] {len(REF_GENES)} gene(s) selected", flush=True)

    # (Optional) small log for clarity
    print(f"[MISFIT] Using {len(REF_GENES)} reference genes from: {ref_dir}", flush=True)

    missing_files = []
    for ref_file in REF_GENES.values():
        ref_path = ref_dir / ref_file
        if ref_path.exists():
            logging.info(f"Found reference: {ref_file}")
        else:
            missing_files.append(ref_file)

    if missing_files:
        logging.error(
            "Missing reference FASTA files:\n" + "\n".join(missing_files)
        )
        raise FileNotFoundError(
            f"Missing reference FASTA files in {ref_dir}:\n" + "\n".join(missing_files)
        )

    for fasta in args.assemblies:
        logging.info(f"Processing assembly: {fasta}")
        detect_mutations(fasta, args.output, str(ref_dir), REF_GENES, args.preset)

    logging.info("Mutation detection completed successfully.")

if __name__ == "__main__":
    main()
