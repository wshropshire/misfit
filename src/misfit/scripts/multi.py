#!/usr/bin/env python3
"""Run MISFIT over a mixed-species cohort, in parallel.

`misfit` itself handles one panel at a time. This command routes every assembly
to the panel for its own organism and runs them across several processes, for
cohorts of hundreds of assemblies spanning several genera. Routing, panel
resolution and calling all come from the same package functions `misfit` uses,
so there is no second copy of that logic to drift.

    misfit-multi assemblies/*.fasta \\
        -o results.tsv --species-map species_map.tsv --jobs 8

    misfit-multi assemblies/*.fasta \\
        -o results.tsv --species "Klebsiella pneumoniae" --jobs 8

The species map is a TSV or CSV with an assembly column
(fasta_file / assembly / sample_id) and a species column
(species_wgs / species / organism).
"""
import argparse
import logging
import os
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from .misfit import (
    OUTPUT_HEADER,
    detect_mutations,
    discover_ref_genes,
    load_gene_list,
    load_species_map,
    reference_root,
    resolve_panel,
    select_genes,
)


def _available_cpus():
    """Cores this process may actually use.

    os.cpu_count() reports the machine's cores, which on a scheduler-managed
    node is the whole box rather than the allocation - asking for 4 cores and
    spawning 128 workers is a good way to have the job killed. Linux exposes
    the real affinity mask; elsewhere fall back to the machine count.
    """
    try:
        return len(os.sched_getaffinity(0))          # Linux, respects cgroups/taskset
    except AttributeError:
        return os.cpu_count() or 4


def _one(job):
    """Call one assembly against one panel, returning the rows it produced."""
    fasta, ref_dir, preset, genes = job
    tmp = tempfile.NamedTemporaryFile("w+", suffix=".tsv", delete=False)
    tmp.close()
    os.unlink(tmp.name)
    try:
        ref_genes = discover_ref_genes(Path(ref_dir))
        if genes:
            ref_genes = {k: v for k, v in ref_genes.items() if k in genes}
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()):
            detect_mutations(fasta, tmp.name, ref_dir, ref_genes, preset)
        with open(tmp.name) as fh:
            return fh.read().split("\n", 1)[1]        # drop the per-file header
    except Exception as exc:                          # one bad assembly must not stop the run
        name = os.path.basename(fasta).split(".")[0]
        return (f"{name}\tERROR\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t--\t"
                f"Error\t{exc}\t--\t--\n")
    finally:
        if os.path.exists(tmp.name):
            os.unlink(tmp.name)


def main():
    ap = argparse.ArgumentParser(
        description="Run MISFIT over a mixed-species cohort, in parallel.")
    ap.add_argument("assemblies", nargs="+")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--species-map", default=None,
                    help="TSV/CSV mapping each assembly to its species")
    ap.add_argument("--species", default=None,
                    help="Single species for every assembly (ignored with --species-map)")
    ap.add_argument("--ref-root", default=None,
                    help="Reference directory root (default: the installed panels)")
    ap.add_argument("--preset", default="asm10", choices=["asm5", "asm10", "asm20"],
                    help="Divergence tolerance when locating a gene (default: asm10). "
                         "Use asm20 for organisms with no dedicated allele in their panel; "
                         "see `misfit --help` for the full explanation.")
    ap.add_argument("--genes", default=None,
                    help="Comma-separated genes of interest, e.g. 'ftsI,ompC'")
    ap.add_argument("--gene-list", default=None,
                    help="File of genes of interest: one per line, or a table with a gene "
                         "column and an optional panel/species column")
    ap.add_argument("--jobs", type=int, default=_available_cpus(),
                    help="Worker processes (default: the cores available to "
                         "this process, not the machine's total)")
    args = ap.parse_args()

    logging.basicConfig(format="%(asctime)s [%(levelname)s] %(message)s", level=logging.INFO)
    ref_root = Path(args.ref_root) if args.ref_root else reference_root()

    if not args.species_map and not args.species:
        ap.error("give either --species-map or --species")

    wanted = {}
    if args.gene_list:
        wanted = load_gene_list(args.gene_list, ref_root)
    if args.genes:
        wanted.setdefault("*", set()).update(
            g.strip() for g in args.genes.split(",") if g.strip())

    smap = load_species_map(args.species_map) if args.species_map else {}
    jobs, skipped = [], defaultdict(list)
    for fasta in args.assemblies:
        base = os.path.basename(fasta)
        species = (smap.get(base) or smap.get(base.split(".")[0])) if smap else args.species
        if not species:
            skipped["not in the species map"].append(base)
            continue
        panel = resolve_panel(species, ref_root)
        if not panel:
            skipped[f"no panel for {species!r}"].append(base)
            continue
        genes = None
        if wanted:
            selected, _ = select_genes(discover_ref_genes(ref_root / panel), wanted, panel)
            if not selected:
                skipped[f"no requested gene is in panel {panel!r}"].append(base)
                continue
            genes = frozenset(selected)
        jobs.append((fasta, str(ref_root / panel), args.preset, genes))

    by_panel = defaultdict(int)
    for _, ref_dir, _, _ in jobs:
        by_panel[os.path.basename(ref_dir)] += 1
    for panel, n in sorted(by_panel.items()):
        logging.info(f"{panel}: {n} assembly(ies)")
    for reason, files in skipped.items():
        logging.warning(f"skipped {len(files)} assembly(ies), {reason}: "
                        f"{', '.join(files[:5])}{' ...' if len(files) > 5 else ''}")
    if not jobs:
        sys.exit("no assemblies could be routed to a reference panel")

    with open(args.output, "w") as out:
        out.write(OUTPUT_HEADER)
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            for i, rows in enumerate(pool.map(_one, jobs, chunksize=1), 1):
                out.write(rows)
                if i % 50 == 0 or i == len(jobs):
                    logging.info(f"  {i}/{len(jobs)}")
    logging.info(f"wrote {args.output}")


if __name__ == "__main__":
    main()
