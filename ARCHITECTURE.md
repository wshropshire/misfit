# How MISFIT works

A map of the code, the two problems that shape its design, and the rough edges
worth knowing about.

## The one-paragraph version

MISFIT answers "how does this isolate's copy of gene X differ from a reference
copy of gene X, described in HGVS?" It finds the gene by aligning a panel of
reference alleles against the assembly with minimap2, extracts the sample's
copy, aligns reference protein to sample protein, and walks that alignment
emitting variants **in reference coordinates**. 

## Data flow

```
assembly.fasta  +  reference/<panel>/          all alleles, all genes
        |
        v
  panel_query()                 one combined FASTA, cached per process
        |
        v
  align_query_to_ref()          ONE minimap2 pass per assembly, preset asm10
        |                       hits then partitioned by gene
        v
  cull_redundant_hits()         pick one per gene: most matching bases, then identity
        |
        v
  analyze()                     the whole variant-calling decision tree
        |
        v
  one TSV row                   type, HGVS protein, HGVS nucleotide, notes
```

minimap2 indexes its *target*, so calling it once per gene re-indexed the whole
assembly for every gene in the panel. Aligning the panel in a single pass returns
identical alignments and is roughly an order of magnitude faster per assembly, the
saving growing with panel size. Verified by comparing per-gene against single-pass
alignment across every panel, and by confirming whole-run output is unchanged field
for field.

## Modules

### `scripts/misfit.py` — CLI and orchestration

`main()` resolves the reference directory from `--species` or `--ref-dir`,
discovers genes by globbing `ref_*.fasta` (the filename *is* the gene name), and
loops over assemblies. `detect_mutations()` does the per-gene work and writes
the row.

Also holds `EXTERNAL_NOTATION`, the lookup from an HGVS string to the name
another tool gives the same event (`p.Tyr334_Asn337dup` -> AMRFinderPlus
`ftsI_N337NYRIN`). Keyed by `(panel, gene, variant)` so a Klebsiella label
cannot leak onto an *E. coli* call.

### `scripts/kleborate_aligner_core.py` — finding the gene

Adapted from Kleborate. `Alignment` parses a PAF line; `align_query_to_ref()`
shells out to minimap2; `cull_redundant_hits()` picks the winner.

**The naming trap.** minimap2 is invoked with the *assembly* as its target and
the *gene panel* as its query. So on the returned object:

| attribute | what it actually is |
|---|---|
| `hit.query_seq`, `hit.query_start`, `hit.query_cov` | the **reference allele** |
| `hit.ref_seq`, `hit.ref_name`, `hit.ref_start` | the **sample's** contig and copy |

`ref_seq` is reverse-complemented when the strand is `-`, so it always reads in
gene orientation. This inversion is the single most confusing thing in the
codebase and is why `detect_mutations()` immediately renames them to
`ref_nucl` / `sample_nucl`.

### `scripts/variants.py` — the actual science

`analyze()` is the decision tree:

1. Keep the **untrimmed** sequences for frame analysis; trim to whole codons
   only for translation. Rounding a CDS down to a multiple of three is exactly
   what hides a single-base indel.
2. If the lengths differ at all, or the query has an internal stop, run
   `find_frame_disruption()`.
3. **Frameshift** -> walk forward to the first amino acid that actually changes
   (HGVS's anchor, and well defined even when the causal indel sits in a
   homopolymer where the nucleotide position is not). Emit `fsTer<n>`, plus any
   variants upstream of the break, which are still translated faithfully.
4. Otherwise align the proteins and walk -> substitutions, indels, stops.
5. **Premature stop** -> discard variants downstream of it; they are not
   translated in vivo.
6. **Truncation** -> query protein under 90% of the expected reference length.
7. Otherwise classify as in-frame indel, missense, or intact.

Supporting pieces:

- `find_frame_disruption()` scores each indel **run as a whole** — a 12-base
  insertion passes through offsets of 1 and 2 on its way to 12, and judging a
  run mid-way turns an in-frame insertion into a phantom frameshift. It also
  requires the frame never to recover; a divergent reference scatters
  compensating gaps that would otherwise read as a frameshift.
- `_left_normalize()` / `_right_normalize()` shift gap runs 5' or 3'. Frame
  analysis uses 5'-most (earliest possible break); HGVS output uses 3'-most, as
  the standard requires.
- `_ins_or_dup()` writes an insertion that exactly copies the sequence just 5'
  of it as a duplication. This is why YRIN is `p.Tyr334_Asn337dup` and YRIK,
  which differs at one residue, stays `p.Ile336_Asn337insLysTyrArgIle`.
- `hgvs_c()` produces the nucleotide-level description with the same
  normalisation, so the two columns always agree about what happened.

### `scripts/manage_refs.py` — curating the reference set

Backs the `misfit-db` console script (`add-allele`, `add-gene`, `add-panel`). The
interesting part is `identify()`, which resolves a reference header from whatever the
user supplied, trying routes in decreasing order of authority: genome annotation, then
a protein accession, then BLAST. Everything except BLAST is authoritative enough to
write; BLAST only reports candidates and exits, so an accession enters a panel only by
an explicit human choice.

Two things worth knowing:

- **NCBI's `FORMAT_TYPE` has no `Tabular` value.** Passing one is silently ignored and
  the response is an empty page — indistinguishable from a search with no hits, which
  is how a perfectly ordinary *E. coli* `ompC` once came back unidentifiable. The code
  uses `XML` and parses `Hit_id` for the versioned accession.
- **The public BLAST queue is much slower than `RTOE` claims.** A 373 aa OmpC query
  measured 866 s against an `RTOE` estimate of 45 s — a factor of ~19. `blast_identify()`
  waits 30 minutes by default, distinguishes a timeout from an empty result, and prints
  the RID so `--blast-rid` can collect the finished job later rather than re-submitting.
  Do not tune this down on the basis of `RTOE`; it is not predictive.

`screen_against_panel()` is the guard that keeps the panels clean, and enforces from the
outside the same invariant `cull_redundant_hits()` depends on: one gene per file.

### `scripts/upstream.py` — EXPERIMENTAL, opt-in

Detects insertions upstream of a gene using the cohort's modal architecture as
the wild type. Works where that architecture is conserved (it called BIDMC169's
ompK35 IS insertion at 777 bp, 20 bp upstream). Fails where it is not — see
Known limitations. Never runs as part of a normal MISFIT run.

## The two problems that shape the design

**Reference coordinates.** The original code numbered variants by alignment
column index. An insertion shifts every downstream column, so an isolate with
the PBP3 duplication had its A413V reported as A417V. The fix is to walk the
alignment maintaining a separate counter of *reference residues consumed*, and
number from that. Only insertions cause the drift; deletions leave reference
numbering intact.

**One gene per reference file.** RefSeq gives paralogous genes the same name. If
two land in one file, both align to their own locus in the assembly and hit
selection picks whichever yields more matching bases — reporting the wrong
locus — silently placing a gene's calls on a distant paralog. Two defences: `cull_redundant_hits()` scores by **absolute matching
bases** (coverage-as-a-fraction-of-itself lets a short allele win by being
short), and the panels are audited so no file holds two different genes.

## Commands

| command | module | role |
|---|---|---|
| `misfit` | `scripts/misfit.py` | call one assembly (or many) against one panel |
| `misfit-multi` | `scripts/multi.py` | route a mixed-species cohort per organism, in parallel |
| `misfit-wide` | `scripts/wide.py` | reshape long-form results to one row per sample |
| `misfit-db` | `scripts/manage_refs.py` | add alleles, genes and panels with provenance |

`multi.py` and `wide.py` both import the package's own routing helpers
(`load_species_map`, `resolve_panel`, `select_genes`) rather than reimplementing them,
so panel resolution cannot drift between the single-panel and cohort paths.

## Panel-building scripts (in `analysis_plan/misfit_results/scripts/`)

Not packaged, which is a wart:

| script | role |
|---|---|
| `annotate_frequency.py` | splits variants into rare / common / fixed-vs-reference |
| `build_panels.py` + `fetch_refs.py` | build reference panels from RefSeq with provenance |
| `select_alleles.py` | pick screened wild-type alleles per phylogroup |
| `phylotype.py` | ezclermont typing |
| `make_provenance.py` | regenerate PROVENANCE.tsv from panel contents |

## Coordinates and the canonical allele

Panels hold several alleles per gene, and those alleles differ in length. A call is
made and numbered against the allele that matched, which is what keeps ordinary
lineage variation out of the results. But a position in one allele is not the same
position in another: MG1655 ompC 195 is 191 in the phylogroup C allele, which carries
a four-residue deletion at 181-184.

So `Mutation_Description` uses matched-allele numbering, and where the canonical
allele — the first record in the gene's file, always the panel's model reference —
numbers a variant differently, `Notes` carries `canonical: p.Arg195Leu`. About 1.4% of
rows need one, concentrated in the length-variable porins.

Downstream aggregation should key on the canonical form where present.

## Rough edges

1. **The `query`/`ref` inversion** in `Alignment` should be renamed outright.
2. **The panel-building scripts live outside the package** (`annotate_frequency.py`,
   `build_panels.py`, `select_alleles.py`, `phylotype.py`). They belong in it, with the
   panel build as a `misfit-db` subcommand. The runner and reshaper have been moved in;
   these have not.
3. **No CI, and the test suite is a script** rather than pytest. It covers the variant
   caller well (12 ground-truth cases) and nothing else.
4. **Partial hits still produce confident calls.** A gene covered at 26% of the
   reference is reported as a truncation rather than flagged unassessable. Treat
   anything under ~90% `Nuc_COV` with suspicion.

## Known limitations

- **Coding sequence only.** A gene whose promoter carries an IS element is
  reported as an intact ORF, correctly but incompletely.
- **Reference divergence shows up as missing data.** minimap2 `asm10` needs
  ~90% nucleotide identity; below that the gene silently reads "Not found",
  which looks like absence. This costs ~80% of gene calls for Enterobacter and
  Citrobacter species other than the panel's own reference species.
- **`upstream.py` needs a conserved upstream architecture.** Where the cohort
  has none (Enterobacter ompF has 17 distinct proximal classes across 80
  isolates) there is no wild type to measure against.

## Running it

```bash
micromamba activate misfit_env

misfit assembly.fasta -o out.tsv --species "K pneumoniae"
misfit-multi assemblies/*.fasta -o cohort.tsv --species-map map.tsv --jobs 8
misfit-wide cohort.tsv -o wide.xlsx --species-map map.tsv
misfit-db add-allele --panel e_coli --gene ompC --sequence new.fasta \
    --protein-accession WP_000865587.1

python tests/test_variants.py
```
