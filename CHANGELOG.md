# Changelog

## v2.0.0

A large update: six organism panels instead of two, HGVS nomenclature throughout,
a corrected coordinate system, and substantially faster runs.

### Fixed — coordinates

**Variant positions were alignment-column indices, not reference positions.** Any
residue downstream of an insertion was renumbered by the length of that insertion.
The practical consequence: an isolate carrying the PBP3 333-loop duplication had
its `A413V` substitution reported as `A417V` — the "A413V aka A417V" ambiguity
familiar from the literature. Positions are now counted in reference residues, so
an insertion no longer shifts what follows it.

The corrected PBP3 calls agree with the independent AMRFinderPlus determinants for
the YRIN, YRIK and A498T alleles.

### Fixed — reference panels

- **A truncated `ompK35` allele** (435 nt / 144 aa) sat alongside the real 1080 nt
  gene. Because hit selection scored coverage as a fraction of each allele's own
  length, the short allele won for most *K. pneumoniae* isolates and made their
  intact porins look frameshifted. Selection now scores by absolute matching bases,
  and the allele is in `reference/quarantine/`.
- **Paralogs sharing a gene name.** RefSeq names several porin and PBP paralogs
  identically; where two landed in one file the wrong locus could win, putting
  calls on a distant paralog rather than the intended gene. Panels are now audited
  so no file holds two different genes; the removed sequences are in
  `reference/quarantine/` with the evidence for each decision.
- **Klebsiella porin naming corrected.** The `ada`/`rcsDB`-syntenic porin is
  `ompK36`, not `ompC`; these genomes carry a separate gene RefSeq calls `ompC`.

### Fixed — calling

- **Mid-codon alignment starts** translated the reference out of frame, inventing
  stop codons in the reference and making ordinary samples look like they had read
  through one. Both sequences are now trimmed to the next codon boundary.
- **`Gene_ID` could be `NA`.** It was re-derived by string similarity; it is now
  taken directly from the alignment, so it always names the allele actually used.
- **Ambiguous bases no longer truncate a gene.** A single `N` used to cut the
  sequence at that point and produce a false truncation.
- **Frameshift detection** now scores each indel run as a whole and requires the
  frame never to recover. Judging runs in isolation turned in-frame insertions into
  phantom frameshifts, which could affect every isolate for a divergent gene.
- Frameshifts anchor to the first changed amino acid, which stays well defined when
  the causal indel sits in a homopolymer and its nucleotide position does not.

### Added — nomenclature

Descriptions follow **HGVS**: three-letter codes, `Ter` for termination, indels
placed 3'-most, and an insertion copying its 5' neighbour written as a duplication.

| event | before | now |
|---|---|---|
| PBP3 333-loop YRIN | `p.P333insYRIN` | `p.Tyr334_Asn337dup` |
| PBP3 333-loop YRIK | `p.P333insYRIK` | `p.Ile336_Asn337insLysTyrArgIle` |
| substitution | `p.A413V` | `p.Ala413Val` |
| frameshift | `p.S90YfsX16` | `p.Ser90TyrfsTer16` |
| stop read-through | `p.*187Q` | `p.Ter187GlnextTer17` |

YRIN is a duplication because it copies the adjacent YRIN exactly; YRIK differs at
one residue and stays an insertion.

New `HGVS_c` column gives the same change in nucleotide coordinates
(`c.1000_1011dup`, `c.268_269del`, `c.1238C>T`). New `Contig`, `Contig_Start`,
`Contig_End` and `Strand` columns locate the gene in the assembly.

`Notes` carries the equivalent notation used by other tools, each verified by
running that tool rather than assumed: `AMRFinderPlus ftsI_N337NYRIN`,
`ftsI_I336IKYRI`, `ftsI_A498T`, `Kleborate OmpK36GD`, `OmpK36TD`. Cross-references
are scoped by panel so a label cannot leak across genera.

Calls are made and numbered against whichever allele fits the isolate best. Where
the panel's canonical allele numbers a variant differently, `Notes` carries the
canonical form (`canonical: p.Arg195Leu`).

### Added — organisms

Six panels, ~900 sequences, every one traceable to a RefSeq record:

| panel | genes | covers |
|---|---|---|
| `e_coli` | 56 | Escherichia, Shigella |
| `k_pneumoniae` | 54 | K. pneumoniae / quasipneumoniae / variicola |
| `k_oxytoca` | 51 | K. oxytoca complex |
| `k_aerogenes` | 54 | K. aerogenes |
| `enterobacter` | 43 | Enterobacter |
| `citrobacter` | 53 | Citrobacter |

Each gene holds a canonical allele plus representatives for the E. coli phylogroups
(assigned by ezclermont, never presumed) and for species within each genus. This
matters because a single reference per genus leaves congeners below minimap2's
identity threshold, so most of their genes return `Not found` — which reads as gene
absence rather than divergence. Per-species alleles recover the great majority of
those genes.

Every candidate allele is screened so a reference can never itself carry a
resistance determinant: intact ORF, expected length, no PBP3 333-loop insertion, no
carbapenemase or ESBL in the source genome, and modal across the screened genomes
of that group. `reference/PROVENANCE.tsv` records the assembly, nucleotide,
locus-tag and protein accession behind every sequence; `reference/UNRESOLVED.tsv`
lists genes deliberately omitted because RefSeq does not identify them.

### Added — interface

- `--list-organisms` lists the available panels.
- `--species` accepts a full binomial (`Enterobacter hormaechei`) or a panel name;
  unknown organisms give a readable error instead of a traceback.
- `--species-map` runs a directory holding several organisms, routing each assembly
  to its own panel.
- `--genes` / `--gene-list` narrow a run to genes of interest instead of the whole
  panel. `--genes` takes a comma-separated list; `--gene-list` takes a file, either
  one gene per line or a table with a gene column and an optional panel/species
  column so a different subset applies to each organism. Names match
  case-insensitively; a name matching nothing is reported with the panel's full
  gene list rather than silently dropped.
- `--preset` sets how much divergence to tolerate when locating a gene:
  `asm5` / `asm10` (default) / `asm20`. `misfit --help` explains when to change it,
  with an example for each. See *Choosing a preset* below.
- `misfit-multi` runs a mixed-species cohort: it routes each assembly to the panel
  for its own organism and works through them in parallel, taking the same routing,
  preset and gene-subset options as `misfit`. It calls the package's own routing
  functions, so there is no second copy to drift.
- `python -m misfit.scripts.upstream` is an **experimental**, opt-in module for
  insertions upstream of a gene. It never runs as part of a normal MISFIT run.

### Added — four commands instead of one

Everything is now installed as a console script, so nothing has to be run out of a
checkout: `misfit` (call one panel), `misfit-multi` (a mixed-species cohort, in
parallel), `misfit-wide` (reshape results to one row per sample) and `misfit-db`
(curate the reference panels). `misfit-multi` was previously `scripts/run_misfit.py`
and `misfit-wide` was `scripts/misfit_wide.py`; both moved into the package, and the
`scripts/` directory is gone.

### Added — reshaping results

`misfit-wide` converts MISFIT's long form (one row per gene per assembly)
to wide form (one row per sample, one column per gene). Rows, columns and cell values
all default sensibly — `Assembly`, `Gene`, `Mutation_Description` — and each is an
option, so the same command produces a table of mutation types, HGVS nucleotide
descriptions or coverage figures.

Given the species map `misfit-multi` already takes, results are split per species,
either as one Excel workbook with a tab per species (the default) or one TSV per
species. This matters for mixed-organism cohorts because panels differ between
organisms, so a single combined table is ragged by construction; each split table
carries only its own species' genes, plus a `Species` column. Assemblies absent from
the map collect in an `unmapped` group rather than being dropped silently.

### Added — extending the panels

`misfit-db` is a second console script for curating the reference set:
`add-allele` puts another allele beside an existing gene, `add-gene` adds a gene to a
panel, and `add-panel` creates an organism panel seeded from an annotated assembly.
Each writes its `PROVENANCE.tsv` row as it goes, so the reference set cannot acquire an
unattributed sequence.

You supply a sequence; the tool resolves the
`gene|protein_accession|nuccore_accession|locus_tag` header from a genome annotation
(`--from-assembly`, `--from-nuccore`), from a protein accession you name
(`--protein-accession`), or — failing those — by BLAST, which only ever *reports*
candidates and exits without writing. Every resolved accession is verified ≥99%
identical to the submitted sequence before it is accepted, and candidates are screened
against the panel for ORF integrity, plausible length, duplication and paralogy.
`--unpublished` records an allele with no public accession.

### Choosing a preset

`asm10` remains the default and is right whenever the panel holds an allele close
to the isolate — the normal case now that panels carry per-species and
per-phylogroup alleles.

`asm20` is worth trying on any cohort containing organisms without a dedicated
allele. Testing on a real mixed-Enterobacterales collection showed it does more
than recover missing genes: where a gene is divergent enough that `asm10` aligns
only part of it, the partial alignment is reported as a **truncation**, which
looks like loss of function but is an alignment artifact. Aligning the full gene
resolves those into their real calls — and because the whole coding sequence
becomes visible, genuine frameshifts and premature stops that the partial
alignment had obscured are found as well. `Not found` and `Truncation` both fall
substantially; specific loss-of-function calls rise.

Nothing was lost in that comparison — no gene called under `asm10` became
`Not found` under `asm20`. If your cohort spans several species, run both and
compare, rather than assuming the default is conservative: on divergent input
`asm10` does not fail safe, it fails vague.

### Performance

minimap2 indexes its target, so calling it once per gene re-indexed the whole
assembly for every gene in the panel. Aligning the panel in a single pass returns
identical alignments and is roughly an order of magnitude faster per assembly; the
saving grows with the size of the panel. Reference files are also parsed once per
process rather than twice per assembly.

The change was validated by comparing per-gene and single-pass alignment across
every panel, and by confirming that whole-run output is unchanged field for field.

### Packaging

Declared the Biopython dependency (previously required but unlisted), removed a
package-data path pointing at a deleted directory, and made the version
single-sourced from `version.py`. `.gitignore` no longer swallows
`src/misfit/build/`. `src/misfit/bin/misfit` is gone, superseded by the `misfit`
console script.

### Known limitations

- Coding sequence only: a gene whose promoter carries an IS element is reported as
  an intact ORF — correctly, but incompletely.
- Divergence below ~90% nucleotide identity reads as `Not found`, which looks like
  absence.
- One locus per gene; a disrupted second copy is not reported.
- Partial hits still produce a call. Treat anything under ~90% `Nuc_COV` with care.

## v0.0.1

Initial release. Porin and cell-wall genes of *E. coli*.
