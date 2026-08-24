# MISFIT

**M**utation **I**dentification in **S**equences and **F**rameshift/**I**ndel **T**racking.

MISFIT answers one question: *how does this isolate's copy of gene X differ from a
reference copy of gene X?* It finds the gene in an assembly, compares it to a panel of
reference alleles for that organism, and reports the difference in HGVS notation.

It is aimed at the chromosomal genes behind beta-lactam and cefiderocol resistance in
Enterobacterales — PBPs, porins and their regulators, efflux systems, AmpC regulation
and iron-uptake receptors.

---

## What it reports

For every assembly x gene:

| column | meaning |
|---|---|
| `Gene_ID` | the reference allele the isolate matched |
| `Contig`, `Contig_Start`, `Contig_End`, `Strand` | where the gene sits in the assembly |
| `Mutation_Type` | WT, Missense, In-frame indel, Frameshift, Premature stop, Truncation, Not found |
| `Mutation_Description` | HGVS protein change, e.g. `p.Tyr334_Asn337dup; p.Ala413Val` |
| `HGVS_c` | the same change in nucleotide coordinates, e.g. `c.1000_1011dup` |
| `Notes` | canonical-allele numbering where it differs, plus cross-references to other tools |

Positions are **reference coordinates**: an insertion in the sample does not renumber
what follows it. An isolate carrying both the PBP3 333-loop duplication and the A413V
substitution reads `p.Tyr334_Asn337dup; p.Ala413Val`, never `p.Ala417Val`.

Descriptions follow HGVS — three-letter codes, `Ter` for termination, indels placed
3'-most, and an insertion that copies its 5' neighbour written as a duplication.

---

## Installation

> [!NOTE]
> Developed and tested on Linux (RHEL 7.9) and macOS. Other environments should work
> but are untested — please open an issue if you hit a platform problem.

```bash
git clone https://github.com/wshropshire/misfit
cd misfit
micromamba env create -f src/misfit/build/misfit_env.yaml   # or conda env create
micromamba activate misfit_env
pip install .
```

`minimap2` must be on `PATH` (the conda environment provides it).

---

## Usage

```bash
misfit assembly.fasta -o results.tsv --species "Klebsiella pneumoniae"
```

A directory holding several organisms is routed per assembly:

```bash
misfit assemblies/*.fasta -o results.tsv --species-map species_map.tsv
```

where `species_map.tsv` is a TSV or CSV with an assembly column
(`fasta_file` / `assembly` / `sample_id`) and a species column
(`species_wgs` / `species` / `organism`):

```
fasta_file                species_wgs
EC1092_contigs.fasta      Escherichia coli
KLPN2485_contigs.fasta    Klebsiella pneumoniae
ENTC1134_contigs.fasta    Enterobacter hormaechei
```

For a mixed-species cohort, `misfit-multi` routes each assembly to its own panel and
runs them in parallel:

```bash
misfit-multi assemblies/*.fasta \
    -o results.tsv --species-map species_map.tsv --jobs 8
```

```bash
misfit --list-organisms
```

**Arguments**

- `assemblies` — one or more assembly FASTA files
- `-o / --output` — output TSV
- `--species` — organism whose panel to use. Accepts a full species name
  (`Enterobacter hormaechei`) or a panel name (`k_pneumoniae`).
- `--species-map` — route a mixed-organism set per assembly; takes precedence over `--species`
- `--preset` — divergence tolerance when locating a gene: `asm5` / `asm10` (default) / `asm20`
- `--genes` — comma-separated genes of interest, instead of the whole panel
- `--gene-list` — a file of genes of interest, optionally per organism
- `--ref-dir` — use an arbitrary directory of reference FASTAs instead
- `--list-organisms` — list the available panels and exit

**Only the matched panel is read.** Genes from other organisms' panels are never used.

### Examining only some genes

By default every gene in the panel is examined. To narrow it:

```bash
misfit isolate.fasta -o out.tsv --species "E coli" --genes ftsI,ompC,ompF
```

Or from a file — one gene per line (`#` comments allowed):

```
ftsI
ompC
```

```bash
misfit isolate.fasta -o out.tsv --species "E coli" --gene-list genes.txt
```

Or a table, so a different subset applies to each organism. A gene with no
species applies to every panel:

```
gene     species
ftsI     Escherichia coli
ompC     Escherichia coli
ompK35   Klebsiella pneumoniae
ompK36   Klebsiella pneumoniae
acrB
```

```bash
misfit assemblies/*.fasta -o out.tsv --species-map map.tsv --gene-list genes.tsv
```

Gene names match case-insensitively. A gene named for another organism is simply
not selected for panels that lack it; a name that matches nothing anywhere is
reported with the list of genes the panel does hold.

### Alignment presets

A gene is located by aligning the panel's alleles to the assembly. `--preset` sets
how far the isolate's copy may have drifted from the nearest allele and still be
found. Set it too strict and the gene reads as `Not found`, which looks like
absence rather than divergence; too loose and alignment gets noisier at gene ends.

| preset | tolerates | when |
|---|---|---|
| `asm10` | ~10% divergence | **Default.** The panel holds an allele close to the isolate — the normal case, since panels carry per-species and per-phylogroup alleles. |
| `asm20` | ~20% | An organism with no dedicated allele in its panel: a congener routed to the genus panel, or a species listed in `reference/UNRESOLVED.tsv`. |
| `asm5` | ~5% | Confirming a near-exact match to a known allele. Too strict for general use. |

```bash
# default — an organism its panel represents
misfit isolate.fasta -o out.tsv --species "Klebsiella pneumoniae"

# a congener with no dedicated allele, routed to the genus panel
misfit isolate.fasta -o out.tsv --species "Citrobacter sedlakii" --preset asm20

# near-exact matching only
misfit isolate.fasta -o out.tsv --species "E coli" --preset asm5
```

On divergent input `asm20` recovers genes `asm10` misses outright, and it does not
make paralog mis-assignment more likely — where a diverged gene is lost, `asm10`
tends to return nothing rather than a safer answer. `asm5` loses genes at ordinary
within-species variation.

If unsure, run the default and re-run only the assemblies whose genes came back
`Not found` with `--preset asm20` to see whether they are recovered.

---

## Reference panels

Six panels, ~900 sequences, every one traceable to a RefSeq record:

| panel | genes | covers |
|---|---|---|
| `e_coli` | 56 | Escherichia, Shigella |
| `k_pneumoniae` | 54 | K. pneumoniae / quasipneumoniae / variicola |
| `k_oxytoca` | 51 | K. oxytoca complex |
| `k_aerogenes` | 54 | K. aerogenes |
| `enterobacter` | 43 | Enterobacter (any species) |
| `citrobacter` | 53 | Citrobacter (any species) |

Each gene file holds several alleles: one canonical allele from the panel's model
reference genome, plus representatives for the phylogroups (E. coli) or species within
the genus. MISFIT compares each isolate to whichever allele fits it best, so ordinary
lineage variation is not reported as mutation. `reference/PROVENANCE.tsv` records the
assembly, nucleotide, locus-tag and protein accession behind every sequence;
`reference/UNRESOLVED.tsv` lists genes deliberately left out because RefSeq does not
identify them.

Alleles are screened so a reference can never carry a resistance determinant: intact
ORF, expected length, no PBP3 333-loop insertion, no carbapenemase or ESBL in the source
genome, and modal across the screened genomes of that group.

---

## One row per sample — `misfit-wide`

MISFIT writes long form: one row per gene per assembly. For a summary table, a
supplementary file, or anything you intend to open in Excel, you usually want one row
per **sample** and one column per **gene**.

```bash
# whole cohort in one table
misfit-wide misfit_raw_calls.tsv -o wide.tsv

# mutation type rather than the HGVS description
misfit-wide misfit_raw_calls.tsv -o types.tsv --value Mutation_Type
```

Rows are keyed on `Assembly`, columns come from `Gene`, and cells hold
`Mutation_Description`. All three are options — `--index`, `--columns`, `--value` — so
`--value HGVS_c`, `--value Nuc_COV` or `--value Notes` give the same shape over a
different field.

### Splitting a mixed-organism run

Pass the same species map `misfit-multi` takes and the output is split per species,
which matters because panels differ between organisms: `ampC` is in some panels and not
others, so a single combined table is unavoidably ragged. Splitting gives each species
only its own genes.

```bash
# one Excel workbook, one tab per species (the default when a map is given)
misfit-wide misfit_raw_calls.tsv -o wide.xlsx \
    --species-map assembly_species_map.tsv

# the same split as one TSV per species in a directory
misfit-wide misfit_raw_calls.tsv -o wide_by_species/ \
    --species-map assembly_species_map.tsv --format tsv
```

Each split table gains a `Species` column immediately after the sample column. Genes
absent from a species' panel are dropped from that tab rather than filling it with
placeholder cells. Assemblies missing from the map are **not** discarded — they collect
in an `unmapped` tab or file, with a warning naming them.

**Other options**

- `--fill` — cell value where a sample has no row for a gene (default `NA`; distinct
  from `--`, which MISFIT itself writes for a gene it looked for and did not find)
- `--join` — separator when one sample/gene pair has several differing rows
  (default `" | "`; deliberately not `"; "`, which MISFIT uses *within* one call to
  separate variants)
- `--species-column` — rename the added column
- `--format` — `excel` or `tsv`; defaults to `excel` with a species map, `tsv` without

Writing `.xlsx` needs `openpyxl`; TSV output needs nothing beyond the standard library.

---

## Extending the panels — `misfit-db`

`misfit-db` adds an allele, a gene or a whole organism panel, and writes the
`PROVENANCE.tsv` row at the same time so nothing enters the reference set unattributed.

```bash
# an allele of a gene already in a panel
misfit-db add-allele --panel "E coli" --gene ompC --sequence my_ompC.fasta \
    --protein-accession WP_000865587.1

# a gene the panel does not yet carry
misfit-db add-gene --panel enterobacter --gene ompX --sequence ompX.fasta \
    --from-nuccore CP049085.2

# a new organism, seeded from an annotated assembly
misfit-db add-panel --panel serratia --aliases '^serratia' \
    --from-assembly GCF_000513215.1 --genes ompC,acrR,ampD
```

`--sequence` takes a FASTA file or a bare nucleotide string.

### Resolving the header

Reference headers are `gene|protein_accession|nuccore_accession|locus_tag`. You supply
the sequence; `misfit-db` works out the rest by whichever route you give it, in
decreasing order of authority:

| route | what it does |
|---|---|
| `--from-assembly GCF_…` | reads the assembly's annotation, finds that gene, compares |
| `--from-nuccore CP…` / `NC_…` | same, from one nucleotide record |
| `--protein-accession WP_…` | fetches the protein and verifies it against your sequence |
| *(nothing)* | BLASTs the translation and **reports candidates without writing** |
| `--unpublished` | records an allele that has no public accession |

Every route that resolves an accession also checks it is ≥99% identical to your
sequence and refuses otherwise, so a typo cannot attach the wrong accession to a
sequence. `--force` overrides; `--dry-run` shows the header and writes nothing.

For a protein accession, a `WP_` multispecies record is the right thing to cite for an
allele — it identifies the protein independently of any one genome. Give `--nuccore`
and `--locus-tag` as well if you want the row to point at a specific genome too.

### If BLAST is slow

Identification with no accession goes through NCBI's public BLAST queue, which is
frequently far slower than the time NCBI itself estimates — a single 373 aa OmpC query
took **14.4 minutes** against an NCBI estimate of 45 seconds. `misfit-db` therefore waits
30 minutes by default (`--blast-timeout`). If it does give up, the job is not lost: it
prints the RID, and `--blast-rid <RID>` collects the finished result without submitting
again.

If you already know what your sequence is, `--protein-accession` skips the queue
entirely and is the faster route by a wide margin.

BLAST identification is always advisory: it prints the candidates and exits without
writing, so an accession only ever enters the panel because you chose it.

### Screening

Before writing, a candidate is checked against the panel: an intact ORF, a length
consistent with the gene's other alleles, not already present, and not a paralog. Real
problems are reported and the write refused. Adding *E. coli* UTI89 `ompC` to the
`e_coli` panel is rejected at 77.3% identity as a paralog rather than an allele.

---

## Experimental

`python -m misfit.scripts.upstream` looks for insertions upstream of a gene, using the
cohort's modal architecture as the wild type. **Not part of a normal run.** It works
where the upstream region is conserved and cannot resolve loci where it is not.

---

## Limitations

- **Coding sequence only.** A gene whose promoter carries an IS element is reported as
  an intact ORF — correctly, but incompletely.
- **Divergence shows up as missing data.** minimap2 needs ~90% nucleotide identity;
  below that a gene reads `Not found`, which looks like absence.
- **One locus per gene.** If a genome carries a disrupted second copy, it is not reported.

`ARCHITECTURE.md` describes how the code works and where the rough edges are;
`CHANGELOG.md` covers what changed in this release.

---

## License

[MIT](LICENSE.txt)

---

## Contributing

Feel free to contribute to the project by submitting issues or pull requests.

---

## Version

misfit v1.0.0

## TODO

- [x] Make reference database generation more flexible and customizable — `misfit-db`
      adds alleles, genes and whole panels, resolving accessions and writing provenance
- [x] Simulate data to validate — `tests/test_variants.py` builds queries by editing a
      real CDS at known codons, so the expected call is known exactly rather than eyeballed
- [ ] Test edge cases where identity/coverage fail to identify best hit. Partly done:
      selection now scores absolute matching bases, which fixed a truncated allele
      winning by being short. Partial hits under ~90% `Nuc_COV` still produce confident calls.
- [ ] Report a gene split by a large insertion as a disruption rather than a truncation
      of whichever fragment aligns better, and anchor the coordinates to the 5' fragment
- [ ] Measure truncation against the matched allele's length, not the longest allele in the file
- [ ] Detect start-loss (`p.Met1?`) as its own category
- [ ] Flag partial hits as unassessable instead of calling them

