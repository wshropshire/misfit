# Mutation Identification in Sequences and Frameshift/ Indel Tracking (MISFIT)

## Overview

The Mutation Identification in Sequences and Frameshift/ Indel Tracking (MISFIT) script detects mutations in reference genes, currently focusing on porin-related and cell wall–associated genes of *E. coli* from whole-genome assemblies. It reports mutation types (e.g., frameshift, in-frame indel, premature stop, missense) with nucleotide and amino acid statistics in a TSV summary file.

---

## **Features**
- Detects mutations in the following genes (currently):
  - `ompC`, `ompF`, `phoE`, `envZ`, `ompR`, `ftsI`, `mrdA`
- Identifies:
  - Frameshifts
  - Premature stop codons
  - In-frame insertions/deletions
  - Missense mutations
  - Truncations or intact genes
- Calculates:
  - Nucleotide and amino acid coverage
  - Amino acid identity
  - Mutation descriptions in HGVS-like notation

---

## Installation

**NOTE** I've only tested this with a Linux, RHEL 7.9 operating system. I am not sure how this will function in other OS environments. The easiest way to install is:

(1) Clone the GitHub repository and then (2) create a conda environment. 

```bash
git clone https://github.com/wshropshire/misfit
cd misfit
# Create a conda environment
conda create -n misfit_env python=3.9
conda activate misfit_env
pip3 install build
# Build a sparse misfit package
python -m build
pip3 install ./dist/misfit-0.0.1.tar.gz
# Use yaml file to build environment with all dependencies - issues with pip/conda install 'path collisions' due to pre-installed python. All dependencies are properly downloaded and identified.
conda env update --name misfit_env --file ./src/misfit/build/misfti_env.yaml
```

---

## **Input Files**
1. **Genome assemblies** in FASTA format (`.fasta` or `.fa`)
2. **Reference sequences** located in a `reference/` directory alongside the script:
   ```
   ref_ompC.fasta
   ref_ompF.fasta
   ref_phoE.fasta
   ref_envZ.fasta
   ref_ompR.fasta
   ref_ftsI.fasta
   ref_mrdA.fasta
   ```

---

## **Usage**
```bash
python detect_mutations.py assemblies/*.fasta -o mutations_summary.tsv
```

**Arguments:**
- `assemblies`: One or more genome assembly FASTA files
- `-o / --output`: Path to output TSV summary file

---

## **Output**
A tab-delimited summary file with the following columns:

| Assembly | Gene | Gene_ID | Nuc_Len | Nuc_ID | Nuc_COV | AA_Len | AA_ID | AA_Cov | Mutation_Type | Mutation_Description |
|---------|------|--------|--------|-------|--------|-------|------|-------|----------------|---------------------|

Example:

```
Assembly   Gene   Gene_ID    Nuc_Len Nuc_ID Nuc_COV AA_Len AA_ID AA_Cov Mutation_Type      Mutation_Description
EC01       ompC   OmpC_Ref1  912     99.5%  100.0%  304    99.0% 100.0% WT                 Recovered intact: no AA diff
EC01       ompF   OmpF_Ref1  1020    98.3%  98.0%   340    95.0%  96.5%  Frameshift         p.Y254EfsX3
```

---

## **Notes**
- For frameshift or truncation mutations, amino acid coverage reflects the aligned portion before disruption.
- Internal stop codons are annotated as `premature_stop` if the remainder of the sequence remains highly similar to the reference.
- The script automatically appends results to an existing output file if rerun with the same `-o` argument.

## License

[MIT License](LICENSE.txt)

---

## Contributing

Feel free to contribute to the project by submitting issues or pull requests.

---

## Version

misfit v0.0.1
