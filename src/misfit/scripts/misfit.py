# Mutation Identification in Sequences and Frameshift/ Indel Tracking (MISFIT)
# misfit.py <assemblies> -o <tsv output>

import argparse
import os
import warnings
import logging
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.Data.CodonTable import TranslationError
from difflib import SequenceMatcher
from pathlib import Path
from .kleborate_aligner_core import align_query_to_ref, cull_redundant_hits
from .version import __version__

# Suppress deprecation warning for pairwise2
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=DeprecationWarning)

REF_GENES = {
    "ompC": "ref_ompC.fasta",
    "ompF": "ref_ompF.fasta",
    "phoE": "ref_phoE.fasta",
    "envZ": "ref_envZ.fasta",
    "ompR": "ref_ompR.fasta",
    "ftsI": "ref_ftsI.fasta",
    "mrdA": "ref_mrdA.fasta"
}

def clean_seq(seq):
    allowed = set("ACGT")
    for b in set(seq.upper()) - allowed:
        seq = seq.split(b)[0]
    return seq[:len(seq) // 3 * 3]

def get_expected_ref_aa_length(ref_path: str) -> int:
    """Return maximum AA length from reference translations, excluding trailing stop codons."""
    max_len = 0
    for record in SeqIO.parse(ref_path, "fasta"):
        try:
            aa_seq = str(record.seq.translate(table="Bacterial", to_stop=False)).rstrip("*")
            max_len = max(max_len, len(aa_seq))
        except Exception:
            continue
    return max_len

def detect_mutation_type(query_nucl: str, ref_nucl: str, gene_name: str = "", ref_aa_len: int = None):
    query_nucl = clean_seq(query_nucl)
    ref_nucl = clean_seq(ref_nucl)

    try:
        query_aa = str(Seq(query_nucl).translate(table="Bacterial", to_stop=False))
        ref_aa = str(Seq(ref_nucl).translate(table="Bacterial", to_stop=False))
    except TranslationError as e:
        return "translation_error", None, str(e)

    print(f"\n→ Checking gene: {gene_name}")
    print(f"→ Query nucleotide length: {len(query_nucl)}")
    print(f"→ Reference nucleotide length: {len(ref_nucl)}")
    print(f"→ Translated query AA length: {len(query_aa)}")
    print(f"→ Translated ref AA length: {len(ref_aa)}")

    if len(query_nucl) % 3 != 0:
        mismatch_pos = next((i for i, (a, b) in enumerate(zip(query_aa, ref_aa)) if a != b), None)
        context = query_aa[max(mismatch_pos - 3, 0): mismatch_pos + 3] if mismatch_pos is not None else "NA"
        return "frameshift", mismatch_pos + 1 if mismatch_pos is not None else None, f"FS{mismatch_pos + 1}: {context}"

    internal_stops = query_aa.count("*") - (1 if query_aa.endswith("*") else 0)
    if internal_stops > 0:
        stop_pos = query_aa.find("*")
        tail_query = query_aa[stop_pos + 1:]
        tail_ref = ref_aa[stop_pos + 1:]
        similarity = SequenceMatcher(None, tail_query, tail_ref).ratio()
        if internal_stops == 1 and similarity >= 0.6:
            aa_ref = ref_aa[stop_pos] if stop_pos < len(ref_aa) else "?"
            return "premature_stop", stop_pos + 1, f"p.{aa_ref}{stop_pos + 1}X"
        else:
            fs_start = next((i for i, (a, b) in enumerate(zip(query_aa, ref_aa)) if a != b), stop_pos)
            aa_ref = ref_aa[fs_start] if fs_start < len(ref_aa) else "?"
            aa_query = query_aa[fs_start] if fs_start < len(query_aa) else "?"
            fs_len = max(1, stop_pos - fs_start)
            return "frameshift", fs_start + 1, f"p.{aa_ref}{fs_start + 1}{aa_query}fsX{fs_len}"

    if ref_aa_len is not None:
        length_ratio = len(query_aa) / ref_aa_len if ref_aa_len > 0 else 0
        if length_ratio < 0.9:
            return "truncation", len(query_aa), f"Truncated protein: {len(query_aa)}/{ref_aa_len} AA"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        alignments = pairwise2.align.globalms(ref_aa, query_aa, 2, -2, -3, -0.5)
    ref_aln, query_aln, score, _, _ = alignments[0]

    print("\n=== Protein Alignment ===")
    print("REF   :", ref_aln)
    print("QUERY :", query_aln)
    print("SCORE :", score)
    print("=========================")

    indels = []
    missense_list = []
    i = 0
    while i < len(ref_aln):
        if query_aln[i] == '-' and ref_aln[i] != '-':
            del_start = i
            deleted = ''
            while i < len(ref_aln) and query_aln[i] == '-':
                deleted += ref_aln[i]
                i += 1
            prev_aa = ref_aln[del_start - 1] if del_start > 0 else "?"
            indels.append(f"p.{prev_aa}{del_start}del{deleted}")
        elif ref_aln[i] == '-' and query_aln[i] != '-':
            ins_start = i
            inserted = ''
            while i < len(query_aln) and ref_aln[i] == '-':
                inserted += query_aln[i]
                i += 1
            prev_aa = query_aln[ins_start - 1] if ins_start > 0 else "?"
            indels.append(f"p.{prev_aa}{ins_start}ins{inserted}")
        elif query_aln[i] != ref_aln[i] and query_aln[i] != '-' and ref_aln[i] != '-':
            missense_list.append(f"p.{ref_aln[i]}{i+1}{query_aln[i]}")
            i += 1
        else:
            i += 1

    detect_mutation_type.ref_aln = ref_aln
    detect_mutation_type.query_aln = query_aln

    if indels and missense_list:
        return "inframe_indel_missense", None, "; ".join(indels + missense_list)
    elif indels:
        return "inframe_indel", None, "; ".join(indels)
    elif missense_list:
        return "missense", None, "; ".join(missense_list)

    return "intact", None, "Recovered intact: no AA diff"

def match_best_hit_to_reference(best_hit_ref_seq, ref_path, debug=False):
    ref_records = list(SeqIO.parse(ref_path, "fasta"))
    best_seq_str = str(best_hit_ref_seq).upper()
    best_match, best_similarity = "NA", 0.0
    for rec in ref_records:
        ref_seq_str = str(rec.seq).upper()
        similarity = SequenceMatcher(None, ref_seq_str, best_seq_str).ratio()
        if similarity > best_similarity:
            best_similarity = similarity
            best_match = rec.id
    if debug:
        print(f"→ Best match: {best_match} (similarity: {best_similarity:.3f})")
    return best_match

def detect_mutations(assembly_file, output_file, ref_dir):
    output_lines = []
    asm_id = os.path.basename(assembly_file).split(".")[0]
    
    for gene, ref_file in REF_GENES.items():
        ref_path = os.path.join(ref_dir, ref_file)
        hits = align_query_to_ref(ref_path, assembly_file, preset='asm10')

        if not hits:
            output_lines.append(f"{asm_id}\t{gene}\t--\t--\t--\t--\tNot found\t--\t--\t--\t--")
            continue

        best_hit = cull_redundant_hits(hits)
        query_nucl = clean_seq(best_hit.query_seq)
        ref_nucl = clean_seq(best_hit.ref_seq)

        try:
            query_aa = str(Seq(query_nucl).translate(table="Bacterial", to_stop=False)).rstrip("*")
            ref_aa = str(Seq(ref_nucl).translate(table="Bacterial", to_stop=False)).rstrip("*")
        except TranslationError:
            query_aa, ref_aa = "", ""

        aa_len = len(query_aa)
        ref_aa_len = len(ref_aa)

        # Estimate max expected length from reference pool
        expected_ref_aa_len = get_expected_ref_aa_length(ref_path)

        # Mutation detection
        mutation_type, position, desc = detect_mutation_type(
            ref_nucl, query_nucl, gene_name=gene, ref_aa_len=expected_ref_aa_len
        )

        # Adjust AA coverage for truncation, premature stop, etc.
        if mutation_type in ["frameshift", "premature_stop", "truncation"] and position is not None:
            aa_cov = (position / expected_ref_aa_len * 100) if expected_ref_aa_len else 0
        else:
            aa_cov = (aa_len / expected_ref_aa_len * 100) if expected_ref_aa_len else 0

        # Use alignment to calculate AA identity
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            alignments = pairwise2.align.globalms(ref_aa, query_aa, 2, -2, -3, -0.5)

        ref_aln, query_aln, score, _, _ = alignments[0]
        aa_matches = sum(1 for r, q in zip(ref_aln, query_aln) if r == q and r != '-' and q != '-')
        aligned_len = sum(1 for r, q in zip(ref_aln, query_aln) if r != '-' and q != '-')
        aa_id = (aa_matches / aligned_len * 100) if aligned_len else 0

        gene_id = match_best_hit_to_reference(best_hit.ref_seq, ref_path, debug=False)
        print(f"→ Final best hit for {gene}: {gene_id}")

        msg = {
            "frameshift": "Frameshift",
            "premature_stop": "Premature stop",
            "translation_error": "Translation error",
            "inframe_indel": "In-frame indel",
            "missense": "Missense",
            "inframe_indel_missense": "In-frame indel and Missense",
            "truncation": "Truncation",
            "intact": "WT"
        }.get(mutation_type, "WT")

        nuc_len = len(best_hit.query_seq)

        line = (
            f"{asm_id}\t{gene}\t{gene_id}\t{nuc_len}\t"
            f"{best_hit.percent_identity:.1f}%\t{best_hit.query_cov:.1f}%\t"
            f"{aa_len}\t{aa_id:.1f}%\t{aa_cov:.1f}%\t{msg}\t{desc}"
        )
        output_lines.append(line)

    write_mode = 'a' if os.path.exists(output_file) else 'w'
    with open(output_file, write_mode) as f:
        if write_mode == 'w':
            f.write("Assembly\tGene\tGene_ID\tNuc_Len\tNuc_ID\tNuc_COV\tAA_Len\tAA_ID\tAA_Cov\tMutation_Type\tMutation_Description\n")
        f.write("\n".join(output_lines) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Detect mutations in porin-related genes in E. coli genomes")
    parser.add_argument("assemblies", nargs='+', help="Input assembly FASTA file(s)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV summary file")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}", help="Show program's version number and exit")
    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        level=logging.INFO
    )

    script_dir = Path(__file__).resolve().parent
    ref_dir = script_dir.parent / "reference"

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
        detect_mutations(fasta, args.output, str(ref_dir))

    logging.info("Mutation detection completed successfully.")

if __name__ == "__main__":
    main()
