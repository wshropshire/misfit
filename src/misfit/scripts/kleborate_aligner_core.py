"""Adapted from Kleborate script (https://github.com/klebgenomics/Kleborate/blob/main/kleborate/shared/alignment.py) 
Uses Minimap2 aligner + culling best hit match to avoid redundancy
"""

import os
import re
import subprocess
import shutil
import sys
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
from .misc import load_fasta, reverse_complement

class Alignment:
    def __init__(self, paf_line, query_seqs=None, ref_seqs=None):
        self.query_name, self.query_length = None, None
        self.query_start, self.query_end = None, None
        self.strand = None
        self.ref_name, self.ref_length = None, None
        self.ref_start, self.ref_end = None, None
        self.matching_bases, self.num_bases = None, None
        self.percent_identity = None
        self.query_cov, self.ref_cov = None, None
        self.cigar, self.alignment_score = None, None
        self.query_seq, self.ref_seq = None, None
        self.ref_header = None
        self.parse_paf_line(paf_line)
        self.set_identity_and_coverages()
        self.set_sequences(query_seqs, ref_seqs)

    def parse_paf_line(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        if len(line_parts) < 11:
            sys.exit('Error: alignment file does not seem to be in PAF format')

        self.query_name = line_parts[0]
        self.query_length = int(line_parts[1])
        self.query_start = int(line_parts[2])
        self.query_end = int(line_parts[3])
        self.strand = line_parts[4]

        self.ref_name = line_parts[5]
        self.ref_length = int(line_parts[6])
        self.ref_start = int(line_parts[7])
        self.ref_end = int(line_parts[8])

        self.matching_bases = int(line_parts[9])
        self.num_bases = int(line_parts[10])

        for part in line_parts:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]
            if part.startswith('AS:i:'):
                self.alignment_score = int(part[5:])

    def set_identity_and_coverages(self):
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases
        self.query_cov = 100.0 * (self.query_end - self.query_start) / self.query_length
        self.ref_cov = 100.0 * (self.ref_end - self.ref_start) / self.ref_length

    def set_sequences(self, query_seqs, ref_seqs):
        if query_seqs is not None:
            self.query_seq = query_seqs[self.query_name][self.query_start:self.query_end]
        if ref_seqs is not None:
            self.ref_seq = ref_seqs[self.ref_name][self.ref_start:self.ref_end]
            self.ref_header = self.ref_name  # default fallback
            for name, seq in ref_seqs.items():
                if self.ref_name == name:
                    self.ref_header = name  # replace with original header key
                    break
            if self.strand == '-':
                self.ref_seq = reverse_complement(self.ref_seq)

    def __repr__(self):
        return f"{self.query_name}:{self.query_start}-{self.query_end}({self.strand}), " \
               f"{self.ref_name}:{self.ref_start}-{self.ref_end} ({self.percent_identity:.3f}%)"

def get_minimap2_path():
    """
    Detects minimap2 binary in local misfit/binaries folder or system PATH.
    """
    # Check local binaries folder (scripts/../binaries)
    local_bin = os.path.join(os.path.dirname(__file__), "..", "binaries", "minimap2")
    local_bin = os.path.abspath(local_bin)
    if os.path.exists(local_bin) and os.access(local_bin, os.X_OK):
        return local_bin

    # Fallback to system PATH
    minimap2_path = shutil.which("minimap2")
    if minimap2_path is None:
        raise RuntimeError("ERROR: minimap2 not found in misfit/binaries or in system PATH!")
    return minimap2_path

def align_query_to_ref(query_filename, ref_filename, ref_index=None, preset='map-ont',
                        min_identity=None, min_query_coverage=None):
    query_seqs = dict(load_fasta(query_filename))
    ref_seqs = dict(load_fasta(ref_filename))
    ref = ref_filename if ref_index is None else ref_index
    with open(os.devnull, 'w') as dev_null:
        minimap2_exec = get_minimap2_path()
        out = subprocess.check_output([
            minimap2_exec, '--end-bonus=10', '--eqx', '-c', '-N', '10', '--secondary=yes', '-x', preset,
            str(ref), str(query_filename)
        ], stderr=dev_null)
    alignments = [Alignment(x, query_seqs=query_seqs, ref_seqs=ref_seqs)
                  for x in out.decode().splitlines()]
    print(f"→ Loaded {len(ref_seqs)} reference sequences: {[k for k in ref_seqs.keys()]}")
    print(f"→ Minimap2 returned {len(alignments)} alignments:")
    for aln in alignments:
        print(f"  Ref: {aln.ref_name}, Identity: {aln.percent_identity:.1f}, Coverage: {aln.query_cov:.1f}")

    if min_identity is not None:
        alignments = [a for a in alignments if a.percent_identity >= min_identity]
    if min_query_coverage is not None:
        alignments = [a for a in alignments if a.query_cov >= min_query_coverage]
    return alignments

def cull_redundant_hits(minimap_hits, w_id=0.6, w_cov=0.4):
    # Score = weighted identity + coverage; default: 60% ID + 40% coverage
    sorted_hits = sorted(
        minimap_hits,
        key=lambda x: (w_id * x.percent_identity + w_cov * x.query_cov),
        reverse=True
    )

    filtered_hits = []
    for h in sorted_hits:
        if not overlapping(h, filtered_hits):
            filtered_hits.append(h)

    if filtered_hits:
        best = filtered_hits[0]
    else:
        best = sorted_hits[0]

    print("\n→ Selected best hit:")
    print(f"  Ref: {best.ref_name}")
    print(f"  Identity: {best.percent_identity:.1f}%")
    print(f"  Coverage: {best.query_cov:.1f}%")
    return best

def overlapping(hit, existing_hits):
    for existing_hit in existing_hits:
        if hit.strand == existing_hit.strand and hit.ref_name == existing_hit.ref_name:
            if hits_overlap(hit, existing_hit):
                return True
    return False

def hits_overlap(a, b):
    return a.ref_start <= b.ref_end and b.ref_start <= a.ref_end and \
           len(range(max(a.ref_start, b.ref_start), min(a.ref_end, b.ref_end))) > 50
