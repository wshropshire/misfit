"""Adapted from Kleborate script (https://github.com/klebgenomics/Kleborate/blob/main/kleborate/shared/alignment.py) 
Uses Minimap2 aligner + culling best hit match to avoid redundancy
"""

import os
import re
import subprocess
import shutil
import sys
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

    # MISFIT runs minimap2 with the gene panel as the QUERY and the assembly as
    # the TARGET, so the PAF-correct field names read backwards here. These
    # aliases name the roles as MISFIT actually uses them.
    @property
    def allele_name(self):
        """The reference allele this alignment used."""
        return self.query_name

    @property
    def allele_cov(self):
        """Percent of that reference allele covered by the alignment."""
        return self.query_cov

    @property
    def sample_contig(self):
        """The assembly contig the gene was found on."""
        return self.ref_name

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
    print(f"→ Query: {len(query_seqs)} reference allele(s) from {os.path.basename(str(query_filename))}")
    print(f"→ Subject: {len(ref_seqs)} contig(s) from {os.path.basename(str(ref_filename))}")
    print(f"→ Minimap2 returned {len(alignments)} alignments:")
    for aln in alignments:
        print(f"  allele={aln.allele_name}  contig={aln.sample_contig}  "
              f"identity={aln.percent_identity:.1f}%  allele_cov={aln.allele_cov:.1f}%  "
              f"matching_bases={aln.matching_bases}")

    if min_identity is not None:
        alignments = [a for a in alignments if a.percent_identity >= min_identity]
    if min_query_coverage is not None:
        alignments = [a for a in alignments if a.query_cov >= min_query_coverage]
    return alignments

def cull_redundant_hits(minimap_hits):
    """Pick the best alignment, scoring by how much sequence actually matched.

    Coverage expressed as a fraction of each reference allele's own length
    lets a short allele win by being short: a 435 bp truncated ompK35 scores
    100% coverage against an isolate whose gene is full length, beating the
    real 1080 bp allele and turning an intact porin into a frameshift call.
    Absolute matching bases has no such bias.

    A consequence worth knowing: identity alone does not decide the winner.
    Identity is matching_bases / alignment_length, so a longer allele matching
    more bases beats a shorter one at marginally higher identity -- 94.2% over
    1131 bases beats 94.4% over 1104. That is deliberate; more of the gene was
    actually matched. The diagnostic prints matching_bases so the choice can be
    checked rather than inferred.
    """
    sorted_hits = sorted(
        minimap_hits,
        key=lambda x: (x.matching_bases, x.percent_identity),
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

    print("\n→ Selected best hit (most matching bases, then identity):")
    print(f"  allele        : {best.allele_name}")
    print(f"  contig        : {best.sample_contig}:{best.ref_start + 1}-{best.ref_end} "
          f"({best.strand})")
    print(f"  matching_bases: {best.matching_bases}   identity: {best.percent_identity:.1f}%   "
          f"allele_cov: {best.allele_cov:.1f}%")
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
