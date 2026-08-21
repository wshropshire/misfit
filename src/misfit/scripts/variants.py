"""Reference-anchored protein variant calling.

Every coordinate produced here is a position in the *reference* protein.
An insertion carried by the sample therefore does not renumber the
substitutions that follow it: an isolate with both the PBP3 333-loop
duplication and the A413V substitution is reported as
``p.P333insYRIN; p.A413V``, never as ``p.A417V``.

Descriptions follow HGVS: three-letter amino-acid codes, ``Ter`` for
termination, indels placed 3'-most, and an insertion that copies the sequence
immediately 5' of it written as a duplication. Positions follow the reference
allele the isolate actually matched; the canonical-allele equivalent is carried
separately in ``VariantResult.canonical``. So the PBP3 333-loop event is
``p.Tyr334_Asn337dup`` and the OmpK36 loop-3 event is ``p.Gly134_Asp135dup``.

Frame analysis still uses 5'-most placement internally, because there the
earliest possible break is what determines the reading frame.
"""

import re
from dataclasses import dataclass, field

from Bio import Align
from Bio.Seq import Seq

CODON = 3
_ZERO = re.match(r"(0)", "0")

AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
       "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
       "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
       "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa"}


def aa3(seq):
    """One-letter to HGVS three-letter."""
    return "".join(AA3.get(c, "Xaa") for c in seq)

# Affine-gap protein scoring; mismatches cost less than opening a gap, so a
# substitution is preferred over an indel unless the evidence is clear.
_AA_ALIGNER = Align.PairwiseAligner(
    mode="global", match_score=2, mismatch_score=-2,
    open_gap_score=-3, extend_gap_score=-0.5,
)
_NT_ALIGNER = Align.PairwiseAligner(
    mode="global", match_score=2, mismatch_score=-3,
    open_gap_score=-6, extend_gap_score=-1,
)


@dataclass
class VariantResult:
    kind: str
    ref_pos: int = None
    description: str = ""
    variants: list = field(default_factory=list)
    query_aa_len: int = 0
    ref_aa_len: int = 0
    aa_identity: float = 0.0
    hgvs_c: str = "c.="
    canonical: list = field(default_factory=list)
    notes: list = field(default_factory=list)


def trim_to_codon(seq):
    return seq[: len(seq) // CODON * CODON]


def translate_cds(nt):
    """Translate to protein. Ambiguous codons become 'X' rather than truncating."""
    return str(Seq(trim_to_codon(nt)).translate(table="Bacterial"))


def _aligned_pair(aligner, ref, qry):
    aln = aligner.align(ref, qry)[0]
    return aln[0], aln[1]


def _left_normalize(ref_aln, qry_aln):
    """Shift every gap run as far 5' as the flanking sequence allows.

    A gap run may move one column left whenever the residue leaving the run
    on the right equals the residue entering it on the left, which leaves
    both ungapped sequences untouched.
    """
    ref_aln, qry_aln = list(ref_aln), list(qry_aln)
    for gapped, other in ((ref_aln, qry_aln), (qry_aln, ref_aln)):
        i = 0
        while i < len(gapped):
            if gapped[i] != "-":
                i += 1
                continue
            start = i
            while i < len(gapped) and gapped[i] == "-":
                i += 1
            end = i
            # Never shift a run across a gap in the opposite sequence.
            while (start > 0 and gapped[start - 1] != "-"
                   and other[start - 1] != "-" and other[end - 1] != "-"
                   and other[start - 1] == other[end - 1]):
                gapped[end - 1] = gapped[start - 1]
                gapped[start - 1] = "-"
                start -= 1
                end -= 1
    return "".join(ref_aln), "".join(qry_aln)


def _right_normalize(ref_aln, qry_aln):
    """Shift every gap run as far 3' as the sequence allows.

    HGVS requires the 3'-most representation of an indel. A run may move one
    column right whenever the residue entering on the left equals the one
    leaving on the right, which leaves both ungapped sequences unchanged.
    """
    ref_aln, qry_aln = list(ref_aln), list(qry_aln)
    for gapped, other in ((ref_aln, qry_aln), (qry_aln, ref_aln)):
        i = len(gapped) - 1
        while i >= 0:
            if gapped[i] != "-":
                i -= 1
                continue
            end = i + 1
            while i >= 0 and gapped[i] == "-":
                i -= 1
            start = i + 1
            while (end < len(gapped) and gapped[end] != "-"
                   and other[start] != "-" and other[end] != "-"
                   and other[start] == other[end]):
                gapped[start] = gapped[end]
                gapped[end] = "-"
                start += 1
                end += 1
    return "".join(ref_aln), "".join(qry_aln)


def position_map(from_aa, to_aa):
    """Map 1-based positions in one allele onto another.

    Multi-allele panels compare each isolate to whichever allele fits best, but
    those alleles differ in length, so a position means different things in
    different ones. Everything is therefore renumbered onto the panel's
    canonical allele before it is reported. Positions with no counterpart --
    residues inside an allele-specific insertion -- map to None.
    """
    if from_aa == to_aa:
        return None                      # identity, nothing to do
    ref_aln, qry_aln = _aligned_pair(_AA_ALIGNER, to_aa, from_aa)
    out = {}
    to_pos = from_pos = 0
    for a, b in zip(ref_aln, qry_aln):
        if a != "-":
            to_pos += 1
        if b != "-":
            from_pos += 1
            out[from_pos] = to_pos if a != "-" else None
    return out


def _ins_or_dup(ref_aa, anchor, inserted, N=None, pos_offset=0):
    """HGVS: an insertion that copies the sequence just 5' of it is a duplication."""
    N = N or (lambda p: p)
    n = len(inserted)
    if anchor >= n and ref_aa[anchor - n:anchor] == inserted:
        first, last = anchor - n + 1, anchor
        if n == 1:
            return f"p.{aa3(ref_aa[first-1])}{N(first + pos_offset)}dup"
        return (f"p.{aa3(ref_aa[first-1])}{N(first + pos_offset)}_"
                f"{aa3(ref_aa[last-1])}{N(last + pos_offset)}dup")
    a, b = anchor, anchor + 1
    return (f"p.{aa3(ref_aa[a-1])}{N(a + pos_offset)}_"
            f"{aa3(ref_aa[b-1]) if b <= len(ref_aa) else 'Ter'}{N(b + pos_offset)}"
            f"ins{aa3(inserted)}")


def call_protein_variants(ref_aa, qry_aa, pos_offset=0, renumber=None):
    """Walk a protein alignment and emit variants in reference coordinates.

    ``pos_offset`` is added to every position, for the case where the
    reference CDS was only partly covered by the alignment. ``renumber`` maps
    matched-allele positions onto the panel's canonical allele; positions with
    no canonical counterpart are reported unmapped and flagged by the caller.
    """
    def C(p):
        """Same position expressed against the panel's canonical allele."""
        if renumber is None:
            return p
        return renumber.get(p - pos_offset, p) or p
    ref_aln, qry_aln = _right_normalize(*_aligned_pair(_AA_ALIGNER, ref_aa, qry_aa))

    substitutions, indels, stops = [], [], []
    ref_pos = pos_offset          # reference residues consumed so far
    qry_seen = 0
    i = 0
    while i < len(ref_aln):
        r, q = ref_aln[i], qry_aln[i]

        if r != "-" and q == "-":
            first = ref_pos + 1
            deleted = ""
            while i < len(ref_aln) and qry_aln[i] == "-" and ref_aln[i] != "-":
                deleted += ref_aln[i]
                ref_pos += 1
                i += 1
            if len(deleted) == 1:
                indels.append((f"p.{aa3(deleted)}{first}del",
                               f"p.{aa3(deleted)}{C(first)}del"))
            else:
                indels.append((f"p.{aa3(deleted[0])}{first}_{aa3(deleted[-1])}{ref_pos}del",
                               f"p.{aa3(deleted[0])}{C(first)}_{aa3(deleted[-1])}{C(ref_pos)}del"))

        elif r == "-" and q != "-":
            inserted = ""
            while i < len(ref_aln) and ref_aln[i] == "-":
                inserted += qry_aln[i]
                qry_seen += 1
                i += 1
            indels.append((_ins_or_dup(ref_aa, ref_pos - pos_offset, inserted, None, pos_offset),
                           _ins_or_dup(ref_aa, ref_pos - pos_offset, inserted, C, pos_offset)))

        else:
            ref_pos += 1
            qry_seen += 1
            if q == "*" and r != "*":
                stops.append((f"p.{aa3(r)}{ref_pos}Ter", f"p.{aa3(r)}{C(ref_pos)}Ter"))
            elif r == "*" and q != "*":
                # read-through of the normal stop: HGVS extension
                tail = qry_aa[qry_seen - 1:]
                k = tail.find("*")
                ext = str(k + 1) if k != -1 else "?"
                stops.append((f"p.Ter{ref_pos}{aa3(q)}extTer{ext}",
                              f"p.Ter{C(ref_pos)}{aa3(q)}extTer{ext}"))
            elif r != q and q != "X" and r != "X":
                substitutions.append((f"p.{aa3(r)}{ref_pos}{aa3(q)}",
                                      f"p.{aa3(r)}{C(ref_pos)}{aa3(q)}"))
            i += 1

    matches = sum(1 for r, q in zip(ref_aln, qry_aln) if r == q and r != "-")
    aligned = sum(1 for r, q in zip(ref_aln, qry_aln) if r != "-" and q != "-")
    identity = 100.0 * matches / aligned if aligned else 0.0
    return substitutions, indels, stops, identity


def find_frame_disruption(ref_nt, qry_nt):
    """Locate the first frame-breaking indel, in reference nucleotide space.

    Must be given untrimmed sequences: rounding a sample CDS down to a whole
    number of codons is exactly what hides a single-base indel.

    Returns ``(ref_nt_pos, offset_before_break)`` where ``ref_nt_pos`` is the
    last reference base still read in frame and ``offset_before_break`` is the
    (in-frame) length difference accumulated ahead of it. Returns
    ``(None, net)`` when the reading frame is never disturbed.
    """
    runs = indel_runs(ref_nt, qry_nt)
    if not runs:
        return None, 0

    # A frameshift is a break the gene never recovers from. A divergent
    # reference can scatter compensating gaps through the alignment, and an
    # early run that is later cancelled out is not a frameshift -- judging a
    # run in isolation is what turns an in-frame 33 bp insertion into a
    # spurious "fsX360".
    if runs[-1]["cum_after"] % CODON == 0:
        return None, runs[-1]["cum_after"]

    for run in runs:
        if run["cum_after"] % CODON != 0:
            return run["ref_pos_before"], run["cum_before"]
    return None, runs[-1]["cum_after"]


def indel_runs(ref_nt, qry_nt, normalize=None):
    """Every indel between the two sequences, in reference coordinates.

    Each run records the reference base it follows, its size, its direction,
    and the running length offset before and after it.
    """
    normalize = normalize or _left_normalize
    ref_aln, qry_aln = normalize(*_aligned_pair(_NT_ALIGNER, ref_nt, qry_nt))
    out = []
    ref_i = 0
    cum = 0
    i = 0
    n = len(ref_aln)
    while i < n:
        if ref_aln[i] == "-":
            before, start, seq = cum, ref_i, ""
            while i < n and ref_aln[i] == "-":
                seq += qry_aln[i]
                cum += 1
                i += 1
            out.append({"kind": "ins", "ref_pos_before": start, "bases": seq,
                        "length": len(seq), "cum_before": before, "cum_after": cum})
        elif qry_aln[i] == "-":
            before, start, seq = cum, ref_i, ""
            while i < n and qry_aln[i] == "-" and ref_aln[i] != "-":
                seq += ref_aln[i]
                cum -= 1
                ref_i += 1
                i += 1
            out.append({"kind": "del", "ref_pos_before": start, "bases": seq,
                        "length": len(seq), "cum_before": before, "cum_after": cum})
        else:
            ref_i += 1
            i += 1
    return out


def hgvs_c(ref_nt, qry_nt, max_report=8, renumber=None, codons=None):
    """Describe the nucleotide differences in HGVS c. notation.

    Coordinates are positions in the reference CDS, indels are placed 3'-most,
    and an insertion copying the bases immediately 5' of it is written as a
    duplication -- the same rules as the protein description, so the two
    columns agree about what happened.

    ``codons`` restricts the substitutions reported to that set of 1-based
    reference codons -- the codons that changed the protein. Synonymous changes
    are left out, so this column lines up one-for-one with the protein
    description instead of listing silent differences alongside them. Indels are
    always reported: they are non-synonymous by definition.
    """
    ref_nt = trim_to_codon(ref_nt)

    def N(p):
        """Nucleotide position, mapped through the codon it belongs to."""
        if not renumber or p < 1:
            return p
        codon, within = (p - 1) // CODON + 1, (p - 1) % CODON
        mapped = renumber.get(codon)
        return p if mapped is None else (mapped - 1) * CODON + within + 1

    parts = []
    for r in indel_runs(ref_nt, qry_nt, normalize=_right_normalize):
        pos, n, bases = r["ref_pos_before"], r["length"], r["bases"]
        if r["kind"] == "ins":
            if pos >= n and ref_nt[pos - n:pos] == bases:
                first = pos - n + 1
                parts.append(f"c.{N(first)}dup" if n == 1 else f"c.{N(first)}_{N(pos)}dup")
            else:
                parts.append(f"c.{N(pos)}_{N(pos+1)}ins{bases}")
        else:
            first, last = pos + 1, pos + n
            parts.append(f"c.{N(first)}del" if n == 1 else f"c.{N(first)}_{N(last)}del")

    ref_aln, qry_aln = _right_normalize(*_aligned_pair(_NT_ALIGNER, ref_nt, qry_nt))
    subs, ref_i = [], 0
    for a, b in zip(ref_aln, qry_aln):
        if a != "-":
            ref_i += 1
        if a != "-" and b != "-" and a != b:
            if codons is not None and ((ref_i - 1) // CODON + 1) not in codons:
                continue                      # synonymous, or outside the reported set
            subs.append(f"c.{N(ref_i)}{a}>{b}")

    shown = parts + subs
    if not shown:
        return "c.="
    if len(shown) > max_report:
        return "; ".join(shown[:max_report]) + f"; +{len(shown) - max_report} more"
    return "; ".join(shown)


def _changed_codons(variants, pos_offset=0):
    """Reference codons named by a set of protein variant strings."""
    out = set()
    for v in variants:
        nums = [int(n) for n in re.findall(r"(\d+)", v)]
        if not nums:
            continue
        first = nums[0] - pos_offset
        last = (nums[1] - pos_offset) if (len(nums) > 1 and "_" in v.split("ins")[0]
                                          .split("del")[0].split("dup")[0]) else first
        out.update(range(min(first, last), max(first, last) + 1))
    return out


def analyze(ref_nt, qry_nt, expected_ref_aa_len=None, ref_nt_offset=0,
            canonical_nt=None):
    """Compare a sample CDS to its reference and describe the difference.

    ``ref_nt_offset`` is the reference nucleotide the supplied ``ref_nt``
    starts at, used to keep positions absolute when only part of the
    reference CDS was covered.

    ``canonical_nt`` is the panel's canonical allele for this gene. Calls are
    made, and numbered, against the allele that actually matched -- that is what
    stops ordinary lineage variation being reported as mutations. Where the
    canonical allele numbers a variant differently, that form is returned in
    ``VariantResult.canonical`` for reporting alongside.
    """
    ref_nt_raw, qry_nt_raw = ref_nt, qry_nt
    ref_nt = trim_to_codon(ref_nt)
    qry_nt = trim_to_codon(qry_nt)
    pos_offset = ref_nt_offset // CODON

    ref_aa = translate_cds(ref_nt).rstrip("*")
    qry_aa_raw = translate_cds(qry_nt)
    qry_aa = qry_aa_raw.rstrip("*")

    renumber = None
    if canonical_nt:
        canon_aa = translate_cds(trim_to_codon(canonical_nt)).rstrip("*")
        renumber = position_map(ref_aa, canon_aa)

    res = VariantResult(kind="intact", ref_aa_len=len(ref_aa), query_aa_len=len(qry_aa))
    if not ref_aa:
        res.kind = "translation_error"
        res.description = "Empty reference translation"
        return res

    # Any length difference, and any internal stop, could be a frameshift.
    # The check runs on the untrimmed sequences on purpose.
    net = len(qry_nt_raw) - len(ref_nt_raw)
    break_nt, offset_before = (None, 0)
    if net != 0 or "*" in qry_aa:
        break_nt, offset_before = find_frame_disruption(ref_nt_raw, qry_nt_raw)
    # Populated for every gene that differs from its reference, including
    # substitution-only genes -- the nucleotide change is as much a result as
    # the protein one. Overridden below for a premature stop.


    if break_nt is not None:
        qoff = offset_before // CODON
        idx = break_nt // CODON

        # HGVS anchors a frameshift to the first amino acid that actually
        # changes. That position is well defined even when the causal indel
        # falls inside a homopolymer, where the nucleotide position is not.
        while (idx < len(ref_aa) and idx + qoff < len(qry_aa_raw)
               and ref_aa[idx] == qry_aa_raw[idx + qoff]):
            idx += 1

        codon = idx + 1
        canon_codon = (renumber.get(codon, codon) or codon) if renumber else codon
        qidx = idx + qoff
        ref_res = ref_aa[idx] if 0 <= idx < len(ref_aa) else "?"
        new_res = qry_aa_raw[qidx] if 0 <= qidx < len(qry_aa_raw) else "?"
        stop_at = qry_aa_raw.find("*", max(qidx, 0))
        # HGVS: fsTer<n> gives the position of the new stop within the new
        # reading frame, counting the first changed amino acid as 1. Where the
        # new frame reaches no stop, HGVS omits the number entirely.
        if stop_at != -1:
            tail = f"fsTer{max(stop_at - qidx + 1, 1)}"
        else:
            tail = "fs"
        fs = f"p.{aa3(ref_res)}{codon + pos_offset}{aa3(new_res)}{tail}"
        fs_canon = f"p.{aa3(ref_res)}{canon_codon + pos_offset}{aa3(new_res)}{tail}"

        # Anything upstream of the break is still translated faithfully, so
        # report it: an isolate can carry the PBP3 duplication and a later
        # frameshift, and both matter.
        upstream, up_indels, up_stops, _ = call_protein_variants(
            ref_aa[:idx], qry_aa_raw[:max(qidx, 0)], pos_offset, renumber
        )
        res.kind = "frameshift"
        res.ref_pos = codon + pos_offset
        res.variants = [v[0] for v in up_indels + upstream + up_stops] + [fs]
        res.canonical = [v[1] for v in up_indels + upstream + up_stops] + [fs_canon]
        res.description = "; ".join(res.variants)
        if stop_at == -1:
            res.notes.append("no stop codon reached within aligned region")
        res.query_aa_len = stop_at if stop_at != -1 else len(qry_aa_raw)
        res.hgvs_c = hgvs_c(ref_nt_raw, qry_nt_raw,
                            codons=_changed_codons(res.variants, pos_offset))
        _, _, _, res.aa_identity = call_protein_variants(ref_aa, qry_aa, pos_offset)
        return res

    subs, indels, stops, identity = call_protein_variants(ref_aa, qry_aa, pos_offset, renumber)
    res.aa_identity = identity

    if stops:
        first = int(re.search(r"(\d+)", stops[0][0]).group(1))
        res.kind = "premature_stop"
        res.ref_pos = first
        # Substitutions after the stop are not translated in vivo.
        kept = [v for v in subs + indels
                if int((re.search(r"(\d+)", v[0].split("ins")[0].split("del")[0].split("dup")[0])
                        or _ZERO).group(1)) < first]
        res.variants = [v[0] for v in stops + kept]
        res.canonical = [v[1] for v in stops + kept]
        res.description = "; ".join(res.variants)
        res.query_aa_len = first - 1
        # Only the codon that became the stop. Everything downstream is not
        # translated, so listing its nucleotide changes is noise.
        res.hgvs_c = hgvs_c(ref_nt_raw, qry_nt_raw, codons={first - pos_offset})
        return res

    if expected_ref_aa_len and len(qry_aa) / expected_ref_aa_len < 0.9:
        res.kind = "truncation"
        res.ref_pos = len(qry_aa)
        res.description = f"Truncated protein: {len(qry_aa)}/{expected_ref_aa_len} AA"
        res.variants = [v[0] for v in subs + indels]
        res.canonical = [v[1] for v in subs + indels]
        res.hgvs_c = hgvs_c(ref_nt_raw, qry_nt_raw,
                            codons=_changed_codons(res.variants, pos_offset))
        return res

    res.variants = [v[0] for v in indels + subs]
    res.canonical = [v[1] for v in indels + subs]
    if indels and subs:
        res.kind = "inframe_indel_missense"
    elif indels:
        res.kind = "inframe_indel"
    elif subs:
        res.kind = "missense"
    else:
        res.kind = "intact"
        res.description = "Recovered intact: no AA diff"
        return res

    res.description = "; ".join(res.variants)
    res.hgvs_c = hgvs_c(ref_nt_raw, qry_nt_raw,
                        codons=_changed_codons(res.variants, pos_offset))
    return res
