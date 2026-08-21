"""Ground-truth tests for reference-anchored variant calling.

Queries are built by editing the real ftsI CDS at known codons, so the
expected call is known exactly rather than eyeballed from an alignment.
"""
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from misfit.scripts.variants import analyze, AA3

REF = "/Users/wcshropshire/Documents/GitHub/misfit/src/misfit/reference/e_coli/ref_ftsI.fasta"
# Any full-length ftsI allele works; the tests only rely on the canonical
# 588 aa numbering, which the assert below pins down.
rec = next(r for r in SeqIO.parse(REF, "fasta")
           if len(r.seq) // 3 - 1 == 588)
CDS = str(rec.seq).upper()
AA = str(Seq(CDS).translate(table="Bacterial")).rstrip("*")

CODONS = [CDS[i:i + 3] for i in range(0, len(CDS), 3)]
assert AA[332] == "P" and AA[412] == "A" and AA[336] == "N", "reference landmarks moved"


def build(edits=(), insert_after=None, insert_codons=(), delete=None, del_bases=None):
    """edits: list of (1-based codon, new AA). insert_after: 1-based codon."""
    cods = list(CODONS)
    for pos, new_aa in edits:
        cods[pos - 1] = {
            "V": "GTT", "T": "ACC", "Y": "TAT", "F": "TTT", "K": "AAA",
            "L": "CTG", "S": "AGC", "*": "TAA", "H": "CAT", "W": "TGG",
        }[new_aa]
    if delete:
        s, e = delete
        cods = cods[: s - 1] + cods[e:]
    if insert_after is not None:
        cods = cods[:insert_after] + list(insert_codons) + cods[insert_after:]
    seq = "".join(cods)
    if del_bases:
        pos, n = del_bases
        seq = seq[: pos - 1] + seq[pos - 1 + n:]
    return seq


def check(name, qry, expect_kind, expect_desc):
    r = analyze(CDS, qry, expected_ref_aa_len=len(AA))
    got = r.description
    ok = (r.kind == expect_kind) and (got == expect_desc)
    print(f"{'PASS' if ok else 'FAIL'}  {name}")
    if not ok:
        print(f"        expected [{expect_kind}] {expect_desc}")
        print(f"        got      [{r.kind}] {got}")
    return ok


results = []
YRIN = CODONS[333:337]          # codons 334-337 = Y R I N

# 1. clean intact
results.append(check("intact", CDS, "intact", "Recovered intact: no AA diff"))

# 2. single missense, no indel anywhere
results.append(check("missense only (A413V)", build(edits=[(413, "V")]),
                     "missense", "p.Ala413Val"))

# 3. the regression: YRIN insertion + a downstream substitution.
#    Alignment-column numbering would call the substitution A417V.
q = build(edits=[(413, "V")], insert_after=333, insert_codons=YRIN)
results.append(check("YRIN insertion + A413V (was A417V)", q,
                     "inframe_indel_missense", "p.Tyr334_Asn337dup; p.Ala413Val"))

# 4. YRIK insertion + A413V, matching isolate EC1092
YRIK = YRIN[:3] + ["AAA"]
q = build(edits=[(413, "V")], insert_after=333, insert_codons=YRIK)
results.append(check("YRIK insertion + A413V", q,
                     "inframe_indel_missense", "p.Ile336_Asn337insLysTyrArgIle; p.Ala413Val"))

# 5. insertion with several downstream substitutions
q = build(edits=[(317, "T"), (413, "V"), (498, "T")], insert_after=333, insert_codons=YRIN)
results.append(check("YRIN + three substitutions spanning the insertion", q,
                     "inframe_indel_missense", "p.Tyr334_Asn337dup; p.Ala317Thr; p.Ala413Val; p.Ala498Thr"))

# 6. the MDA allele, no indel
results.append(check("MDA allele A317T/L369F/S408Y/A498T",
                     build(edits=[(317, "T"), (369, "F"), (408, "Y"), (498, "T")]),
                     "missense", "p.Ala317Thr; p.Leu369Phe; p.Ser408Tyr; p.Ala498Thr"))

# 7. in-frame deletion, and a substitution downstream of it.
#    A deletion must NOT shift reference numbering either.
q = build(edits=[(413, "V")], delete=(200, 202))
results.append(check("3-codon deletion + downstream A413V", q,
                     "inframe_indel_missense",
                     f"p.{AA3[AA[199]]}200_{AA3[AA[201]]}202del; p.Ala413Val"))

# 8. frameshift from a 1-base deletion at codon 29 (acrR-like)
q = build(del_bases=(85, 1))
r = analyze(CDS, q, expected_ref_aa_len=len(AA))
ok = r.kind == "frameshift" and r.ref_pos == 29
print(f"{'PASS' if ok else 'FAIL'}  frameshift at codon 29 -> got [{r.kind}] pos={r.ref_pos} {r.description}")
results.append(ok)

# 9. frameshift downstream of an in-frame insertion: position must stay
#    in reference coordinates (codon 400, not 404)
q = build(insert_after=333, insert_codons=YRIN)
q = q[: (333 + 4) * 3 + (400 - 334) * 3] + q[(333 + 4) * 3 + (400 - 334) * 3 + 1:]
r = analyze(CDS, q, expected_ref_aa_len=len(AA))
ok = r.kind == "frameshift" and r.ref_pos == 400
print(f"{'PASS' if ok else 'FAIL'}  frameshift at codon 400 downstream of YRIN -> "
      f"got [{r.kind}] pos={r.ref_pos} {r.description}")
results.append(ok)

# 10. premature stop, in reference coordinates, downstream of an insertion
q = build(edits=[(450, "*")], insert_after=333, insert_codons=YRIN)
r = analyze(CDS, q, expected_ref_aa_len=len(AA))
ok = r.kind == "premature_stop" and r.ref_pos == 450
print(f"{'PASS' if ok else 'FAIL'}  premature stop at 450 downstream of YRIN -> "
      f"got [{r.kind}] pos={r.ref_pos} {r.description}")
results.append(ok)

# 11. determinism: the same event must get the same name every run
q = build(insert_after=333, insert_codons=YRIN)
names = {analyze(CDS, q, expected_ref_aa_len=len(AA)).description for _ in range(5)}
ok = len(names) == 1
print(f"{'PASS' if ok else 'FAIL'}  deterministic naming -> {names}")
results.append(ok)

# 12. ambiguous bases must not silently truncate the protein
q = CDS[:600] + "NNN" + CDS[603:]
r = analyze(CDS, q, expected_ref_aa_len=len(AA))
ok = r.kind in ("missense", "intact") and r.query_aa_len == len(AA)
print(f"{'PASS' if ok else 'FAIL'}  ambiguous codon kept in frame -> "
      f"[{r.kind}] aa_len={r.query_aa_len} (ref {len(AA)}) {r.description}")
results.append(ok)

print(f"\n{sum(results)}/{len(results)} passed")
sys.exit(0 if all(results) else 1)
