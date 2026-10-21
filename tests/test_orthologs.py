"""Ortholog route of misfit-db (no network): the reciprocal-best-hit check must
pick the true ortholog and refuse a paralog that is merely the closest homolog."""
import random

from misfit.scripts.manage_refs import find_ortholog

AA = "ACDEFGHIKLMNPQRSTVWY"
rnd = random.Random(7)


def mutate(p, frac):
    p = list(p)
    for i in rnd.sample(range(len(p)), int(len(p) * frac)):
        p[i] = rnd.choice([a for a in AA if a != p[i]])
    return "".join(p)


def cds(name, aa):
    return {"gene": name, "protein_id": f"WP_{name}", "nuccore": "NZ_X", "locus_tag": name,
            "product": name, "seq": "ATG" * len(aa), "aa": aa}


X = "M" + "".join(rnd.choice(AA) for _ in range(400))      # the gene
Y = mutate(X, 0.55)                                         # a paralog (~45% identical)
seed_proteome = [cds("geneX", X), cds("geneY", Y)] + \
    [cds(f"o{i}", "M" + "".join(rnd.choice(AA) for _ in range(300))) for i in range(30)]
seed = {"panel": "e_coli", "ref_id": "geneX|WP_geneX|NZ_X|geneX", "aa": X, "seq": "",
        "assembly": "GCF_seed"}


def unnamed(c):
    return {**c, "gene": ""}


def test_picks_true_ortholog_over_paralog():
    target = [unnamed(cds("t1", mutate(X, 0.25))), unnamed(cds("t2", mutate(Y, 0.25)))]
    cand, ev = find_ortholog("geneX", target, seed, seed_proteome, min_id=60)
    assert cand is not None and cand["protein_id"] == "WP_t1", ev
    assert "reciprocal best hit" in ev


def test_refuses_paralog_when_ortholog_missing():
    target = [unnamed(cds("t2", mutate(Y, 0.15)))]   # only the paralog's ortholog present
    cand, ev = find_ortholog("geneX", target, seed, seed_proteome, min_id=30)
    assert cand is None and "reciprocal" in ev, ev


def test_identity_floor():
    target = [unnamed(cds("t1", mutate(X, 0.5)))]
    cand, ev = find_ortholog("geneX", target, seed, seed_proteome, min_id=60)
    assert cand is None and "identity" in ev, ev


def test_refseq_name_wins():
    target = [cds("geneX", mutate(X, 0.3)), unnamed(cds("t1", mutate(X, 0.05)))]
    cand, ev = find_ortholog("geneX", target, seed, seed_proteome)
    assert cand["gene"] == "geneX" and ev == "refseq_gene_name"


def test_long_protein_does_not_win():
    """A much longer protein must not beat the true ortholog through gappy alignment."""
    long = "".join(rnd.choice(AA) for _ in range(3400))
    target = [unnamed(cds("big", long)), unnamed(cds("t1", mutate(X, 0.25)))]
    cand, ev = find_ortholog("geneX", target, seed, seed_proteome, min_id=60)
    assert cand is not None and cand["protein_id"] == "WP_t1", ev


# ---- synteny route, alternative starts, panel collisions ------------------------------
from misfit.scripts.manage_refs import (find_neighbor, check_orf_annotated, _collision,
                                        _write_gene)

A = "M" + "".join(rnd.choice(AA) for _ in range(110))      # anchor (ramA-like)
R = "M" + "".join(rnd.choice(AA) for _ in range(190))      # neighbour (ramR-like)
anchor = {"panel": "k_pneumoniae", "ref_id": "ramA|WP_A|NZ_K|A", "aa": A, "seq": "",
          "assembly": "GCF_k"}


def gcds(name, aa, product, nuc="NZ_T"):
    return {**cds(name, aa), "gene": "", "product": product, "nuccore": nuc}


def genome(neigh_product="TetR/AcrR family transcriptional regulator", gap=1):
    filler = [gcds(f"f{i}", "M" + "".join(rnd.choice(AA) for _ in range(150)),
                   "hypothetical protein") for i in range(6)]
    return (filler[:3] + [gcds("ramA_t", mutate(A, 0.1), "AraC family regulator")]
            + filler[3:3 + gap - 1]
            + [gcds("ramR_t", mutate(R, 0.4), neigh_product)] + filler[3 + gap - 1:])


def test_neighbor_found_next_to_anchor():
    c, ev = find_neighbor("ramR", genome(), anchor, "TetR|RamR")
    assert c is not None and c["locus_tag"] == "ramR_t", ev
    assert "neighbour of ramA" in ev and "1 gene(s) away" in ev


def test_neighbor_respects_window_and_product():
    c, ev = find_neighbor("ramR", genome(gap=4), anchor, "TetR|RamR", window=2)
    assert c is None and "within 2 genes" in ev
    c, ev = find_neighbor("ramR", genome("LysR family regulator"), anchor, "TetR|RamR")
    assert c is None


def test_neighbor_needs_anchor():
    g = [x for x in genome() if x["locus_tag"] != "ramA_t"]
    c, ev = find_neighbor("ramR", g, anchor, "TetR|RamR")
    assert c is None and "anchor ramA not found" in ev


def test_find_ortholog_uses_synteny_without_seed():
    c, ev = find_ortholog("ramR", genome(), None, None, anchor=anchor, product_regex="TetR")
    assert c["locus_tag"] == "ramR_t" and ev.startswith("neighbour")


def test_alternative_start_accepted():
    nt = "ATT" + "GCT" * 50 + "TAA"
    assert check_orf_annotated(nt) == []
    assert check_orf_annotated("AAA" + "GCT" * 50 + "TAA")      # not a start in table 11


def test_short_protein_is_not_a_collision(tmp_path):
    """arnE (111 aa) used to 'collide' with emrA (390 aa) under the gap-free identity."""
    from misfit.scripts.manage_refs import Seq
    big = "M" + "".join(rnd.choice(AA) for _ in range(390))
    nt = "ATG" + "".join({"A": "GCT"}.get(a, "GCT") for a in big[1:]) + "TAA"
    (tmp_path / "ref_emrA.fasta").write_text(f">emrA|WP|NZ|x\n{nt}\n")
    small = str(Seq(nt[:333]).translate(table="Bacterial"))
    assert _collision(tmp_path, small, "arnE") is None
    assert _collision(tmp_path, str(Seq(nt).translate(table="Bacterial")).rstrip("*"),
                      "emrX") is not None


def test_allele_identity_threshold(tmp_path):
    pdir = tmp_path / "p"; pdir.mkdir()
    nt = "ATG" + "".join(rnd.choice(["GCT", "AAA", "CTG", "GGC", "TTC"]) for _ in range(200)) + "TAA"
    (pdir / "ref_g.fasta").write_text(f">g|WP_1|NZ_1|l1\n{nt}\n")
    (tmp_path / "PROVENANCE.tsv").write_text("")
    from misfit.scripts.manage_refs import translate
    aa = translate(nt)
    mut = mutate(aa, 0.27)                     # ~73% identical
    back = {"A": "GCT", "K": "AAA", "L": "CTG", "G": "GGC", "F": "TTC"}
    mnt = "ATG" + "".join(back.get(a, "GCT") for a in mut[1:]) + "TAA"
    cand = {"protein_id": "WP_2", "nuccore": "NZ_2", "locus_tag": "l2", "product": "x",
            "seq": mnt, "aa": translate(mnt)}
    assert not _write_gene(tmp_path, "p", "g", cand, "GCF_2", "o", "ev", True, False)
    assert _write_gene(tmp_path, "p", "g", cand, "GCF_2", "o", "ev", True, False,
                       allele_min_id=60)
