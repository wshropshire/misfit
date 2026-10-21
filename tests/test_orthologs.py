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
