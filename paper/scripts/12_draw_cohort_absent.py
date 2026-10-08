#!/usr/bin/env python3
"""Fifty permuted Absent pairs from the auto cohort's Existing/Established.

Takes genes and phenotypes that the new run called Existing or Established,
reassigns the phenotypes (same derangement as scripts/10_bench_shuffle.py),
and keeps a pair only when it also fails the figure-2 Absent tests: not a
real AoU/BRAVA association, not in GenCC, not a MONDO neighbour of a curated
disease of that gene, the phenotype names no gene, and BioLORD cosine below
0.60 against every disease GenCC records for the gene.

  python3 scripts/12_draw_cohort_absent.py
"""
from __future__ import annotations

import collections
import csv
import os
import random
import re
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("PESTO_PROJECT_ROOT", ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "efo"))

from pesto import ot_shortlist  # noqa: E402

from importlib import import_module  # noqa: E402
shuffle = import_module("10_bench_shuffle")

GENCC = os.path.join(ROOT, "data", "raw", "gencc-submissions.tsv")
from pesto.config import APP_ROOT as _PESTO  # noqa: E402
ALIASES = os.path.join(_PESTO, "aliases.tsv")
SRC = os.path.join(ROOT, "results", "all_runs_auto.tsv")
OUT = os.path.join(ROOT, "results", "bench_fig1c_absent50.tsv")
THRESHOLD = 0.60
N = 50
SEED = 2026
ENZYME = (
    "dehydrogenase", "synthase", "synthetase", "phosphatase", "bisphosphatase",
    "transferase", "reductase", "oxidase", "kinase", "lyase", "hydratase",
    "carboxylase", "acetyltransferase",
)


def norm(s):
    return re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip()


def load_gene_names():
    gencc = pd.read_csv(GENCC, sep="\t", dtype=str, low_memory=False)
    symbols = set(gencc.gene_symbol.dropna().astype(str))
    approved_name, alias_of = {}, {}
    gene_aliases = collections.defaultdict(set)
    with open(ALIASES, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            approved = (row.get("Approved symbol") or "").strip()
            if not approved:
                continue
            approved_name[approved.upper()] = (row.get("Approved name") or "").strip()
            alias_of[approved.upper()] = approved.upper()
            for cell in (row.get("Previous symbols") or "",
                         row.get("Alias symbols") or ""):
                for tok in re.split(r"[,;]", cell):
                    tok = tok.strip()
                    if tok:
                        alias_of.setdefault(tok.upper(), approved.upper())
                        gene_aliases[approved.upper()].add(tok)
    return symbols, approved_name, alias_of, gene_aliases


def names_a_gene(phenotype, gene, symbols, approved_name, alias_of, gene_aliases):
    """Whether the label carries a gene symbol, its own or anyone's.

    Copied from scripts/draw_figure2_benchmark.py so this draw does not import
    that script (it parses argv and runs on import).
    """
    pheno_n = norm(phenotype)
    if "deficien" in pheno_n and any(e in pheno_n for e in ENZYME):
        return True
    tokens = re.split(r"[^A-Za-z0-9]+", str(phenotype))

    def token_is_a_gene(t):
        if t.islower():
            return False
        if t in symbols:
            return True
        return t in alias_of and (any(c.isdigit() for c in t) or len(t) >= 4)

    if any(token_is_a_gene(t) for t in tokens):
        return True
    if gene:
        own = {gene} | gene_aliases.get(gene.upper(), set())
        if any(t in own and not t.islower() for t in tokens):
            return True
    name = approved_name.get((gene or "").upper(), "")
    if not name:
        return False
    stop = {"gene", "protein", "the", "of", "and", "a"}
    need = [w for w in norm(name).split() if w not in stop and not w.isdigit()]
    have = [w for w in norm(phenotype).split() if w not in stop]
    if not need:
        return False

    def hit(word):
        if word in have:
            return True
        if len(word) < 6:
            return False
        return any(h.startswith(word[:6]) and word.startswith(h[:6])
                   for h in have if len(h) >= 6)

    return all(hit(w) for w in need)


def main():
    pool = []
    with open(SRC, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["source"] in ("AoU", "BRAVA") and r["verdict"] in (
                    "Existing", "Established"):
                pool.append((r["gene"], r["phenotype"]))
    if len(pool) < N:
        raise SystemExit(f"only {len(pool)} Existing/Established pairs")

    gencc = pd.read_csv(GENCC, sep="\t", dtype=str, low_memory=False)
    gencc = gencc.dropna(subset=["gene_symbol", "disease_title", "disease_curie"])

    curated = collections.defaultdict(set)
    curies_of_gene = collections.defaultdict(set)
    title_curie = {}
    for gene, title, original, curie in zip(
            gencc.gene_symbol, gencc.disease_title,
            gencc.get("disease_original_title", [""] * len(gencc)),
            gencc.disease_curie):
        curated[norm(gene)].add(str(title))
        if isinstance(original, str) and original:
            curated[norm(gene)].add(original)
        cid = str(curie).replace(":", "_")
        curies_of_gene[norm(gene)].add(cid)
        title_curie[norm(title)] = cid
        if isinstance(original, str) and original:
            title_curie.setdefault(norm(original), cid)
    known = {(g, norm(d)) for g, ds in curated.items() for d in ds}

    import relate  # noqa: E402
    graph = relate.Graph()
    RELATED = (relate.SAME, relate.DESCENDANT, relate.ANCESTOR)
    for tid, term in graph.terms.items():
        label = term.get("label") or ""
        if label:
            title_curie.setdefault(norm(label), tid)

    def canon(curie):
        return graph.canonical(str(curie).replace(":", "_"))

    symbols, approved_name, alias_of, gene_aliases = load_gene_names()
    enc = ot_shortlist.encoder("biolord")
    if enc is None:
        raise SystemExit("BioLORD encoder required for the cosine cut")

    real = shuffle.real_pairs()
    real_n = {(g, norm(p)) for g, p in real}
    with open(SRC, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            real_n.add((r["gene"], norm(r["phenotype"])))
    rng = random.Random(SEED)
    kept, seen_g, seen_p = [], set(), set()
    attempts = 0
    refused = collections.Counter()
    while len(kept) < N and attempts < 80:
        attempts += 1
        by_gene = {}
        for g, p in pool:
            if g not in seen_g:
                by_gene.setdefault(g, []).append(p)
        need = N - len(kept)
        if len(by_gene) < need:
            break
        genes = rng.sample(sorted(by_gene), min(len(by_gene), max(need * 4, need)))
        draw = [(g, rng.choice(by_gene[g])) for g in genes]
        gs = [g for g, _ in draw]
        ps = [p for _, p in draw]
        real_exact = set(real)
        for g, p in real_n:
            for s in ps:
                if norm(s) == p:
                    real_exact.add((g, s))
        try:
            order = shuffle.derange(gs, ps, real_exact, rng)
        except SystemExit:
            continue
        cands = [(gs[i], ps[order[i]], ps[i]) for i in range(len(gs))]
        texts = sorted({p for _, p, _ in cands}
                       | {d for g, _, _ in cands for d in curated[norm(g)]})
        vec = enc.encode(texts, max_length=64)
        idx = {t: i for i, t in enumerate(texts)}

        def nearest(gene, pheno):
            sims = []
            for d in curated[norm(gene)]:
                if pheno in idx and d in idx:
                    sims.append((float(vec[idx[pheno]] @ vec[idx[d]]), d))
            return max(sims) if sims else (0.0, "")

        for gene, pheno, real_pheno in cands:
            if len(kept) == N:
                break
            if gene in seen_g or pheno in seen_p:
                continue
            if (gene, pheno) in real or (gene, norm(pheno)) in real_n:
                refused["real_cohort"] += 1
                continue
            if (norm(gene), norm(pheno)) in known:
                refused["gencc_exact"] += 1
                continue
            if names_a_gene(pheno, gene, symbols, approved_name,
                            alias_of, gene_aliases):
                refused["names_a_gene"] += 1
                continue
            qid = title_curie.get(norm(pheno))
            if qid:
                related = False
                for c in sorted(curies_of_gene[norm(gene)]):
                    if graph.relate(canon(qid), canon(c)) in RELATED:
                        related = True
                        break
                if related:
                    refused["mondo"] += 1
                    continue
            cos, near = nearest(gene, pheno)
            if cos >= THRESHOLD:
                refused["cosine"] += 1
                continue
            kept.append({
                "gene": gene,
                "phenotype": pheno,
                "phenotype_real": real_pheno,
                "gencc_class": "Absent",
                "max_cosine": f"{cos:.3f}",
                "nearest_curated_disease": near,
            })
            seen_g.add(gene)
            seen_p.add(pheno)

    if len(kept) < N:
        raise SystemExit(
            f"only {len(kept)} Absent pairs after {attempts} draws "
            f"(refused {dict(refused)})")

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(kept[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(kept[:N])
    print(f"wrote {N} pairs to {os.path.relpath(OUT, ROOT)} "
          f"({attempts} derangements)")
    print("refused", dict(refused))
    print("nearest cosine",
          min(float(r["max_cosine"]) for r in kept),
          "–", max(float(r["max_cosine"]) for r in kept))


if __name__ == "__main__":
    main()
