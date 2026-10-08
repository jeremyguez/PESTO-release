"""Thirty more random Absents, same rules as Figure 2a's 30 and the extra 40.

Skips genes and phenotypes already in the Figure 2a table or in the two
extra-20 draws. Seed 2029 so the draw is independent of 2027 and 2028.

Usage:
  python3 scripts/100_draw_fig2a_absent30.py
"""
from __future__ import annotations

import collections
import csv
import os
import random
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("PESTO_PROJECT_ROOT", str(ROOT))
sys.path.insert(0, str(ROOT / "efo"))

GENCC = ROOT / "data" / "raw" / "gencc-submissions.tsv"
BENCH = ROOT / "benchmark" / "gencc_clingen_g2p" / "pesto.tsv"
ALIASES = ROOT / "src" / "pesto" / "aliases.tsv"
SKIP = (
    ROOT / "results" / "fig2a_absent_extra20.tsv",
    ROOT / "results" / "fig2a_absent_extra20b.tsv",
)
OUT = ROOT / "results" / "fig2a_absent_extra30.tsv"
N = 30
THRESHOLD = 0.60
MAX_DESCENDANTS = 100
SEED = 2029
ENZYME = (
    "dehydrogenase", "synthase", "synthetase", "phosphatase", "bisphosphatase",
    "transferase", "reductase", "oxidase", "kinase", "lyase", "hydratase",
    "carboxylase", "acetyltransferase",
)
HEADINGS = {"congenital nervous system disorder", "thrombotic disease"}


def norm(s):
    return re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip()


def main():
    rng = random.Random(SEED)
    gencc = pd.read_csv(GENCC, sep="\t", dtype=str, low_memory=False)
    gencc = gencc.dropna(subset=["gene_symbol", "disease_title",
                                 "classification_title", "disease_curie"])

    classes_of = collections.defaultdict(set)
    titles_of = collections.defaultdict(collections.Counter)
    for gene, title, curie, cls in zip(
            gencc.gene_symbol, gencc.disease_title, gencc.disease_curie,
            gencc.classification_title):
        key = (gene, curie)
        classes_of[key].add(cls)
        titles_of[key][str(title)] += 1

    def title_of(key):
        return min(titles_of[key].items(), key=lambda kv: (-kv[1], kv[0]))[0]

    SYMBOLS = set(gencc.gene_symbol.dropna().astype(str))
    approved_name = {}
    alias_of = {}
    gene_aliases = collections.defaultdict(set)
    with ALIASES.open(encoding="utf-8") as fh:
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

    def names_a_gene(phenotype, gene=""):
        pheno_n = norm(phenotype)
        if "deficien" in pheno_n and any(e in pheno_n for e in ENZYME):
            return True
        tokens = re.split(r"[^A-Za-z0-9]+", str(phenotype))

        def token_is_a_gene(t):
            if t.islower():
                return False
            if t in SYMBOLS:
                return True
            return t in alias_of and (any(c.isdigit() for c in t) or len(t) >= 4)

        if any(token_is_a_gene(t) for t in tokens):
            return True
        if gene:
            own = {gene} | gene_aliases.get(gene.upper(), set())
            if any(t in own and not t.islower() for t in tokens):
                return True
        if not gene:
            return False
        name = approved_name.get(gene.upper(), "")
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

    from pesto import ot_shortlist
    import relate

    enc = ot_shortlist.encoder("biolord")
    if enc is None:
        sys.exit("biolord encoder missing")
    graph = relate.Graph()
    RELATED = (relate.SAME, relate.DESCENDANT, relate.ANCESTOR)

    def canon(curie):
        return graph.canonical(str(curie).replace(":", "_"))

    def descendants(curie):
        tid = canon(curie)
        return len(graph.descendants(tid)) if tid in graph.terms else 0

    def is_category(curie):
        return descendants(curie) >= MAX_DESCENDANTS

    curated = collections.defaultdict(set)
    for col in ("disease_title", "disease_original_title"):
        for gene, disease in zip(gencc.gene_symbol, gencc[col]):
            if isinstance(disease, str):
                curated[norm(gene)].add(disease)
    known = {(g, norm(d)) for g, ds in curated.items() for d in ds}
    curies_of_gene = collections.defaultdict(set)
    for gene, curie in zip(gencc.gene_symbol, gencc.disease_curie):
        curies_of_gene[norm(gene)].add(canon(curie))

    spent, seen_p = set(), set()
    with BENCH.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            spent.add(r["gene"])
            if r.get("phenotype"):
                seen_p.add(r["phenotype"])
    for path in SKIP:
        with path.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                spent.add(r["gene"])
                if r.get("phenotype"):
                    seen_p.add(r["phenotype"])

    gene_pool = sorted(SYMBOLS)
    pheno_pool = sorted({title_of(k): k[1] for k in classes_of
                         if not is_category(k[1])
                         and title_of(k) not in HEADINGS
                         and not names_a_gene(title_of(k))}.items())
    print(f"pool {len(gene_pool)} genes, {len(pheno_pool)} phenotypes; "
          f"spent {len(spent)} genes, {len(seen_p)} phenotypes from Fig 2a + extra 40")

    candidates, seen_g = [], set(spent)
    while len(candidates) < N * 8 and len(seen_p) < len(pheno_pool):
        gene = rng.choice(gene_pool)
        pheno, curie = rng.choice(pheno_pool)
        if gene in seen_g or pheno in seen_p:
            continue
        if (norm(gene), norm(pheno)) in known or (gene, curie) in classes_of:
            continue
        if any(graph.relate(canon(curie), c) in RELATED
               for c in sorted(curies_of_gene[norm(gene)])):
            continue
        candidates.append((gene, pheno, curie))
        seen_g.add(gene)
        seen_p.add(pheno)

    need = [(g, p) for g, p, _ in candidates]
    ctexts = sorted({p for _, p in need}
                    | {d for g, _ in need for d in curated[norm(g)]})
    print(f"{len(candidates)} candidates, {len(ctexts)} labels to embed")
    cvec = enc.encode(ctexts, max_length=64)
    cidx = {t: i for i, t in enumerate(ctexts)}

    def nearest(gene, pheno):
        sims = [(float(cvec[cidx[pheno]] @ cvec[cidx[d]]), d)
                for d in sorted(curated[norm(gene)]) if d in cidx]
        return max(sims) if sims else (0.0, "")

    absent, refused = [], 0
    for gene, pheno, _ in candidates:
        if len(absent) == N:
            break
        if nearest(gene, pheno)[0] >= THRESHOLD:
            refused += 1
            continue
        absent.append((gene, pheno))
    if len(absent) < N:
        sys.exit(f"only {len(absent)} controls survive at cosine {THRESHOLD}")
    print(f"drew {len(absent)}, refused {refused} on cosine")

    fields = ["gene", "phenotype", "gencc_class", "max_cosine",
              "nearest_curated_disease"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t")
        w.writeheader()
        for gene, pheno in absent:
            cos, d = nearest(gene, pheno)
            w.writerow({"gene": gene, "phenotype": pheno, "gencc_class": "Absent",
                        "max_cosine": f"{cos:.3f}", "nearest_curated_disease": d})
            print(f"  {gene:12} {pheno[:52]:52}  cos {cos:.2f}  {d[:40]}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
