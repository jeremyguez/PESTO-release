#!/usr/bin/env python3
"""Fifty Absent pairs: ordinary HGNC gene × random AoU/BRAVA trait.

The gene is an HGNC approved symbol that looks like a name people use
(protein-coding, no LINC, orf, olfactory, immunoglobulin or tRNA loci).
The phenotype is drawn from every trait in the two cohort tables. A pair
is kept only when it fails the figure-2 Absent tests.

  python3 scripts/13_draw_random_absent.py
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
from pesto.config import AOU_FILE, BRAVA_FILE  # noqa: E402

from importlib import import_module  # noqa: E402
draw12 = import_module("12_draw_cohort_absent")


def real_pairs():
    """Every association that genuinely exists in either cohort."""
    real = set()
    with open(AOU_FILE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            real.add((r["gene_symbol"], r["description"]))
    with open(BRAVA_FILE) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            real.add((r["external_gene_name"], r["phenotype_full"]))
            real.add((r["external_gene_name"], r["phenotype"]))
    return real

GENCC = os.path.join(ROOT, "data", "raw", "gencc-submissions.tsv")
from pesto.config import APP_ROOT as _PESTO  # noqa: E402
ALIASES = os.path.join(_PESTO, "aliases.tsv")
OUT = os.path.join(ROOT, "results", "bench_fig1c_absent50_hgnc.tsv")
THRESHOLD = 0.60
N = 50
SEED = 2026
OK_SYM = re.compile(r"^[A-Z][A-Z0-9]{1,9}(?:-[A-Z0-9]{1,4})?$")
WEIRD_SYM = re.compile(
    r"^(?:OR\d|KIR\d|LILR|IGH|IGK|IGL|TRA[VDJC]|TRB|TRG|TRD|FAM\d|KIAA\d|"
    r"LOC\d|C\d+ORF|LINC|MIR|RNU|SNOR|SNAR|PIRC|ERV|MT-)|"
    r"(?:-AS\d*|-IT\d*|-DT)$|orf\d",
    re.I)
WEIRD_NAME = (
    "antisense", "long intergenic", "long non-coding", "long noncoding",
    "microrna", "small nuclear", "small nucleolar", "pseudogene",
    "readthrough", "uncharacterized", "non-protein coding",
    "olfactory receptor", "putative", "transfer rna", "trna ",
    "pirna", "small nf90", "immunoglobulin",
)


def ordinary_symbols():
    """HGNC approved symbols that read as ordinary protein-coding genes."""
    out = []
    with open(ALIASES, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sym = (row.get("Approved symbol") or "").strip()
            name = (row.get("Approved name") or "").strip()
            if not sym or not OK_SYM.match(sym) or WEIRD_SYM.search(sym):
                continue
            if sym.count("-") >= 2:
                continue
            low = name.lower()
            if any(w in low for w in WEIRD_NAME):
                continue
            out.append(sym)
    return sorted(set(out))


def cohort_phenotypes():
    """One wording per trait, preferring the longer AoU/BRAVA label."""
    by_norm = {}
    with open(AOU_FILE, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            p = (r.get("description") or "").strip()
            if p:
                by_norm.setdefault(draw12.norm(p), p)
                if len(p) > len(by_norm[draw12.norm(p)]):
                    by_norm[draw12.norm(p)] = p
    with open(BRAVA_FILE, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for key in ("phenotype_full", "phenotype"):
                p = (r.get(key) or "").strip()
                if not p:
                    continue
                k = draw12.norm(p)
                if k not in by_norm or len(p) > len(by_norm[k]):
                    by_norm[k] = p
    return sorted(by_norm.values())


def main():
    genes = ordinary_symbols()
    phenos = cohort_phenotypes()
    if len(genes) < N or len(phenos) < N:
        raise SystemExit(f"pool too small: {len(genes)} genes, {len(phenos)} traits")

    gencc = pd.read_csv(GENCC, sep="\t", dtype=str, low_memory=False)
    gencc = gencc.dropna(subset=["gene_symbol", "disease_title", "disease_curie"])
    curated = collections.defaultdict(set)
    curies_of_gene = collections.defaultdict(set)
    title_curie = {}
    for gene, title, original, curie in zip(
            gencc.gene_symbol, gencc.disease_title,
            gencc.get("disease_original_title", [""] * len(gencc)),
            gencc.disease_curie):
        curated[draw12.norm(gene)].add(str(title))
        if isinstance(original, str) and original:
            curated[draw12.norm(gene)].add(original)
        cid = str(curie).replace(":", "_")
        curies_of_gene[draw12.norm(gene)].add(cid)
        title_curie[draw12.norm(title)] = cid
        if isinstance(original, str) and original:
            title_curie.setdefault(draw12.norm(original), cid)

    known = {(g, draw12.norm(d)) for g, ds in curated.items() for d in ds}
    import relate  # noqa: E402
    graph = relate.Graph()
    RELATED = (relate.SAME, relate.DESCENDANT, relate.ANCESTOR)
    # Only GenCC titles get a curie, as in figure 2. Indexing every EFO
    # label would map Height/weight onto the graph and force a walk for
    # every curated gene.

    def canon(curie):
        return graph.canonical(str(curie).replace(":", "_"))

    related_memo = {}

    def is_related(qid, cid):
        key = (qid, cid)
        if key not in related_memo:
            related_memo[key] = graph.relate(canon(qid), canon(cid)) in RELATED
        return related_memo[key]

    symbols, approved_name, alias_of, gene_aliases = draw12.load_gene_names()
    print(f"pool {len(genes)} ordinary genes × {len(phenos)} traits", flush=True)
    enc = ot_shortlist.encoder("biolord")
    if enc is None:
        raise SystemExit("BioLORD encoder required for the cosine cut")

    real = real_pairs()
    real_n = {(g, draw12.norm(p)) for g, p in real}
    auto = os.path.join(ROOT, "results", "all_runs_auto.tsv")
    if os.path.exists(auto):
        with open(auto, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                real_n.add((r["gene"], draw12.norm(r["phenotype"])))

    rng = random.Random(SEED)
    seen_g, seen_p = set(), set()
    cheap, refused = [], collections.Counter()
    tries = 0
    while len(cheap) < N * 8 and tries < 20000:
        tries += 1
        gene = rng.choice(genes)
        pheno = rng.choice(phenos)
        if gene in seen_g or pheno in seen_p:
            continue
        if (gene, pheno) in real or (gene, draw12.norm(pheno)) in real_n:
            refused["real_cohort"] += 1
            continue
        if (draw12.norm(gene), draw12.norm(pheno)) in known:
            refused["gencc_exact"] += 1
            continue
        if draw12.names_a_gene(pheno, gene, symbols, approved_name,
                               alias_of, gene_aliases):
            refused["names_a_gene"] += 1
            continue
        qid = title_curie.get(draw12.norm(pheno))
        if qid and any(is_related(qid, c)
                       for c in curies_of_gene[draw12.norm(gene)]):
            refused["mondo"] += 1
            continue
        seen_g.add(gene)
        seen_p.add(pheno)
        cheap.append((gene, pheno))
    print(f"{len(cheap)} candidates after {tries} draws, {dict(refused)}",
          flush=True)

    texts = sorted({p for _, p in cheap}
                   | {d for g, _ in cheap for d in curated[draw12.norm(g)]})
    vec = enc.encode(texts, max_length=64)
    idx = {t: i for i, t in enumerate(texts)}

    def nearest(gene, pheno):
        sims = []
        for d in curated[draw12.norm(gene)]:
            if pheno in idx and d in idx:
                sims.append((float(vec[idx[pheno]] @ vec[idx[d]]), d))
        return max(sims) if sims else (0.0, "")

    kept = []
    for gene, pheno in cheap:
        if len(kept) == N:
            break
        cos, near = nearest(gene, pheno)
        if cos >= THRESHOLD:
            refused["cosine"] += 1
            continue
        kept.append({
            "gene": gene,
            "phenotype": pheno,
            "gencc_class": "Absent",
            "max_cosine": f"{cos:.3f}",
            "nearest_curated_disease": near,
        })
    if len(kept) < N:
        raise SystemExit(
            f"only {len(kept)} Absent pairs (refused {dict(refused)})")

    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(kept[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(kept[:N])
    print(f"wrote {N} pairs to {os.path.relpath(OUT, ROOT)}")
    print(f"pool {len(genes)} ordinary genes × {len(phenos)} traits")
    print("refused", dict(refused))
    print("nearest cosine",
          min(float(r["max_cosine"]) for r in kept),
          "–", max(float(r["max_cosine"]) for r in kept))


if __name__ == "__main__":
    main()
