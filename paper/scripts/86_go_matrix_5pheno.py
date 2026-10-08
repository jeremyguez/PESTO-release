"""GO enrichment across five phenotypes on the 400-gene universe.

Drops developmental disorder: 295 of 400 genes reach it, so its background
is 105 genes and any contrast against it is thin by construction.

Effect size is the log2 odds ratio with a Haldane-Anscombe 0.5 correction,
not the fold enrichment. Fold is capped at N / n_reaching, which runs from
1.68 (ASD) to 7.84 (bipolar), so fold cannot be compared across columns.
The odds ratio is what Fisher's exact test tests and has no such ceiling.

Terms are picked per phenotype: FDR < 0.05, sorted by log2 OR, kept greedily
only when the reaching genes overlap an already-kept term by less than
OVERLAP. The union of those picks becomes the matrix, scored in every column.

Writes results/go_matrix_5pheno.tsv for scripts/87_figure_go_matrix.R.

Usage: python3 scripts/86_go_matrix_5pheno.py
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "fig3", ROOT / "scripts" / "42_fig3_cd_bars.py")
m = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)

SCORES = ROOT / "results" / "fig3_universe400_scores.tsv"
OUT = ROOT / "results" / "go_matrix_5pheno.tsv"
W = {"Novel": 1, "Hypothesized": 2, "Existing": 3, "Established": 4}

PHENOS = [
    ("Autism Spectrum Disorder", "ASD"),
    ("epilepsy", "epilepsy"),
    ("schizophrenia", "SCZ"),
    ("bipolar disorder", "bipolar"),
    ("type 2 diabetes", "T2D"),
]
FDR = 0.05
OVERLAP = 0.5
TOP_PER_PHENO = 5
MIN_K = 5


def log2_or(k, K, n_fg, N):
    """Haldane-Anscombe corrected, so a fully reaching term stays finite."""
    a, b = k + 0.5, K - k + 0.5
    c, d = n_fg - k + 0.5, N - K - n_fg + k + 0.5
    return float(np.log2((a / b) / (c / d)))


def enrich(score):
    universe = set(score)
    fg = {g for g in universe if score[g] >= 3}
    bg = universe - fg
    rows = []
    for name, genes in m.read_gmt(universe):
        hit = genes & fg
        if not hit:
            continue
        k, K, N = len(hit), len(genes), len(universe)
        _, p = fisher_exact(
            [[k, len(fg) - k], [len(bg & genes), len(bg) - len(bg & genes)]],
            alternative="greater")
        rows.append({
            "item": name, "k": k, "K": K, "n_fg": len(fg), "N": N,
            "p": p, "lor": log2_or(k, K, len(fg), N), "hit": hit,
        })
    for r, q in zip(rows, m.bh([r["p"] for r in rows])):
        r["q"] = q
    return {r["item"]: r for r in rows}


def pick(by_item):
    """Strongest by odds ratio, skipping terms that repeat earlier gene sets."""
    cand = [r for r in by_item.values()
            if r["q"] < FDR and r["k"] >= MIN_K]
    cand.sort(key=lambda r: -r["lor"])
    kept = []
    for r in cand:
        if any(len(r["hit"] & s["hit"]) / min(len(r["hit"]), len(s["hit"]))
               >= OVERLAP for s in kept):
            continue
        kept.append(r)
        if len(kept) == TOP_PER_PHENO:
            break
    return kept


def main():
    raw = pd.read_csv(SCORES, sep="\t")
    raw = raw[raw.lit_mean.isin(W)]

    tables, chosen = {}, {}
    for pheno, short in PHENOS:
        sub = raw[raw.phenotype == pheno]
        score = {g.upper(): W[v] for g, v in zip(sub.gene, sub.lit_mean)}
        tables[short] = enrich(score)
        chosen[short] = pick(tables[short])
        n_fg = next(iter(tables[short].values()))["n_fg"]
        print(f"{short:8} reaching {n_fg:3}/400  "
              f"FDR<{FDR}: {sum(r['q'] < FDR for r in tables[short].values()):3}"
              f"  picked {len(chosen[short])}")

    order, seen = [], set()
    for _, short in PHENOS:
        for r in chosen[short]:
            if r["item"] not in seen:
                seen.add(r["item"])
                order.append((r["item"], short))

    rows = []
    for item, source in order:
        for _, short in PHENOS:
            r = tables[short].get(item)
            if r is None:
                continue
            rows.append({
                "item": item, "source": source, "phenotype": short,
                "k": r["k"], "K": r["K"], "n_fg": r["n_fg"],
                "lor": r["lor"], "q": r["q"], "p": r["p"],
                "frac": r["k"] / r["K"],
            })
    out = pd.DataFrame(rows)
    out.to_csv(OUT, sep="\t", index=False)
    print(f"\nwrote {OUT}  terms={len(order)}  rows={len(out)}")

    wide = out.pivot(index="item", columns="phenotype", values="lor")
    wide = wide[[s for _, s in PHENOS]].reindex([i for i, _ in order])
    print("\nlog2 odds ratio\n")
    print(wide.round(2).to_string())


if __name__ == "__main__":
    main()
