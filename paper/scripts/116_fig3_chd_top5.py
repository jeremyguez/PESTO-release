"""Figure-3 top-5 parent GO on the 550-gene table, with CHD.

Same rules as scripts/113_fig3_top5_parent.py (mean > 2.5, no DD, no CC,
parent-preferring overlap, five terms per phenotype). Adds congenital
heart disease to the signature and clustering columns. Drops ZNF710, the
51st extra101 random gene (already out of Fig. 2b), so random draws are
50+50+50+50. Does not overwrite the 551-gene meancall table.

  python3 scripts/116_fig3_chd_top5.py
  Rscript scripts/117_figure3_chd_top5.R
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet, dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
SPEC108 = importlib.util.spec_from_file_location(
    "s108", ROOT / "scripts" / "108_fig3_go_nocc.py")
s108 = importlib.util.module_from_spec(SPEC108)
SPEC108.loader.exec_module(s108)
SPEC110 = importlib.util.spec_from_file_location(
    "s110", ROOT / "scripts" / "110_fig3_nocc_meancall.py")
s110 = importlib.util.module_from_spec(SPEC110)
SPEC110.loader.exec_module(s110)
SPEC113 = importlib.util.spec_from_file_location(
    "s113", ROOT / "scripts" / "113_fig3_top5_parent.py")
s113 = importlib.util.module_from_spec(SPEC113)
SPEC113.loader.exec_module(s113)

PAIRS = ROOT / "results" / "fig3_universe551_chd_scores.tsv"
UNIVERSE = ROOT / "results" / "fig3_universe550_meancall_scores.tsv"
OUT = ROOT / "results" / "fig3_chd_nocc_mean_top5"
CHD_RUNS = ROOT / "results" / "fig3_chd" / "runs" / "flow"
DROP = {"ZNF710"}

SIG_PHENOS = [
    ("Autism Spectrum Disorder", "ASD"),
    ("epilepsy", "epilepsy"),
    ("congenital heart disease", "CHD"),
    ("bipolar disorder", "bipolar"),
    ("schizophrenia", "SCZ"),
    ("type 2 diabetes", "T2D"),
    ("Alzheimer's disease", "AD"),
]
CLUST_PHENOS = list(SIG_PHENOS)


def assemble():
    s110.PAIRS = PAIRS
    s110.UNIVERSE = UNIVERSE
    if CHD_RUNS not in s110.RUN_DIRS:
        s110.RUN_DIRS.append(CHD_RUNS)
    s110.assemble()
    raw = pd.read_csv(UNIVERSE, sep="\t")
    raw["gene"] = raw.gene.str.upper()
    kept = raw[~raw.gene.isin(DROP)]
    kept.to_csv(UNIVERSE, sep="\t", index=False)
    print(f"dropped {sorted(DROP)}  {kept.gene.nunique()} genes  "
          f"{len(kept)} rows")


def score_of(raw, pheno):
    sub = raw[raw.phenotype == pheno]
    return {gene: (4 if mean > 2.5 else 1)
            for gene, mean in zip(sub.gene, sub["lit_score"])}


def go():
    raw = pd.read_csv(UNIVERSE, sep="\t")
    raw["gene"] = raw.gene.str.upper()
    n_genes = raw.gene.nunique()
    tables, chosen, swaps_of = {}, {}, {}
    for pheno, short in SIG_PHENOS:
        enrich = s108.enrich(score_of(raw, pheno))
        tables[short] = enrich
        chosen[short], swaps_of[short] = s113.pick_parent(enrich)
        n_fg = next(iter(enrich.values()))["n_fg"]
        n_sig = sum(r["q"] < s108.s106.FDR and r["k"] >= s108.s106.MIN_K_PICK
                    for r in enrich.values())
        print(f"{short:8} reaching {n_fg:3}/{n_genes}  "
              f"FDR<{s108.s106.FDR}: {n_sig:3}  kept {len(chosen[short])}")
        for r in chosen[short]:
            print(f"         OR={2 ** r['lor']:.1f}  q={r['q']:.3g}  "
                  f"K={r['K']:3}  {r['item'].split(' (GO:')[0]}")
        for child, parent in swaps_of[short]:
            print(f"         replaced {child.split(' (GO:')[0]}  →  "
                  f"{parent.split(' (GO:')[0]}")

    shorts = [s for _, s in SIG_PHENOS]
    rows = []
    for _, short in SIG_PHENOS:
        for rank, r in enumerate(chosen[short], start=1):
            yid = f"{short}\t{r['item']}"
            for col in shorts:
                cell = tables[col].get(r["item"])
                if cell is None:
                    continue
                rows.append({
                    "nominator": short,
                    "rank": rank,
                    "item": r["item"],
                    "yid": yid,
                    "source": short,
                    "phenotype": col,
                    "lor": cell["lor"],
                    "q": cell["q"],
                    "k": cell["k"],
                    "K": cell["K"],
                    "n_fg": cell["n_fg"],
                    "p": cell["p"],
                    "frac": cell["k"] / cell["K"],
                })
    OUT.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "go_top5.tsv", sep="\t", index=False)
    out.to_csv(OUT / "go_signatures.tsv", sep="\t", index=False)

    clust = {}
    for pheno, short in CLUST_PHENOS:
        enrich = s108.enrich(score_of(raw, pheno))
        clust[short] = enrich
        n_fg = next(iter(enrich.values()))["n_fg"]
        n_sig = sum(r["q"] < s108.s106.FDR for r in enrich.values())
        print(f"{short:8} cluster reaching {n_fg:3}/{n_genes}  "
              f"FDR<{s108.s106.FDR}: {n_sig:3}")

    shared = set.intersection(*(set(t) for t in clust.values()))
    shared = {i for i in shared
              if next(iter(clust.values()))[i]["K"] >= s108.s106.MIN_K_CLUST}
    items = sorted(shared)
    cshorts = [s for _, s in CLUST_PHENOS]
    mat = np.array([[clust[s][i]["lor"] for i in items] for s in cshorts])
    rho = spearmanr(mat, axis=1).statistic
    corr = pd.DataFrame(rho, index=cshorts, columns=cshorts)
    corr.to_csv(OUT / "pheno_go_correlation.tsv", sep="\t")
    print(f"\nclustering on {len(items)} terms with K >= {s108.s106.MIN_K_CLUST}")
    print(corr.round(3).to_string())

    dist = 1 - corr.values
    np.fill_diagonal(dist, 0.0)
    link = linkage(squareform(dist, checks=False), method="average")
    coph = cophenet(link, squareform(dist, checks=False))[0]
    leaves = dendrogram(link, no_plot=True, labels=cshorts)["ivl"]
    print(f"\naverage-linkage order: {' '.join(leaves)}")
    print(f"cophenetic correlation: {coph:.3f}")

    n_terms = out.drop_duplicates(["nominator", "item"]).shape[0]
    print(f"wrote {OUT.relative_to(ROOT)}/go_signatures.tsv  {n_terms} terms")


def main():
    assemble()
    go()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
