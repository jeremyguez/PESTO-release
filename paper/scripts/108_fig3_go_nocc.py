"""GO signatures on the 509-gene table, without cellular-component terms.

Nucleus and other localizations ate the ASD process terms under the overlap
filter. This run tests only Biological Process and Molecular Function, so
Benjamini-Hochberg is over those tests alone. Same 509-gene scores, no DD,
same pick rules as scripts/106_fig3_random50b.py. Does not overwrite
results/fig3_random50b/ or figures/figure3_ad_r50.*

  python3 scripts/108_fig3_go_nocc.py
  Rscript scripts/109_figure3_ad_r50_nocc.R
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
SPEC106 = importlib.util.spec_from_file_location(
    "s106", ROOT / "scripts" / "106_fig3_random50b.py")
s106 = importlib.util.module_from_spec(SPEC106)
SPEC106.loader.exec_module(s106)

GMT = (
    ROOT / "data" / "go" / "GO_Biological_Process_2026.gmt",
    ROOT / "data" / "go" / "GO_Molecular_Function_2026.gmt",
)
OUT = ROOT / "results" / "fig3_random50b_nocc"


def read_gmt(universe):
    terms = []
    for path in GMT:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                name = f[0]
                if name.endswith(" CC"):
                    continue
                genes = {g.upper() for g in f[2:] if g} & universe
                if len(genes) >= 5:
                    terms.append((name, genes))
    return terms


def enrich(score):
    universe = set(score)
    fg = {g for g in universe if score[g] >= 3}
    bg = universe - fg
    rows = []
    for name, genes in read_gmt(universe):
        hit = genes & fg
        if not hit:
            continue
        k, K, N = len(hit), len(genes), len(universe)
        _, p = fisher_exact(
            [[k, len(fg) - k], [len(bg & genes), len(bg) - len(bg & genes)]],
            alternative="greater")
        rows.append({
            "item": name, "k": k, "K": K, "n_fg": len(fg), "N": N,
            "p": p, "lor": s106.g.log2_or(k, K, len(fg), N), "hit": hit,
            "members": frozenset(genes),
        })
    for r, q in zip(rows, s106.g.m.bh([r["p"] for r in rows])):
        r["q"] = q
    return {r["item"]: r for r in rows}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s106.OUT = OUT
    s106.SIG = OUT / "go_signatures.tsv"
    s106.CORR = OUT / "pheno_go_correlation.tsv"
    s106.HEAT = OUT / "pheno_heatmap.tsv"
    s106.g.enrich = enrich
    s106.go()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
