"""Top-5 GO terms per phenotype, overlap filter preferring GO parents.

Candidates are FDR < 0.05 with k >= 5, walked in descending odds ratio.
If two terms overlap (reaching-gene Jaccard-style min-overlap >= 0.5):
the GMT gene-set superset is treated as the parent and kept, replacing a
child already stored; with no inclusion, the first (higher OR) stays.
Writes go_signatures.tsv for the figure-3 blocks panel (five terms per
phenotype, source = nominator). Does not overwrite the nocc_mean tables.

  python3 scripts/113_fig3_top5_parent.py
  Rscript scripts/114_figure3_top5_parent.R
"""
from __future__ import annotations

import importlib.util
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPEC106 = importlib.util.spec_from_file_location(
    "s106", ROOT / "scripts" / "106_fig3_random50b.py")
s106 = importlib.util.module_from_spec(SPEC106)
SPEC106.loader.exec_module(s106)
SPEC108 = importlib.util.spec_from_file_location(
    "s108", ROOT / "scripts" / "108_fig3_go_nocc.py")
s108 = importlib.util.module_from_spec(SPEC108)
SPEC108.loader.exec_module(s108)
SPEC110 = importlib.util.spec_from_file_location(
    "s110", ROOT / "scripts" / "110_fig3_nocc_meancall.py")
s110 = importlib.util.module_from_spec(SPEC110)
SPEC110.loader.exec_module(s110)

UNIVERSE = s110.UNIVERSE
OUT = ROOT / "results" / "fig3_random50b_nocc_mean_top5"
N_KEEP = 5


def overlaps(a, b):
    n = min(len(a["hit"]), len(b["hit"]))
    if not n:
        return False
    return len(a["hit"] & b["hit"]) / n >= s106.OVERLAP


def pick_parent(by_item, n=N_KEEP):
    cand = [r for r in by_item.values()
            if r["q"] < s106.FDR and r["k"] >= s106.MIN_K_PICK]
    cand.sort(key=lambda r: (-r["lor"], r["q"]))
    kept = []
    swaps = []
    for r in cand:
        child_idx = []
        reject = False
        for i, s in enumerate(kept):
            if not overlaps(r, s):
                continue
            if r["members"] > s["members"]:
                child_idx.append(i)
            else:
                reject = True
                break
        if reject:
            continue
        if child_idx:
            for i in sorted(child_idx, reverse=True):
                child = kept.pop(i)
                swaps.append((child["item"], r["item"]))
            kept.append(r)
        elif len(kept) < n:
            kept.append(r)
    kept.sort(key=lambda r: (-r["lor"], r["q"]))
    return kept, swaps


def score_of(raw, pheno):
    sub = raw[raw.phenotype == pheno]
    return {gene: (4 if mean > 2.5 else 1)
            for gene, mean in zip(sub.gene, sub["lit_score"])}


def main():
    raw = pd.read_csv(UNIVERSE, sep="\t")
    raw["gene"] = raw.gene.str.upper()
    n_genes = raw.gene.nunique()
    tables, chosen, swaps_of = {}, {}, {}
    for pheno, short in s106.SIG_PHENOS:
        enrich = s108.enrich(score_of(raw, pheno))
        tables[short] = enrich
        chosen[short], swaps_of[short] = pick_parent(enrich)
        n_fg = next(iter(enrich.values()))["n_fg"]
        n_sig = sum(r["q"] < s106.FDR and r["k"] >= s106.MIN_K_PICK
                    for r in enrich.values())
        print(f"{short:8} reaching {n_fg:3}/{n_genes}  "
              f"FDR<{s106.FDR}: {n_sig:3}  kept {len(chosen[short])}")
        for r in chosen[short]:
            print(f"         OR={2 ** r['lor']:.1f}  q={r['q']:.3g}  "
                  f"K={r['K']:3}  {r['item'].split(' (GO:')[0]}")
        for child, parent in swaps_of[short]:
            print(f"         replaced {child.split(' (GO:')[0]}  →  "
                  f"{parent.split(' (GO:')[0]}")

    shorts = [s for _, s in s106.SIG_PHENOS]
    rows = []
    sig = []
    for _, short in s106.SIG_PHENOS:
        for rank, r in enumerate(chosen[short], start=1):
            yid = f"{short}\t{r['item']}"
            for col in shorts:
                cell = tables[col].get(r["item"])
                if cell is None:
                    continue
                rec = {
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
                }
                rows.append(rec)
                sig.append(rec)
    OUT.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "go_top5.tsv", sep="\t", index=False)
    pd.DataFrame(sig).to_csv(OUT / "go_signatures.tsv", sep="\t", index=False)
    corr_src = ROOT / "results" / "fig3_random50b_nocc_mean" / "pheno_go_correlation.tsv"
    shutil.copy(corr_src, OUT / "pheno_go_correlation.tsv")
    n_terms = out.drop_duplicates(["nominator", "item"]).shape[0]
    print(f"wrote {OUT.relative_to(ROOT)}/go_signatures.tsv  {n_terms} terms")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
