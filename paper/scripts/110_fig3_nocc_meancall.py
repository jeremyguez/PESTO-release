"""Figure-3 nocc scores from mean_call, reaching if mean > 2.5.

Re-reads the hundred-point distributions in the flow cache (no new model
calls). mean_call is pesto.cli.mean_call (round of the 1-4 average).
GO uses Biological Process + Molecular Function only, foreground = mean > 2.5.
Does not overwrite the argmax nocc figure.

  python3 scripts/110_fig3_nocc_meancall.py
  Rscript scripts/111_figure3_ad_r50_nocc_mean.R

  python3 scripts/110_fig3_nocc_meancall.py --pick-or
  Rscript scripts/111_figure3_ad_r50_nocc_mean.R --pick-or
"""
from __future__ import annotations

import importlib.util
import json
import sys
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

PAIRS = ROOT / "results" / "fig3_universe509_ad_scores.tsv"
UNIVERSE = ROOT / "results" / "fig3_universe509_meancall_scores.tsv"
OUT = ROOT / "results" / "fig3_random50b_nocc_mean"
PICK_BY = "or" if "--pick-or" in sys.argv else "q"
if PICK_BY == "or":
    OUT = ROOT / "results" / "fig3_random50b_nocc_mean_or"
RUN_DIRS = [
    ROOT / "benchmark" / "asc_ref299" / "runs" / "flow",
    ROOT / "benchmark" / "asd_extra101" / "runs" / "flow",
    ROOT / "benchmark" / "scz_extra101" / "runs" / "flow",
    ROOT / "benchmark" / "extra50_random" / "runs" / "flow",
    ROOT / "benchmark" / "asd" / "runs" / "flow",
    ROOT / "results" / "fig3_four_phenos" / "runs" / "flow",
    ROOT / "results" / "fig3_random50" / "runs" / "flow",
    ROOT / "results" / "fig3_alzheimer" / "runs" / "flow",
    ROOT / "results" / "fig3_random50b" / "runs" / "flow",
]
WEIGHT = {"Established": 4, "Existing": 3, "Hypothesized": 2, "Novel": 1}
ORDER = ("Novel", "Hypothesized", "Existing", "Established")
FP = "__ff136250f372.json"


def mean_of(dist):
    total = sum(dist.get(k, 0) for k in WEIGHT) or 1
    return sum(WEIGHT[k] * dist.get(k, 0) for k in WEIGHT) / total


def mean_call(rank):
    return ORDER[max(1, min(4, round(rank))) - 1]


def load_dist():
    found = {}
    for folder in RUN_DIRS:
        if not folder.exists():
            continue
        for path in folder.glob(f"*{FP}"):
            data = json.loads(path.read_text())
            dist = (data.get("verdict") or {}).get("distribution")
            if not dist:
                continue
            found[(str(data["gene"]).upper(), data["phenotype"])] = {
                "dist": dist,
                "argmax": data["verdict"].get("call") or "",
            }
    return found


def assemble():
    raw = pd.read_csv(PAIRS, sep="\t")
    raw["gene"] = raw.gene.str.upper()
    cache = load_dist()
    rows = []
    miss = 0
    for r in raw.itertuples():
        rec = cache.get((r.gene, r.phenotype))
        if rec is None:
            miss += 1
            continue
        rank = mean_of(rec["dist"])
        call = mean_call(rank)
        rows.append({
            "gene": r.gene,
            "phenotype": r.phenotype,
            "lit_score": rank,
            "mean_call": call,
            "argmax": rec["argmax"],
            "lit_mean": call,
            "reaching": bool(rank > 2.5),
            "source": getattr(r, "source", ""),
            "set": getattr(r, "set", ""),
        })
    if miss:
        raise SystemExit(f"missing distributions for {miss} pairs")
    out = pd.DataFrame(rows).sort_values(["gene", "phenotype"])
    out.to_csv(UNIVERSE, sep="\t", index=False)
    print(f"wrote {UNIVERSE.relative_to(ROOT)}  {len(out)} rows  "
          f"{out.gene.nunique()} genes")
    disagree = (out["mean_call"] != out["argmax"]).sum()
    print(f"  mean_call != argmax: {disagree}/{len(out)}")
    print(f"  mean>2.5 != mean_call in {{Existing,Established}}: "
          f"{((out['lit_score'] > 2.5) != out['mean_call'].isin(('Existing', 'Established'))).sum()}")
    phenos = [
        "Autism Spectrum Disorder", "epilepsy", "schizophrenia",
        "bipolar disorder", "type 2 diabetes", "Alzheimer's disease",
    ]
    for p in phenos:
        sub = out[out.phenotype == p]
        n_arg = sub["argmax"].isin(("Existing", "Established")).sum()
        n_mean = int(sub["reaching"].sum())
        print(f"  {p:<28} argmax {n_arg:3}  mean>2.5 {n_mean:3}")


def go():
    raw = pd.read_csv(UNIVERSE, sep="\t")
    raw["gene"] = raw.gene.str.upper()
    n_genes = raw.gene.nunique()

    def score_of(pheno):
        sub = raw[raw.phenotype == pheno]
        return {gene: (4 if mean > 2.5 else 1)
                for gene, mean in zip(sub.gene, sub["lit_score"])}

    s106.g.enrich = s108.enrich
    s106.OUT = OUT
    s106.SIG = OUT / "go_signatures.tsv"
    s106.CORR = OUT / "pheno_go_correlation.tsv"
    s106.HEAT = OUT / "pheno_heatmap.tsv"

    def go_mean():
        sig_tables, chosen = {}, {}
        for pheno, short in s106.SIG_PHENOS:
            enrich = s108.enrich(score_of(pheno))
            sig_tables[short] = enrich
            chosen[short] = s106.pick(enrich, by=PICK_BY)
            n_fg = next(iter(enrich.values()))["n_fg"]
            n_sig = sum(r["q"] < s106.FDR and r["k"] >= s106.MIN_K_PICK
                        for r in enrich.values())
            head = chosen[short][0]["item"] if chosen[short] else "-- none --"
            print(f"{short:8} reaching {n_fg:3}/{n_genes}  "
                  f"FDR<{s106.FDR}: {n_sig:3}  picked {len(chosen[short])}  "
                  f"top: {head[:44]}")
            cation = [r for r in enrich.values() if "0098655" in r["item"]]
            if cation and short == "epilepsy":
                r = cation[0]
                print(f"         cation k={r['k']} K={r['K']} n_fg={r['n_fg']} "
                      f"p={r['p']:.4e} q={r['q']:.4f} lor={r['lor']:.2f}")

        selected = {r["item"] for picks in chosen.values() for r in picks}
        order = {short: i for i, (_, short) in enumerate(s106.SIG_PHENOS)}
        block = {}
        for item in selected:
            owners = [(sig_tables[s].get(item), s) for _, s in s106.SIG_PHENOS]
            owners = [(r, s) for r, s in owners
                      if r is not None and r["q"] < s106.FDR]
            best = max(owners, key=lambda rs: rs[0]["lor"])
            block[item] = (best[1], best[0]["q"], best[0]["lor"])

        rows = []
        for item in sorted(selected,
                           key=lambda i: (order[block[i][0]],
                                          -block[i][2] if PICK_BY == "or"
                                          else block[i][1])):
            source = block[item][0]
            for _, short in s106.SIG_PHENOS:
                r = sig_tables[short].get(item)
                if r is None:
                    continue
                rows.append({
                    "item": item, "source": source, "phenotype": short,
                    "k": r["k"], "K": r["K"], "n_fg": r["n_fg"],
                    "lor": r["lor"], "q": r["q"], "p": r["p"],
                    "frac": r["k"] / r["K"],
                })
        OUT.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(s106.SIG, sep="\t", index=False)
        print(f"wrote {s106.SIG.relative_to(ROOT)}  terms={len(selected)}")

        tables = {}
        for pheno, short in s106.CLUST_PHENOS:
            enrich = s108.enrich(score_of(pheno))
            tables[short] = enrich
            n_fg = next(iter(enrich.values()))["n_fg"]
            n_sig = sum(r["q"] < s106.FDR for r in enrich.values())
            print(f"{short:8} cluster reaching {n_fg:3}/{n_genes}  "
                  f"FDR<{s106.FDR}: {n_sig:3}")

        # Reuse the rest of 106.go by writing SIG first then calling from
        # clustering onward: simpler to keep the linkage block here.
        import numpy as np
        from scipy.cluster.hierarchy import cophenet, dendrogram, linkage
        from scipy.spatial.distance import squareform
        from scipy.stats import spearmanr

        shared = set.intersection(*(set(t) for t in tables.values()))
        shared = {i for i in shared
                  if next(iter(tables.values()))[i]["K"] >= s106.MIN_K_CLUST}
        items = sorted(shared)
        shorts = [s for _, s in s106.CLUST_PHENOS]
        mat = np.array([[tables[s][i]["lor"] for i in items] for s in shorts])
        rho = spearmanr(mat, axis=1).statistic
        corr = pd.DataFrame(rho, index=shorts, columns=shorts)
        corr.to_csv(s106.CORR, sep="\t")
        print(f"\nclustering on {len(items)} terms with K >= {s106.MIN_K_CLUST}")
        print(corr.round(3).to_string())

        dist = 1 - corr.values
        np.fill_diagonal(dist, 0.0)
        link = linkage(squareform(dist, checks=False), method="average")
        coph = cophenet(link, squareform(dist, checks=False))[0]
        leaves = dendrogram(link, no_plot=True, labels=shorts)["ivl"]
        print(f"\naverage-linkage order: {' '.join(leaves)}")
        print(f"cophenetic correlation: {coph:.3f}")

        shown = list(dict.fromkeys(pd.read_csv(s106.SIG, sep="\t").item))
        heat = []
        for item in shown:
            for short in shorts:
                r = tables[short].get(item)
                if r is None:
                    continue
                heat.append({
                    "item": item, "phenotype": short, "lor": r["lor"],
                    "q": r["q"], "k": r["k"], "K": r["K"], "n_fg": r["n_fg"],
                })
        pd.DataFrame(heat).to_csv(s106.HEAT, sep="\t", index=False)
        print(f"wrote {s106.CORR.relative_to(ROOT)}\n"
              f"wrote {s106.HEAT.relative_to(ROOT)}")

    go_mean()


def main():
    if PICK_BY == "q":
        assemble()
    elif not UNIVERSE.exists():
        raise SystemExit(f"missing {UNIVERSE}; run without --pick-or first")
    print(f"pick by {PICK_BY}  → {OUT.relative_to(ROOT)}")
    go()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
