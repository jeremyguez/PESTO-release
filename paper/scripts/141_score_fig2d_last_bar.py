"""Fig. 2d, 2023–2025 bar only: does an arm beat the published PESTO and the
closed-book model on average precision?

Rebuilds the bar exactly as scripts/21_figure2.R (year_bin_prep) does: the 40
Limited pairs first reported in 2023–2025 (results/prior_year10) against the
100 Absent pairs of Fig. 2a, score = (4 Est + 3 Exi + 2 Hyp + 1 Nov) / 100,
average precision with ties grouped, 2000 bootstrap draws for the interval
and a one-sided paired permutation (20000 swaps). The published numbers are
reproduced first as a check, then each arm table given on the command line is
scored on the same 140 pairs and compared, paired, with both.

  python3 scripts/141_score_fig2d_last_bar.py results/arm_currentv2-cheap_fig2d_ap_2023_2025.tsv
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
BOOT_B, PERM_B = 2000, 20000


def score(est, exi, hyp, nov):
    return (4 * float(est or 0) + 3 * float(exi or 0)
            + 2 * float(hyp or 0) + float(nov or 0)) / 100


def read_tsv(path):
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def key(gene, phenotype):
    return gene.upper(), phenotype


def limited():
    know = {}
    for r in read_tsv(RESULTS / "prior_year10" / "knowledge.tsv"):
        if r["year"] and r["lit_score"] and 2023 <= int(r["year"]) <= 2025:
            know[key(r["gene"], r["phenotype"])] = float(r["lit_score"])
    pesto = {}
    for r in json.loads((RESULTS / "prior_year10" / "pesto.json").read_text()):
        k = key(r["gene"], r["phenotype"])
        if k in know:
            d = r["distribution"]
            pesto[k] = score(d.get("Established"), d.get("Existing"),
                             d.get("Hypothesized"), d.get("Novel"))
    assert len(know) == len(pesto) == 40, (len(know), len(pesto))
    return pesto, know


def absent():
    pesto, know = {}, {}
    p30 = {key(r["gene"], r["phenotype"]): score(r["p_established"], r["p_existing"],
                                                 r["p_hypothesized"], r["p_novel"])
           for r in read_tsv(ROOT / "benchmark" / "gencc_clingen_g2p" / "pesto.tsv")
           if r["gencc_class"] == "Absent"}
    k30 = {key(r["gene"], r["phenotype"]): float(r["lit_score"])
           for r in read_tsv(RESULTS / "ablation_fig2" / "knowledge.tsv")
           if int(r["draw"]) == 1 and r["lit_score"] not in ("", "NA")}
    blocks = [(p30, k30)]
    for d in ("fig2a_absent_extra20", "fig2a_absent_extra20b", "fig2a_absent_extra30"):
        p = {key(r["gene"], r["phenotype"]): float(r["lit_score"])
             for r in read_tsv(RESULTS / d / "pesto.tsv")
             if not r.get("error") and r["lit_score"] not in ("", "NA")}
        k = {key(r["gene"], r["phenotype"]): float(r["lit_score"])
             for r in read_tsv(RESULTS / d / "knowledge.tsv")
             if r["lit_score"] not in ("", "NA")}
        blocks.append((p, k))
    for p, k in blocks:
        for kk in p.keys() & k.keys():
            pesto[kk], know[kk] = p[kk], k[kk]
    assert len(pesto) == len(know) == 100, (len(pesto), len(know))
    return pesto, know


def ap(s, y):
    o = np.argsort(-s, kind="stable")
    s, y = s[o], y[o]
    ends = np.flatnonzero(np.r_[s[1:] != s[:-1], True])
    tp = np.cumsum(y)[ends]
    recall = tp / y.sum()
    return float(np.sum(np.diff(np.r_[0, recall]) * tp / (ends + 1)))


def boot_ci(s, y, rng):
    n = len(y)
    vals = [ap(s[i], y[i]) for i in (rng.integers(0, n, n) for _ in range(BOOT_B))]
    return np.quantile(vals, [0.025, 0.975])


def paired_perm(a, b, y, rng):
    obs = ap(a, y) - ap(b, y)
    ge = 0
    for _ in range(PERM_B):
        flip = rng.random(len(y)) < 0.5
        if ap(np.where(flip, b, a), y) - ap(np.where(flip, a, b), y) >= obs:
            ge += 1
    return (1 + ge) / (1 + PERM_B)


def arm_scores(path):
    out = {}
    for r in read_tsv(path):
        if r.get("error"):
            continue
        out[key(r["gene"], r["phenotype"])] = score(
            r["p_established"], r["p_existing"], r["p_hypothesized"], r["p_novel"])
    return out


def main():
    lim_p, lim_k = limited()
    abs_p, abs_k = absent()
    pairs = sorted(lim_p) + sorted(abs_p)
    y = np.array([1] * 40 + [0] * 100)
    methods = {"PESTO (published)": {**lim_p, **abs_p},
               "Model knowledge": {**lim_k, **abs_k}}
    for path in sys.argv[1:]:
        s = arm_scores(path)
        missing = [p for p in pairs if p not in s]
        if missing:
            sys.exit(f"{path}: {len(missing)} of 140 pairs missing, e.g. {missing[:3]}")
        methods[Path(path).stem] = s

    vec = {m: np.array([s[p] for p in pairs]) for m, s in methods.items()}
    print(f"2023–2025 bar: 40 Limited vs 100 Absent, chance = {40 / 140:.3f}\n")
    for m, v in vec.items():
        lo, hi = boot_ci(v, y, np.random.default_rng(1))
        print(f"  {m:45s} AP {ap(v, y):.3f}  (95% CI {lo:.3f}–{hi:.3f})")
    print()
    names = list(vec)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            hi, lo = (a, b) if ap(vec[a], y) >= ap(vec[b], y) else (b, a)
            p = paired_perm(vec[hi], vec[lo], y, np.random.default_rng(1))
            print(f"  {hi} > {lo}: one-sided paired permutation P = {p:.4f}")


if __name__ == "__main__":
    main()
