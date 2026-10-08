#!/usr/bin/env python3
"""Data for the draft Extended Data Fig. 2: PESTO-titles beside PESTO and the
model's internal knowledge on the 2023–2025 bar of Fig. 2d.

Answers Kyle's comment that PESTO-titles is quoted in the text but drawn
nowhere. Same 140 pairs as Fig. 2d (40 Limited first reported 2023–2025
against the 100 Absent of Fig. 2a), same score, AP, bootstrap and one-sided
paired permutation as scripts/141_score_fig2d_last_bar.py, whose loaders are
reused. No model call.

Writes
  results/ed2_titles_ap.tsv        AP and 95% CI per method
  results/ed2_titles_tests.tsv     paired permutation P per comparison
  results/ed2_titles_verdicts.tsv  verdict counts on the 40 Limited pairs
  results/ed2_titles_cost.tsv      USD per pair on the Fig. 2a,b pairs

  python3 scripts/143_ed2_titles_data.py
"""
from __future__ import annotations

import csv
import importlib.util
import json
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

_spec = importlib.util.spec_from_file_location(
    "fig2d", ROOT / "scripts" / "141_score_fig2d_last_bar.py")
fig2d = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fig2d)

VERDICTS = ["Novel", "Hypothesized", "Existing", "Established"]
METHODS = ["Model knowledge", "PESTO-titles", "PESTO"]


def verdict_of(score):
    return VERDICTS[max(1, min(4, round(score))) - 1]


def titles_scores():
    lim = {}
    for r in json.loads((RESULTS / "prior_year10" / "titles_2023_2025.json").read_text()):
        d = r.get("distribution") or {}
        if d:
            lim[fig2d.key(r["gene"], r["phenotype"])] = fig2d.score(
                d.get("Established"), d.get("Existing"), d.get("Hypothesized"), d.get("Novel"))
    absent = {}
    for path, keep in ((RESULTS / "ablation_fig2" / "titles.tsv", "gencc_absent"),
                       (RESULTS / "fig2a_absent_extra70" / "titles.tsv", None)):
        for r in fig2d.read_tsv(path):
            if keep and r.get("set") != keep:
                continue
            if r.get("error") or r.get("lit_score") in (None, "", "NA"):
                continue
            absent[fig2d.key(r["gene"], r["phenotype"])] = float(r["lit_score"])
    return lim, absent


def write(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("wrote", path.relative_to(ROOT))


def main():
    lim_p, lim_k = fig2d.limited()
    abs_p, abs_k = fig2d.absent()
    lim_t, abs_t = titles_scores()
    limited = sorted(lim_p)
    pairs = limited + sorted(abs_p)
    y = np.array([1] * len(limited) + [0] * len(abs_p))
    scores = {"Model knowledge": {**lim_k, **abs_k},
              "PESTO-titles": {**lim_t, **abs_t},
              "PESTO": {**lim_p, **abs_p}}
    for m, s in scores.items():
        missing = [p for p in pairs if p not in s]
        assert not missing, (m, missing[:3])
    vec = {m: np.array([s[p] for p in pairs]) for m, s in scores.items()}

    ap_rows = []
    for m in METHODS:
        lo, hi = fig2d.boot_ci(vec[m], y, np.random.default_rng(1))
        ap_rows.append({"method": m, "ap": round(fig2d.ap(vec[m], y), 4),
                        "lo": round(float(lo), 4), "hi": round(float(hi), 4)})
        print(f"{m:16s} AP {ap_rows[-1]['ap']:.3f} ({lo:.3f}–{hi:.3f})")
    write(RESULTS / "ed2_titles_ap.tsv", ap_rows)

    tests = []
    for a, b in (("PESTO", "PESTO-titles"), ("PESTO-titles", "Model knowledge"),
                 ("PESTO", "Model knowledge")):
        p = fig2d.paired_perm(vec[a], vec[b], y, np.random.default_rng(7))
        tests.append({"better": a, "worse": b, "p_one_sided": round(p, 5)})
        print(f"{a} > {b}: P = {p:.4f}")
    write(RESULTS / "ed2_titles_tests.tsv", tests)

    verdict_rows = []
    for m in METHODS:
        counts = Counter(verdict_of(scores[m][p]) for p in limited)
        for v in VERDICTS:
            verdict_rows.append({"method": m, "verdict": v, "n": counts.get(v, 0)})
    write(RESULTS / "ed2_titles_verdicts.tsv", verdict_rows)

    cost = fig2d.read_tsv(RESULTS / "fig2_ed5_cost.tsv")
    cost_rows = []
    for m, col in (("Model knowledge", "prior_usd"), ("PESTO-titles", "titles_usd"),
                   ("PESTO", "current_usd")):
        v = np.array([float(r[col]) for r in cost if r[col] not in ("", "NA")])
        cost_rows.append({"method": m, "n": len(v), "mean_usd": round(v.mean(), 4),
                          "sd_usd": round(v.std(ddof=1), 4)})
    write(RESULTS / "ed2_titles_cost.tsv", cost_rows)


if __name__ == "__main__":
    main()
