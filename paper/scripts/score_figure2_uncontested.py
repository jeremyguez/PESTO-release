#!/usr/bin/env python3
"""Combine the current-arm run and the Open Targets run into one table.

Literature is the rounded mean of the hundred points, not the argmax: a split
45/40/12/3 is Existing, which is what figure 2a already plotted. Open Targets
returns one band. `both` is the more established of the two, Novel being the
weakest claim.

Usage:
  python3 scripts/score_figure2_uncontested.py
"""
from __future__ import annotations

import csv
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.join(ROOT, "benchmark", "gencc_uncontested")
LIT = os.path.join(HERE, "results", "arm_current.json")
OT = os.path.join(HERE, "results", "ot.tsv")
BENCH = os.path.join(HERE, "pairs.tsv")
OUT = os.path.join(HERE, "results", "scored.tsv")

VERDICTS = ["Novel", "Hypothesized", "Existing", "Established"]
WEIGHT = {"Established": 4, "Existing": 3, "Hypothesized": 2, "Novel": 1}


def mean_call(dist):
    if not dist:
        return "", None
    # Keys as stored by Distribution.as_dict.
    pts = {k.capitalize(): dist.get(k) or dist.get(k.lower()) or 0
           for k in ("established", "existing", "hypothesized", "novel")}
    # Tolerate already-capitalised keys.
    for k in WEIGHT:
        if k in dist:
            pts[k] = dist[k]
    total = sum(pts.values()) or 1
    rank = sum(WEIGHT[k] * pts[k] for k in WEIGHT) / total
    idx = max(1, min(4, round(rank)))
    return VERDICTS[idx - 1], round(rank, 3)


def stronger(a, b):
    if a not in WEIGHT:
        return b if b in WEIGHT else ""
    if b not in WEIGHT:
        return a
    return a if WEIGHT[a] >= WEIGHT[b] else b


lit = {}
if os.path.exists(LIT):
    for r in json.load(open(LIT, encoding="utf-8")):
        call, rank = mean_call(r.get("distribution") or {})
        lit[(r["gene"], r["phenotype"])] = {
            "lit_argmax": r.get("verdict") or "",
            "lit_mean": call,
            "lit_rank": "" if rank is None else rank,
            "lit_error": r.get("error") or "",
            "read": r.get("read") or "",
        }

ot = {}
if os.path.exists(OT):
    for r in csv.DictReader(open(OT, encoding="utf-8"), delimiter="\t"):
        ot[(r["gene"], r["phenotype"])] = r

fields = ["gene", "phenotype", "gencc_class",
          "lit_mean", "lit_argmax", "lit_rank",
          "ot_verdict", "ot_score", "both",
          "read", "lit_error", "ot_error"]
rows = []
for r in csv.DictReader(open(BENCH, encoding="utf-8"), delimiter="\t"):
    key = (r["gene"], r["phenotype"])
    L, O = lit.get(key, {}), ot.get(key, {})
    lm, ov = L.get("lit_mean", ""), O.get("ot_verdict", "")
    rows.append({
        "gene": r["gene"], "phenotype": r["phenotype"],
        "gencc_class": r["gencc_class"],
        "lit_mean": lm, "lit_argmax": L.get("lit_argmax", ""),
        "lit_rank": L.get("lit_rank", ""),
        "ot_verdict": ov, "ot_score": O.get("ot_score", ""),
        "both": stronger(lm, ov),
        "read": L.get("read", ""),
        "lit_error": L.get("lit_error", ""),
        "ot_error": O.get("error", ""),
    })

with open(OUT, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

n_lit = sum(1 for r in rows if r["lit_mean"])
n_ot = sum(1 for r in rows if r["ot_verdict"] and r["ot_verdict"] != "Error")
print(f"wrote {os.path.relpath(OUT, ROOT)}, {n_lit} literature, {n_ot} Open Targets")
