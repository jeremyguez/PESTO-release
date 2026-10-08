"""ROC of the returned label and the ASD-candidate false-Established rate.

Two things the four-colour bars of figure 2 cannot show. First, how well the
label separates the GenCC classes a curator called strong (Definitive,
Strong) from the ones they did not (Moderate, Limited): a ROC over the four
label levels, one curve per reading prompt, AUCs compared by DeLong on the
same 400 pairs. Second, how often each prompt calls Established a gene the
ASC curators listed only as a candidate, which cannot be Established by
construction: an exact rate with a paired mid-p McNemar against current.

Both are computed on the labels, not on the 100-point distribution, because
the label is what the pipeline returns. On the continuous score the arms are
indistinguishable, and that number is printed here so the choice is visible.

  PESTO_PROJECT_ROOT=$PWD python3 scripts/127_fig2_roc_candidate.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy import optimize
from scipy import stats as spstats

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("PESTO_PROJECT_ROOT", str(ROOT))

RESULTS = ROOT / "results"
ORDER = ("Novel", "Hypothesized", "Existing", "Established")
CALLRANK = {c: i + 1 for i, c in enumerate(ORDER)}
GENCC = ("Definitive", "Strong", "Moderate", "Limited")
STRONG = ("Definitive", "Strong")
WEAK = ("Moderate", "Limited")
ARMS = ("current", "score", "current-bare", "current-open")
ASD_CANDIDATE = {"ID known", "autism candidate"}
ROC_TSV = RESULTS / "fig2_roc_points.tsv"
AUC_TSV = RESULTS / "fig2_roc_auc.tsv"
CAND_TSV = RESULTS / "fig2_asd_candidate.tsv"


def load_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def as_call(x):
    """The mean_call of figure 2: round the 1-4 mean of the distribution."""
    return ORDER[max(1, min(4, round(float(x)))) - 1]


def lit4(dist):
    p = {k: float(dist.get(k) or dist.get(k.lower()) or 0) for k in ORDER}
    total = sum(p.values()) or 1
    return (4 * p["Established"] + 3 * p["Existing"]
            + 2 * p["Hypothesized"] + p["Novel"]) / total


def key_of(gene, phenotype):
    return gene.upper(), (phenotype or "").strip().lower()


# ------------------------------------------------------------------ loading

def gencc_labels():
    """The 400 curated GenCC pairs under each of the four prompts."""
    call, cont, cls = {a: {} for a in ARMS}, {a: {} for a in ARMS}, {}
    for r in load_tsv(RESULTS / "fig2a_extra40" / "combined_n100.tsv"):
        k = key_of(r["gene"], r["phenotype"])
        cls[k] = r["gencc_class"]
        call["current"][k] = as_call(r["current"])
        cont["current"][k] = float(r["current"])
        call["score"][k] = as_call(r["score"])
        cont["score"][k] = float(r["score"])
    for arm, name in (("current-bare", "bare"), ("current-open", "open")):
        path = RESULTS / f"arm_{arm}_fig2ab.json"
        for r in json.loads(path.read_text()):
            if r.get("panel") != "2a" or r.get("error") or not r.get("distribution"):
                continue
            if r["group"] not in GENCC:
                continue
            k = key_of(r["gene"], r["phenotype"])
            call[arm][k] = as_call(lit4(r["distribution"]))
            cont[arm][k] = lit4(r["distribution"])
    keys = sorted(set.intersection(*(set(call[a]) for a in ARMS)))
    return keys, cls, call, cont


def asd_candidate_labels():
    """The 56 ASC candidate genes: a curator's floor, not a grey zone."""
    call = {a: {} for a in ARMS}
    for r in load_tsv(RESULTS / "bench_asd_current.tsv"):
        if r.get("classification") not in ASD_CANDIDATE:
            continue
        if r["lit_mean"] not in CALLRANK:
            continue
        call["current"][key_of(r["gene"], r["phenotype"])] = r["lit_mean"]
    for r in load_tsv(RESULTS / "arm_current-score_fig2ab.tsv"):
        if r.get("panel") != "2b" or r.get("group") != "ASD candidate":
            continue
        if r.get("error") or not r.get("novelty"):
            continue
        k = key_of(r["gene"], r["phenotype"])
        call["score"][k] = as_call(1 + 3 * (100 - float(r["novelty"])) / 100)
    for arm in ("current-bare", "current-open"):
        path = RESULTS / f"arm_{arm}_fig2ab.json"
        for r in json.loads(path.read_text()):
            if r.get("panel") != "2b" or r.get("group") != "ASD candidate":
                continue
            if r.get("error") or not r.get("distribution"):
                continue
            call[arm][key_of(r["gene"], r["phenotype"])] = as_call(
                lit4(r["distribution"]))
    keys = sorted(set.intersection(*(set(call[a]) for a in ARMS)))
    return keys, call


# ------------------------------------------------------------------- ROC

def roc_points(scores, positive):
    """TPR and FPR at every cut of an ordinal score, plus both corners."""
    pos, neg = scores[positive], scores[~positive]
    pts = [(0.0, 0.0)]
    for cut in sorted(set(scores), reverse=True):
        pts.append((float(np.mean(neg >= cut)), float(np.mean(pos >= cut))))
    if pts[-1] != (1.0, 1.0):
        pts.append((1.0, 1.0))
    return pts


def psi(x, y):
    d = x[:, None] - y[None, :]
    return (d > 0) * 1.0 + (d == 0) * 0.5


def delong(scores, positive):
    """AUCs, their SEs and the covariance-aware test of Delong 1988."""
    m, n = int(positive.sum()), int((~positive).sum())
    aucs, v10, v01 = [], [], []
    for s in scores:
        p = psi(s[positive], s[~positive])
        aucs.append(float(p.mean()))
        v10.append(p.mean(axis=1))
        v01.append(p.mean(axis=0))
    cov = np.cov(np.array(v10), ddof=1) / m + np.cov(np.array(v01), ddof=1) / n
    cov = np.atleast_2d(cov)
    return np.array(aucs), np.sqrt(np.diag(cov)), cov


def delong_test(aucs, cov, i, j):
    var = cov[i, i] + cov[j, j] - 2 * cov[i, j]
    if var <= 0:
        return float("nan"), float("nan")
    z = (aucs[i] - aucs[j]) / np.sqrt(var)
    return z, 2 * spstats.norm.sf(abs(z))


# ---------------------------------------------------------------- proportions

def mid_p_mcnemar(n01, n10):
    """Exact conditional McNemar minus half its point mass.

    Fagerland, Lydersen & Laake 2013 measured the exact conditional test as
    over-conservative on matched pairs; mid-p holds its level and is what a
    handful of discordant pairs can support.
    """
    disc = n01 + n10
    if disc == 0:
        return 1.0, 1.0
    exact = spstats.binomtest(n10, disc, 0.5).pvalue
    return exact, exact - spstats.binom.pmf(n10, disc, 0.5)


def tango_ci(n01, n10, n, level=0.95):
    """Score interval for a paired difference in proportions.

    Returned for the other arm minus current, so that it brackets the excess
    rate reported next to it rather than its negative.
    """
    z = spstats.norm.ppf(1 - (1 - level) / 2)

    def score(delta):
        a = 2 * n
        b = -n01 - n10 + (2 * n - n01 + n10) * delta
        c = -n10 * delta * (1 - delta)
        disc = b * b - 4 * a * c
        if disc < 0:
            return np.inf
        q = (-b + np.sqrt(disc)) / (2 * a)
        var = n * (2 * q + delta * (1 - delta))
        if var <= 0:
            return np.inf
        return (n01 - n10 - n * delta) / np.sqrt(var)

    lo = optimize.brentq(lambda d: score(d) - z, -0.9999, 0.9999)
    hi = optimize.brentq(lambda d: score(d) + z, -0.9999, 0.9999)
    return -hi, -lo


def main():
    keys, cls, call, cont = gencc_labels()
    if len(keys) != 400:
        print(f"WARNING: {len(keys)} paired GenCC pairs, expected 400",
              file=sys.stderr)
    klass = np.array([cls[k] for k in keys])
    positive = np.isin(klass, STRONG)
    print(f"GenCC {len(keys)} pairs, {int(positive.sum())} Definitive+Strong "
          f"against {int((~positive).sum())} Moderate+Limited")

    labels = [np.array([CALLRANK[call[a][k]] for k in keys], float) for a in ARMS]
    aucs, ses, cov = delong(labels, positive)
    conts = [np.array([cont[a][k] for k in keys], float) for a in ARMS]
    caucs, _, _ = delong(conts, positive)

    rows, auc_rows = [], []
    for i, arm in enumerate(ARMS):
        z, p = delong_test(aucs, cov, 0, i) if i else (float("nan"), float("nan"))
        auc_rows.append({
            "arm": arm, "auc": round(aucs[i], 4), "se": round(ses[i], 4),
            "auc_lo": round(max(0, aucs[i] - 1.96 * ses[i]), 4),
            "auc_hi": round(min(1, aucs[i] + 1.96 * ses[i]), 4),
            "dauc_vs_current": round(aucs[0] - aucs[i], 4),
            "delong_z": round(z, 3), "delong_p": p,
            "auc_continuous": round(caucs[i], 4),
        })
        for fpr, tpr in roc_points(labels[i], positive):
            rows.append({"arm": arm, "fpr": round(fpr, 6), "tpr": round(tpr, 6)})
    write_tsv(ROC_TSV, rows)
    write_tsv(AUC_TSV, auc_rows)
    print("\nAUC of the label, Definitive+Strong vs Moderate+Limited")
    for r in auc_rows:
        tail = ("" if r["arm"] == "current"
                else f"  dAUC={r['dauc_vs_current']:+.3f}  "
                     f"DeLong Z={r['delong_z']:.2f} P={r['delong_p']:.3g}")
        print(f"  {r['arm']:13s} {r['auc']:.3f} (SE {r['se']:.3f})"
              f"   continuous {r['auc_continuous']:.3f}{tail}")

    ckeys, ccall = asd_candidate_labels()
    n = len(ckeys)
    print(f"\nASD candidate {n} genes, none of which can be Established")
    est = {a: np.array([ccall[a][k] == "Established" for k in ckeys]) for a in ARMS}
    cand_rows = []
    for arm in ARMS:
        hits = int(est[arm].sum())
        lo, hi = spstats.binomtest(hits, n).proportion_ci(method="exact")
        if arm == "current":
            exact = mid = float("nan")
            n01 = n10 = 0
            dlo = dhi = float("nan")
        else:
            n01 = int((est["current"] & ~est[arm]).sum())
            n10 = int((~est["current"] & est[arm]).sum())
            exact, mid = mid_p_mcnemar(n01, n10)
            dlo, dhi = tango_ci(n01, n10, n)
        cand_rows.append({
            "arm": arm, "n": n, "established": hits,
            "rate": round(100 * hits / n, 2),
            "rate_lo": round(100 * lo, 2), "rate_hi": round(100 * hi, 2),
            "lost_vs_current": n01, "gained_vs_current": n10,
            "mcnemar_exact_p": exact, "mcnemar_midp": mid,
            "diff_pts": round(100 * (hits - int(est["current"].sum())) / n, 2),
            "diff_lo": round(100 * dlo, 2), "diff_hi": round(100 * dhi, 2),
        })
        tail = ("" if arm == "current" else
                f"   +{cand_rows[-1]['diff_pts']:.1f} pts "
                f"[{cand_rows[-1]['diff_lo']:.1f}, {cand_rows[-1]['diff_hi']:.1f}]"
                f"  mid-p={mid:.4g} (exact {exact:.4g})")
        print(f"  {arm:13s} {hits:2d}/{n}  {100*hits/n:5.1f}% "
              f"[{100*lo:.1f}, {100*hi:.1f}]{tail}")
    write_tsv(CAND_TSV, cand_rows)
    print(f"\nwrote {ROC_TSV.relative_to(ROOT)}, {AUC_TSV.relative_to(ROOT)} "
          f"and {CAND_TSV.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
