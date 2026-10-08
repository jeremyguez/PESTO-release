"""Add 30 GenCC pairs per class to Figure 2a and re-read them two ways.

Keeps the original 120 curated pairs. Draws 30 more Definitive / Strong /
Moderate / Limited from the same ClinGen+G2P uncontested pool, spending genes
already used in those 120 and in the 100 Absent of Figure 2a. Seed 2027 so the
extra 30 are a new shuffle, not the tail of seed 2026.

  python3 scripts/123_fig2a_extra30.py --draw
  PESTO_PROJECT_ROOT=$PWD python3 scripts/run_arm.py --arm current \
      --bench results/fig2a_extra30/pairs.tsv --workers 20 \
      --runs-dir results/fig2a_extra30/runs/flow \
      --out results/fig2a_extra30/pesto.json
  PESTO_PROJECT_ROOT=$PWD   python3 scripts/123_fig2a_extra30.py --draw40
  PESTO_PROJECT_ROOT=$PWD python3 scripts/run_arm.py --arm current \
      --bench results/fig2a_extra40/pairs.tsv --workers 20 \
      --runs-dir results/fig2a_extra40/runs/flow \
      --out results/fig2a_extra40/pesto.json
  PESTO_PROJECT_ROOT=$PWD python3 scripts/123_fig2a_extra30.py --score40
  PESTO_PROJECT_ROOT=$PWD python3 scripts/123_fig2a_extra30.py --stats
  Rscript scripts/122_figure2_current_vs_score.R
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import threading
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("PESTO_PROJECT_ROOT", str(ROOT))

ORIG = ROOT / "benchmark" / "gencc_clingen_g2p" / "pairs.tsv"
ORIG_PESTO = ROOT / "benchmark" / "gencc_clingen_g2p" / "pesto.tsv"
KEEP = ROOT / "results" / "bench_figure2_clingen_g2p_keep.tsv"
N60 = ROOT / "results" / "bench_figure2_clingen_g2p_n60.tsv"
OUT = ROOT / "results" / "fig2a_extra30"
PAIRS = OUT / "pairs.tsv"
SCORE_JSON = OUT / "score.json"
SCORE_TSV = OUT / "score.tsv"
COMBINED = OUT / "combined_n60.tsv"
OUT40 = ROOT / "results" / "fig2a_extra40"
KEEP100 = ROOT / "results" / "bench_figure2_clingen_g2p_keep100.tsv"
N100 = ROOT / "results" / "bench_figure2_clingen_g2p_n100.tsv"
COMBINED100 = ROOT / "results" / "fig2a_extra40" / "combined_n100.tsv"
CURRENT_FP = "ff136250f372"
CLASSES = ["Definitive", "Strong", "Moderate", "Limited"]
WORKERS = 20
NOVELTY = re.compile(
    r"<novelty>\s*([0-9]+(?:\.[0-9]+)?)\s*</novelty>", re.I)
LOCK = threading.Lock()
RANK = {"Limited": 1, "Moderate": 2, "Strong": 3, "Definitive": 4}


def load_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def write_keep():
    rows = load_tsv(ORIG)
    have = {(r["gene"].upper(), r["phenotype"].strip().lower()) for r in rows}
    fields = ["gene", "phenotype", "gencc_class", "disease_curie", "submitters",
              "n_submitters", "max_cosine", "nearest_curated_disease"]
    extra = []
    for folder in ("fig2a_absent_extra20", "fig2a_absent_extra20b",
                   "fig2a_absent_extra30"):
        for r in load_tsv(ROOT / "results" / folder / "pesto.tsv"):
            key = (r["gene"].upper(), r["phenotype"].strip().lower())
            if key in have or r.get("error"):
                continue
            have.add(key)
            extra.append({
                "gene": r["gene"], "phenotype": r["phenotype"],
                "gencc_class": "Absent", "disease_curie": "",
                "submitters": "", "n_submitters": "0",
                "max_cosine": "", "nearest_curated_disease": "",
            })
    keep = [{k: r.get(k, "") for k in fields} for r in rows] + extra
    write_tsv(KEEP, keep, fields)
    n_abs = sum(1 for r in keep if r["gencc_class"] == "Absent")
    n_cur = len(keep) - n_abs
    print(f"keep {n_cur} curated + {n_abs} absent -> {KEEP.relative_to(ROOT)}")
    return n_abs


def draw():
    n_abs = write_keep()
    cmd = [
        sys.executable, str(ROOT / "scripts" / "draw_figure2_benchmark.py"),
        "--n", "60", "--absent", str(n_abs),
        "--submitters", "ClinGen,G2P",
        "--keep", str(KEEP),
        "--seed", "2027",
        "--out", str(N60),
        "--funnel", str(ROOT / "results" / "bench_figure2_clingen_g2p_n60_funnel.tsv"),
    ]
    print(" ".join(cmd), flush=True)
    subprocess.check_call(cmd, cwd=ROOT)
    orig_genes = {r["gene"].upper() for r in load_tsv(ORIG)
                  if r["gencc_class"] in CLASSES}
    orig_keys = {(r["gene"].upper(), r["phenotype"].strip().lower())
                 for r in load_tsv(ORIG) if r["gencc_class"] in CLASSES}
    n60 = load_tsv(N60)
    pool, dropped = [], []
    kept_orig = 0
    for r in n60:
        if r["gencc_class"] not in CLASSES:
            continue
        key = (r["gene"].upper(), r["phenotype"].strip().lower())
        if key in orig_keys:
            kept_orig += 1
            continue
        if r["gene"].upper() in orig_genes:
            continue
        pool.append(r)
    for r in load_tsv(ORIG):
        if r["gencc_class"] not in CLASSES:
            continue
        key = (r["gene"].upper(), r["phenotype"].strip().lower())
        in_n60 = any(
            (x["gene"].upper(), x["phenotype"].strip().lower()) == key
            and x["gencc_class"] in CLASSES for x in n60)
        if not in_n60:
            dropped.append(f"{r['gene']} / {r['phenotype'][:50]} ({r['gencc_class']})")
    rng = __import__("random").Random(2027)
    extra = []
    by = defaultdict(list)
    for r in pool:
        by[r["gencc_class"]].append(r)
    for cls in CLASSES:
        cand = sorted(by[cls], key=lambda r: (r["gene"].upper(), r["phenotype"]))
        rng.shuffle(cand)
        extra.extend(cand[:30])
        print(f"  {cls:12s} pool {len(cand):3d}  take {min(30, len(cand))}")
    counts = Counter(r["gencc_class"] for r in extra)
    write_tsv(PAIRS, extra, ["gene", "phenotype", "gencc_class", "disease_curie",
                             "submitters", "n_submitters"])
    print(f"original curated still in n60: {kept_orig}/120")
    print(f"new pairs {len(extra)}  {dict(counts)}")
    if dropped:
        print(f"dropped from keep ({len(dropped)}); originals stay in the analysis")
    print(f"wrote {PAIRS.relative_to(ROOT)}")
    if any(counts[c] != 30 for c in CLASSES):
        return 1
    return 0


def spent_genes():
    genes = set()
    for r in load_tsv(ORIG):
        genes.add(r["gene"].upper())
    for r in load_tsv(PAIRS):
        genes.add(r["gene"].upper())
    for folder in ("fig2a_absent_extra20", "fig2a_absent_extra20b",
                   "fig2a_absent_extra30"):
        for r in load_tsv(ROOT / "results" / folder / "pesto.tsv"):
            if r.get("gene"):
                genes.add(r["gene"].upper())
    return genes


def write_keep100():
    fields = ["gene", "phenotype", "gencc_class", "disease_curie", "submitters",
              "n_submitters", "max_cosine", "nearest_curated_disease"]
    rows = []
    have = set()

    def add(r, cls=None):
        key = (r["gene"].upper(), r["phenotype"].strip().lower())
        if key in have:
            return
        have.add(key)
        rows.append({
            "gene": r["gene"], "phenotype": r["phenotype"],
            "gencc_class": cls or r["gencc_class"],
            "disease_curie": r.get("disease_curie", ""),
            "submitters": r.get("submitters", ""),
            "n_submitters": r.get("n_submitters", ""),
            "max_cosine": r.get("max_cosine", ""),
            "nearest_curated_disease": r.get("nearest_curated_disease", ""),
        })

    for r in load_tsv(ORIG):
        add(r)
    for r in load_tsv(PAIRS):
        add(r)
    for folder in ("fig2a_absent_extra20", "fig2a_absent_extra20b",
                   "fig2a_absent_extra30"):
        for r in load_tsv(ROOT / "results" / folder / "pesto.tsv"):
            if r.get("error"):
                continue
            add(r, "Absent")
    write_tsv(KEEP100, rows, fields)
    n_abs = sum(1 for r in rows if r["gencc_class"] == "Absent")
    n_cur = len(rows) - n_abs
    print(f"keep {n_cur} curated + {n_abs} absent -> {KEEP100.relative_to(ROOT)}")
    return n_abs


def draw40():
    n_abs = write_keep100()
    cmd = [
        sys.executable, str(ROOT / "scripts" / "draw_figure2_benchmark.py"),
        "--n", "100", "--absent", str(n_abs),
        "--submitters", "ClinGen,G2P",
        "--keep", str(KEEP100),
        "--seed", "2028",
        "--out", str(N100),
        "--funnel", str(ROOT / "results" / "bench_figure2_clingen_g2p_n100_funnel.tsv"),
    ]
    print(" ".join(cmd), flush=True)
    subprocess.check_call(cmd, cwd=ROOT)
    spent = spent_genes()
    n100 = load_tsv(N100)
    pool = []
    for r in n100:
        if r["gencc_class"] not in CLASSES:
            continue
        if r["gene"].upper() in spent:
            continue
        pool.append(r)
    rng = __import__("random").Random(2028)
    extra = []
    by = defaultdict(list)
    for r in pool:
        by[r["gencc_class"]].append(r)
    for cls in CLASSES:
        cand = sorted(by[cls], key=lambda r: (r["gene"].upper(), r["phenotype"]))
        rng.shuffle(cand)
        extra.extend(cand[:40])
        print(f"  {cls:12s} pool {len(cand):3d}  take {min(40, len(cand))}")
    counts = Counter(r["gencc_class"] for r in extra)
    write_tsv(OUT40 / "pairs.tsv", extra,
              ["gene", "phenotype", "gencc_class", "disease_curie",
               "submitters", "n_submitters"])
    print(f"new pairs {len(extra)}  {dict(counts)}")
    print(f"wrote {(OUT40 / 'pairs.tsv').relative_to(ROOT)}")
    if any(counts[c] != 40 for c in CLASSES):
        return 1
    return 0


def parse_novelty(text):
    m = NOVELTY.search(text or "")
    if not m:
        return None
    value = float(m.group(1))
    if value < 0 or value > 100:
        return None
    return value


def lit4_from_p(r):
    return (4 * float(r["p_established"]) + 3 * float(r["p_existing"])
            + 2 * float(r["p_hypothesized"]) + 1 * float(r["p_novel"])) / 100


def lit4_from_dist(d):
    if not d:
        return None
    return (4 * float(d.get("Established") or d.get("established") or 0)
            + 3 * float(d.get("Existing") or d.get("existing") or 0)
            + 2 * float(d.get("Hypothesized") or d.get("hypothesized") or 0)
            + 1 * float(d.get("Novel") or d.get("novel") or 0)) / 100


def score4(nov):
    return 1 + 3 * (100 - float(nov)) / 100


def score(out_dir=None):
    from pesto.flow import arms
    from pesto.flow.steps import NoveltyScoreRead
    from pesto.flow.store import slug
    from pesto.flow.types import Corpus, Query, FIELDS as ARTICLE_FIELDS

    current = arms.get("current")
    if current.fingerprint() != CURRENT_FP:
        print(f"ERROR: current fingerprint moved to {current.fingerprint()}",
              file=sys.stderr)
        return 1
    folder = Path(out_dir) if out_dir else OUT
    pairs = folder / "pairs.tsv"
    score_json = folder / "score.json"
    score_tsv = folder / "score.tsv"
    rows = load_tsv(pairs)
    flow_dir = folder / "runs" / "flow"
    have = json.loads(score_json.read_text()) if score_json.exists() else []
    done_keys = {(r["gene"].upper(), r["phenotype"].strip().lower())
                 for r in have if r.get("novelty") is not None and not r.get("error")}
    todo = [r for r in rows
            if (r["gene"].upper(), r["phenotype"].strip().lower()) not in done_keys]
    print(f"{len(rows)} extra pairs, {len(done_keys)} scored, {len(todo)} to run",
          flush=True)
    reader = NoveltyScoreRead("opus5")

    def work(row):
        query = Query(row["gene"], row["phenotype"])
        out = {"gene": query.gene, "phenotype": query.phenotype,
               "gencc_class": row["gencc_class"]}
        path = flow_dir / f"{slug(query)}__{CURRENT_FP}.json"
        try:
            raw = json.loads(path.read_text())
            corpus = Corpus.of(raw.get("articles") or [], ARTICLE_FIELDS)
            text, _ = reader.run(query, corpus)
            novelty = parse_novelty(text)
            if novelty is None:
                text, _ = reader.run(query, corpus)
                novelty = parse_novelty(text)
            if novelty is None:
                out["error"] = "PARSE_FAIL"
                out["read"] = len(corpus)
                return out
            out.update({"novelty": novelty, "evidence": 100.0 - novelty,
                        "read": len(corpus)})
        except Exception as exc:
            out["error"] = f"{type(exc).__name__}: {exc}"
        return out

    done = list(have)
    started = time.time()
    n_new = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = [pool.submit(work, r) for r in todo]
        for fut in as_completed(futs):
            rec = fut.result()
            with LOCK:
                done.append(rec)
                score_json.write_text(json.dumps(done, indent=1))
                write_tsv(score_tsv, done, ["gene", "phenotype", "gencc_class",
                                            "novelty", "evidence", "read", "error"])
                n_new += 1
                shown = rec.get("error") or f"n={rec.get('novelty')}"
                print(f"  {n_new:3d}/{len(todo)}  {rec['gene']:10s} {shown:16s}  "
                      f"{time.time()-started:5.0f}s", flush=True)
    n_bad = sum(1 for r in done if r.get("error"))
    print(f"errors {n_bad}/{len(done)}")
    return 1 if n_bad else 0


def orig_current():
    pesto = {(r["gene"].upper(), r["phenotype"].strip().lower()): r
             for r in load_tsv(ORIG_PESTO)}
    out = []
    for r in load_tsv(ORIG):
        if r["gencc_class"] not in CLASSES:
            continue
        key = (r["gene"].upper(), r["phenotype"].strip().lower())
        p = pesto[key]
        out.append({
            "gene": r["gene"], "phenotype": r["phenotype"],
            "gencc_class": r["gencc_class"], "source": "original",
            "current": lit4_from_p(p),
        })
    return out


def extra_current():
    out = []
    for folder, source in ((OUT, "extra30"), (OUT40, "extra40")):
        pesto_json = folder / "pesto.json"
        pairs = folder / "pairs.tsv"
        if not pesto_json.exists() or not pairs.exists():
            continue
        meta = {(r["gene"].upper(), r["phenotype"].strip().lower()): r
                for r in load_tsv(pairs)}
        for r in json.loads(pesto_json.read_text()):
            key = (r["gene"].upper(), r["phenotype"].strip().lower())
            if r.get("error") or key not in meta:
                continue
            out.append({
                "gene": r["gene"], "phenotype": r["phenotype"],
                "gencc_class": meta[key]["gencc_class"], "source": source,
                "current": lit4_from_dist(r.get("distribution") or {}),
            })
    return out


def orig_score():
    out = []
    for r in load_tsv(ROOT / "results" / "arm_current-score_fig2ab.tsv"):
        if r.get("panel") != "2a" or r.get("group") not in CLASSES:
            continue
        if r.get("error") or r.get("novelty") in ("", None):
            continue
        out.append({
            "gene": r["gene"], "phenotype": r["phenotype"],
            "gencc_class": r["group"], "source": "original",
            "score": score4(r["novelty"]), "novelty": float(r["novelty"]),
        })
    return out


def extra_score():
    out = []
    for folder, source in ((OUT, "extra30"), (OUT40, "extra40")):
        path = folder / "score.tsv"
        if not path.exists():
            continue
        for r in load_tsv(path):
            if r.get("error") or r.get("novelty") in ("", None):
                continue
            out.append({
                "gene": r["gene"], "phenotype": r["phenotype"],
                "gencc_class": r["gencc_class"], "source": source,
                "score": score4(r["novelty"]), "novelty": float(r["novelty"]),
            })
    return out


def cliffs(a, b):
    a, b = list(a), list(b)
    gt = sum(x > y for x in a for y in b)
    lt = sum(x < y for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def perm_gap(dA, dB, nperm=20000, rng=None):
    import numpy as np
    rng = rng or np.random.default_rng(1)
    obs = dA.mean() - dB.mean()
    epsA = rng.choice([-1.0, 1.0], (nperm, len(dA)))
    epsB = rng.choice([-1.0, 1.0], (nperm, len(dB)))
    stat = (epsA * dA).mean(1) - (epsB * dB).mean(1)
    p1 = (np.sum(stat >= obs) + 1) / (nperm + 1)
    p2 = (np.sum(np.abs(stat) >= abs(obs)) + 1) / (nperm + 1)
    return obs, p1, p2


def perm_spearman(rk, cu, sc, nperm=20000, rng=None):
    import numpy as np
    from scipy import stats
    rng = rng or np.random.default_rng(1)

    def spear(x, y):
        return stats.spearmanr(x, y).statistic

    obs = spear(rk, cu) - spear(rk, sc)
    h1 = h2 = 0
    for _ in range(nperm):
        flip = rng.random(len(rk)) < 0.5
        a = np.where(flip, sc, cu)
        b = np.where(flip, cu, sc)
        d = spear(rk, a) - spear(rk, b)
        h1 += d >= obs
        h2 += abs(d) >= abs(obs)
    return obs, (h1 + 1) / (nperm + 1), (h2 + 1) / (nperm + 1)


def stats_main():
    import numpy as np
    from scipy import stats as spstats

    cur = orig_current() + extra_current()
    sco = orig_score() + extra_score()
    sco_map = {(r["gene"].upper(), r["phenotype"].strip().lower()): r for r in sco}
    rows = []
    for r in cur:
        key = (r["gene"].upper(), r["phenotype"].strip().lower())
        s = sco_map.get(key)
        if s is None:
            continue
        rows.append({
            "gene": r["gene"], "phenotype": r["phenotype"],
            "gencc_class": r["gencc_class"], "source": r["source"],
            "current": r["current"], "score": s["score"],
            "novelty": s["novelty"],
        })
    write_tsv(COMBINED, [r for r in rows if r["source"] != "extra40"])
    write_tsv(COMBINED100, rows)
    ORDER = ("Novel", "Hypothesized", "Existing", "Established")

    def as_verdict(x):
        return ORDER[max(1, min(4, round(x))) - 1]

    for source, folder in (("extra30", OUT), ("extra40", OUT40)):
        calls = [{"gene": r["gene"], "phenotype": r["phenotype"],
                  "gencc_class": r["gencc_class"],
                  "lit_mean": as_verdict(r["current"]),
                  "score_mean": as_verdict(r["score"]),
                  "novelty": r["novelty"]}
                 for r in rows if r["source"] == source]
        if calls:
            write_tsv(folder / "current_calls.tsv", calls)
    print(f"combined {len(rows)}  {Counter((r['gencc_class'], r['source']) for r in rows)}")
    by = defaultdict(lambda: {"current": [], "score": []})
    for r in rows:
        by[r["gencc_class"]]["current"].append(r["current"])
        by[r["gencc_class"]]["score"].append(r["score"])
    print("\nmeans")
    for cls in CLASSES:
        c = np.array(by[cls]["current"])
        s = np.array(by[cls]["score"])
        print(f"  {cls:12s} n={len(c):3d}  current {c.mean():.3f}  score {s.mean():.3f}")
    pairs = [("Definitive", "Strong"), ("Strong", "Moderate"), ("Moderate", "Limited")]
    print("\nadjacent current Δ > score Δ")
    for a, b in pairs:
        ca, cb = np.array(by[a]["current"]), np.array(by[b]["current"])
        sa, sb = np.array(by[a]["score"]), np.array(by[b]["score"])
        dA = ca - sa
        dB = cb - sb
        gap, p1, p2 = perm_gap(dA, dB)
        mw_c = spstats.mannwhitneyu(ca, cb, alternative="greater", method="asymptotic")
        mw_s = spstats.mannwhitneyu(sa, sb, alternative="greater", method="asymptotic")
        print(f"  {a[:3]}-{b[:3]}  curΔ={ca.mean()-cb.mean():.3f} (P1={mw_c.pvalue:.3g})  "
              f"scoΔ={sa.mean()-sb.mean():.3f} (P1={mw_s.pvalue:.3g})  "
              f"gap={gap:.3f}  P1={p1:.4g} P2={p2:.4g}  "
              f"δcur={cliffs(ca,cb):.3f} δsco={cliffs(sa,sb):.3f}")
    top = np.concatenate([by[c]["current"] for c in ("Definitive", "Strong", "Moderate")])
    top_s = np.concatenate([by[c]["score"] for c in ("Definitive", "Strong", "Moderate")])
    lim, lim_s = np.array(by["Limited"]["current"]), np.array(by["Limited"]["score"])
    gap, p1, p2 = perm_gap(top - top_s, lim - lim_s)
    print(f"\nDef+Strong+Mod vs Limited  curΔ={top.mean()-lim.mean():.3f}  "
          f"scoΔ={top_s.mean()-lim_s.mean():.3f}  gap={gap:.3f}  "
          f"P1={p1:.4g} P2={p2:.4g}")
    rk = np.array([RANK[r["gencc_class"]] for r in rows])
    cu = np.array([r["current"] for r in rows])
    sc = np.array([r["score"] for r in rows])
    rho_c = spstats.spearmanr(rk, cu)
    rho_s = spstats.spearmanr(rk, sc)
    d, p1, p2 = perm_spearman(rk, cu, sc)
    print(f"\nSpearman  current ρ={rho_c.statistic:.3f} (P={rho_c.pvalue:.3g})  "
          f"score ρ={rho_s.statistic:.3f} (P={rho_s.pvalue:.3g})  "
          f"Δρ={d:.3f}  P1={p1:.4g} P2={p2:.4g}")
    print(f"wrote {COMBINED.relative_to(ROOT)}")
    if (OUT40 / "pesto.json").exists():
        print(f"wrote {COMBINED100.relative_to(ROOT)}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--draw", action="store_true")
    ap.add_argument("--draw40", action="store_true")
    ap.add_argument("--score", action="store_true")
    ap.add_argument("--score40", action="store_true")
    ap.add_argument("--stats", action="store_true")
    args = ap.parse_args()
    if args.draw:
        return draw()
    if args.draw40:
        return draw40()
    if args.score:
        return score(OUT)
    if args.score40:
        return score(OUT40)
    if args.stats:
        return stats_main()
    ap.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
