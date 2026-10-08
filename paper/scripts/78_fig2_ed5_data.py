"""Pair costs and Novel rates for Extended Data 2.

Three methods: current, titles (thinking on), prior. The Novel rate is taken
over every pair of figure 2 a class holds present in a reference: A not
Absent, B not « not in ASC », and the forty Limited pairs of Fig. 2d
(2023–2025 first-report bin). The costs stay on A+B, the only pairs whose
runs are reconstructible here. Writes results/fig2_ed5_novel.tsv and
results/fig2_ed5_cost.tsv.
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

from pesto.config import price_usd
from pesto.flow.store import slug
from pesto.flow.types import Query
from pesto.harness import BATCH, DISMISS
from pesto.services.pubmed_service import format_articles_for_llm
from pesto.utils.helpers import load_prompt

ROOT = Path(__file__).resolve().parents[1]
AB = ROOT / "results" / "thinking_off_titles_ab.tsv"
YEAR10 = ROOT / "results" / "prior_year10"
KNOW = ROOT / "results" / "ablation_fig2" / "knowledge.tsv"
TITLES_COST = ROOT / "results" / "ablation_fig2" / "titles.cost.tsv"
PESTO_A = ROOT / "benchmark" / "gencc_clingen_g2p" / "pesto.tsv"
PESTO_B = ROOT / "results" / "bench_asd_current.tsv"
TITLES_RUNS = ROOT / "results" / "ablation_fig2" / "runs" / "titles"
FLOW_DIRS = [
    ROOT / "benchmark" / "gencc_clingen_g2p" / "runs" / "flow",
    ROOT / "benchmark" / "asd" / "runs" / "flow",
    ROOT / "benchmark" / "asd_extra101" / "runs" / "flow",
    ROOT / "benchmark" / "extra50_random" / "runs" / "flow",
]
E50_COST = ROOT / "benchmark" / "extra50_random" / "cost.tsv"
E50_RUNS = ROOT / "benchmark" / "extra50_random" / "runs" / "flow"
FP = "ff136250f372"
TITLES_FP = "1ec6267b715a"
OK = {"Novel", "Hypothesized", "Existing", "Established"}
BOOT = 100_000
SEED = 2026
OUT_NOVEL = ROOT / "results" / "fig2_ed5_novel.tsv"
OUT_TEST = ROOT / "results" / "fig2_ed5_novel_tests.tsv"
OUT_COST = ROOT / "results" / "fig2_ed5_cost.tsv"


def load(path):
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def key(gene, phenotype):
    return (gene.upper().strip(), (phenotype or "").strip().lower())


def year_bin(year):
    y = int(year)
    if y <= 2010:
        return None
    if y <= 2013:
        return "2011-2013"
    if y <= 2016:
        return "2014-2016"
    if y <= 2019:
        return "2017-2019"
    if y <= 2022:
        return "2020-2022"
    return "2023-2025"


def load_2d_limited():
    """Forty Limited pairs of Fig. 2d, 2023–2025 bin."""
    know = load(YEAR10 / "knowledge.tsv")
    pesto = {key(r["gene"], r["phenotype"]): r
             for r in json.loads((YEAR10 / "pesto.json").read_text())}
    titles = {key(r["gene"], r["phenotype"]): r
              for r in json.loads(
                  (YEAR10 / "titles_2023_2025.json").read_text())}
    rows = []
    for r in know:
        if not r.get("year") or year_bin(r["year"]) != "2023-2025":
            continue
        k = key(r["gene"], r["phenotype"])
        p, t = pesto.get(k), titles.get(k)
        if not p or not t:
            raise SystemExit(f"missing 2d scores for {r['gene']}")
        cur, tit, pri = p.get("verdict"), t.get("verdict"), r.get("verdict")
        if cur not in OK or tit not in OK or pri not in OK:
            raise SystemExit(f"bad verdict {r['gene']}")
        rows.append({"current": cur, "titles": tit, "prior": pri,
                     "panel": "D", "group": "Limited 2023-2025",
                     "gene": r["gene"], "phenotype": r["phenotype"]})
    if len(rows) != 40:
        raise SystemExit(f"expected 40 Limited 2023-2025, got {len(rows)}")
    return rows


def run_path(gene, phenotype, folder, fp):
    return folder / f"{slug(Query(gene, phenotype))}__{fp}.json"


def load_run(gene, phenotype, folders, fp):
    for folder in folders:
        path = run_path(gene, phenotype, folder, fp)
        if path.exists():
            return json.loads(path.read_text())
    return None


def articles(run):
    arts = run.get("articles") or []
    if isinstance(arts, dict):
        arts = list(arts.values())
    return arts


def selected(run):
    arts = articles(run)
    n = int((run.get("counts") or {}).get("read") or 0)
    if n and len(arts) >= n:
        return arts[:n]
    return arts


def pr_prompt(run, abstracts):
    from pesto.flow.steps import ProbabilityRead
    tpl = ProbabilityRead().template()
    return tpl.format(
        gene_name=run["gene"], phenotype=run["phenotype"],
        articles_list=format_articles_for_llm(selected(run),
                                             use_abstracts=abstracts),
        alias_note="", dismiss_note=DISMISS)


def band_prompts(run, abstracts):
    tpl = load_prompt("evidence_band_prompt")
    syns = (run.get("terms") or {}).get("synonyms") or []
    names = "\n".join(f"  {s}" for s in syns) or f"  {run['phenotype']}"
    arts = articles(run)
    found = int((run.get("counts") or {}).get("found") or 0)
    scale = max(found, len(arts), 1) / max(len(arts), 1)
    out = []
    for i in range(0, len(arts), BATCH):
        chunk = arts[i:i + BATCH]
        out.append(tpl.format(
            gene_name=run["gene"], phenotype=run["phenotype"],
            synonym_note=names,
            articles_list=format_articles_for_llm(chunk,
                                                 use_abstracts=abstracts)))
    return out, scale, max(found, len(arts), 1)


def calibrate():
    e50 = load(E50_COST)
    in_rats, out_rats = [], []
    for r in e50:
        if r["step"] != "ProbabilityRead":
            continue
        run = load_run(r["gene"], r["phenotype"], [E50_RUNS], FP)
        if not run:
            continue
        prompt = pr_prompt(run, True)
        raw = run.get("raw") or ""
        inn, out = int(r["input_tokens"]), int(r["output_tokens"])
        if inn and prompt:
            in_rats.append(inn / len(prompt))
        if out and raw:
            out_rats.append(out / max(len(raw), 1))
    return (float(np.median(in_rats)) if in_rats else 0.25,
            float(np.median(out_rats)) if out_rats else 1.5)


def tokens(texts, rat):
    return sum(len(t) for t in texts) * rat


def mcnemar_p(only_cur, only_oth):
    n = only_cur + only_oth
    if n == 0:
        return 1.0
    return float(binomtest(only_oth, n, 0.5, alternative="greater").pvalue)


def boot_p(cur, oth, rng):
    n = len(cur)
    d0 = oth.mean() - cur.mean()
    more = 0
    for _ in range(BOOT):
        ix = rng.integers(0, n, n)
        if oth[ix].mean() - cur[ix].mean() <= 0:
            more += 1
    return d0, more / BOOT


def main():
    ab = load(AB)
    pesto = {}
    for path, col in (
        (PESTO_A, "lit_mean"),
        (PESTO_B, "lit_mean"),
        (ROOT / "benchmark" / "asd" / "pesto.tsv", "lit_mean"),
        (ROOT / "benchmark" / "asd_extra101" / "pesto.tsv", "lit_mean"),
        (ROOT / "benchmark" / "extra50_random" / "pesto.tsv", "lit_mean"),
    ):
        for r in load(path):
            pesto[key(r["gene"], r["phenotype"])] = r[col]
    prior = {}
    for r in load(KNOW):
        if str(r.get("draw", "1")) == "1":
            prior[key(r["gene"], r["phenotype"])] = (r["verdict"], float(r["usd"]))
    tcost = defaultdict(dict)
    for r in load(TITLES_COST):
        tcost[key(r["gene"], r["phenotype"])][r["step"]] = r

    in_rat, out_rat = calibrate()
    print(f"opus PR calibration  in {in_rat:.3f}  out {out_rat:.3f} tok/char")

    present, costs = [], []
    n_cur_run = n_tit_run = 0
    for r in ab:
        k = key(r["gene"], r["phenotype"])
        cur = pesto[k]
        tit = r["titles_verdict"]
        pri, pri_usd = prior[k]
        if cur not in OK or tit not in OK or pri not in OK:
            raise SystemExit(f"bad verdict {r['gene']}")
        if r["group"] not in {"Absent", "not in ASC"}:
            present.append({"current": cur, "titles": tit, "prior": pri,
                            "panel": r["panel"], "group": r["group"],
                            "gene": r["gene"], "phenotype": r["phenotype"]})

        steps = tcost.get(k, {})
        syn = float(steps.get("Synonyms", {}).get("usd") or 0)
        br = float(steps.get("Broadened", {}).get("usd") or 0)
        tit_pr = float(steps.get("ProbabilityRead", {}).get("usd") or 0)

        trun = load_run(r["gene"], r["phenotype"], [TITLES_RUNS], TITLES_FP)
        crun = load_run(r["gene"], r["phenotype"], FLOW_DIRS, FP)
        tit_band = 0.0
        if trun:
            n_tit_run += 1
            prompts, scale, n_found = band_prompts(trun, False)
            tit_band = price_usd("haiku", tokens(prompts, in_rat) * scale,
                                 15 * n_found)
        cur_usd = syn + br
        if crun:
            n_cur_run += 1
            prompts, scale, n_found = band_prompts(crun, True)
            cur_usd += price_usd("haiku", tokens(prompts, in_rat) * scale,
                                 15 * n_found)
            prompt = pr_prompt(crun, True)
            raw = crun.get("raw") or ""
            cur_usd += price_usd("opus5", len(prompt) * in_rat,
                                 max(len(raw), 1) * out_rat)
        costs.append({
            "gene": r["gene"], "phenotype": r["phenotype"],
            "panel": r["panel"],
            "prior_usd": f"{pri_usd:.6f}",
            "titles_usd": f"{syn + br + tit_pr + tit_band:.6f}",
            "current_usd": f"{cur_usd:.6f}" if crun else "",
        })

    n_ab_present = len(present)
    extra = load_2d_limited()
    present.extend(extra)

    print(f"A+B {len(ab)}  present {n_ab_present} + D "
          f"{len(extra)} = {len(present)}  "
          f"titles runs {n_tit_run}  current runs {n_cur_run}")

    rng = np.random.default_rng(SEED)
    novel_rows, test_rows = [], []
    print(f"{'method':<10} novel/n")
    for method in ("current", "titles", "prior"):
        y = np.array([p[method] == "Novel" for p in present], dtype=float)
        novel_rows.append({
            "method": method, "n": len(present),
            "novel": int(y.sum()), "rate": f"{y.mean():.4f}",
        })
        print(f"  {method:<10} {int(y.sum())}/{len(present)} {y.mean():.1%}")
    cur = np.array([p["current"] == "Novel" for p in present], dtype=float)
    for other in ("titles", "prior"):
        oth = np.array([p[other] == "Novel" for p in present], dtype=float)
        only_c = int(((cur == 1) & (oth == 0)).sum())
        only_o = int(((cur == 0) & (oth == 1)).sum())
        delta, p_boot = boot_p(cur, oth, rng)
        test_rows.append({
            "vs": other, "delta": f"{delta:+.4f}",
            "only_current": only_c, "only_other": only_o,
            "mcnemar_p": f"{mcnemar_p(only_c, only_o):.4g}",
            "bootstrap_p": f"{p_boot:.4g}",
        })
        print(f"  current vs {other}: Δ {delta:+.3f}  "
              f"McNemar {mcnemar_p(only_c, only_o):.4g}  "
              f"boot {p_boot:.4g}")

    with OUT_NOVEL.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["method", "n", "novel", "rate"], delimiter="\t")
        w.writeheader()
        w.writerows(novel_rows)
    with OUT_TEST.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["vs", "delta", "only_current", "only_other",
                                "mcnemar_p", "bootstrap_p"], delimiter="\t")
        w.writeheader()
        w.writerows(test_rows)
    with OUT_COST.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["gene", "phenotype", "panel",
                                "prior_usd", "titles_usd", "current_usd"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(costs)
    for name, col in (("prior", "prior_usd"), ("titles", "titles_usd"),
                      ("current", "current_usd")):
        xs = [float(r[col]) for r in costs if r[col]]
        print(f"  {name:<10} n={len(xs)}  mean ${np.mean(xs):.3f}  "
              f"sd ${np.std(xs, ddof=1):.3f}")
    print(f"wrote {OUT_NOVEL.name}, {OUT_TEST.name}, {OUT_COST.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
