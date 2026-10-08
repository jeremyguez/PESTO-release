"""Prior + current on 20 extra Figure 2a-style random Absents.

Does not run titles. Writes results/fig2a_absent_extra20/.

Usage:
  python3 scripts/86_score_fig2a_absent20.py
"""
from __future__ import annotations

import csv
import os
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("PESTO_PROJECT_ROOT", str(ROOT))

from pesto.cost import capturing, from_trace, load_table, with_total, write_table
from pesto.flow import arms, run as flow_run
from pesto.flow.steps import Argmax
from pesto.flow.types import Query
from pesto.harness import MODELS
from pesto.services.llm_service import call_llm_with_usage
from pesto.utils.helpers import load_prompt

SRC = ROOT / "results" / "fig2a_absent_extra20.tsv"
OUT = ROOT / "results" / "fig2a_absent_extra20"
WORKERS = 20
PRIOR_FP = "2b18abb989c5"
CURRENT_FP = "ff136250f372"
PROMPT = load_prompt("knowledge_probabilities")
LOCK = threading.Lock()
KNOW_FIELDS = [
    "gene", "phenotype", "verdict", "lit_score",
    "p_established", "p_existing", "p_hypothesized", "p_novel",
    "justification", "input_tokens", "output_tokens", "usd", "calls",
]
FLOW_FIELDS = [
    "gene", "phenotype", "verdict", "lit_score",
    "found", "read", "cached", "error", "usd", "bands_usd",
]


def lit_score(dist):
    if not dist:
        return None
    return (4 * dist.get("Established", 0) + 3 * dist.get("Existing", 0)
            + 2 * dist.get("Hypothesized", 0) + 1 * dist.get("Novel", 0)) / 100


def load_tsv(path):
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def save_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def step_usd(spent, name):
    return sum(float(r["usd"]) for r in (spent or []) if r.get("step") == name)


def pairs():
    out = []
    for r in load_tsv(SRC):
        out.append({"gene": r["gene"], "phenotype": r["phenotype"]})
    return out


def knowledge_one(gene, phenotype):
    text = PROMPT.format(gene_name=gene, phenotype=phenotype)
    with capturing() as cap:
        resp = call_llm_with_usage(
            text, MODELS["opus5"], 0.0, agent_name="knowledge_prior") or {}
    raw = resp.get("text", "") or ""
    verdict = Argmax().run(Query(gene, phenotype), raw)
    spent = cap.summary(MODELS["opus5"])
    dist = verdict.distribution.as_dict() if verdict.distribution else {}
    return {
        "gene": gene, "phenotype": phenotype,
        "verdict": verdict.call, "lit_score": lit_score(dist),
        "p_established": dist.get("Established", ""),
        "p_existing": dist.get("Existing", ""),
        "p_hypothesized": dist.get("Hypothesized", ""),
        "p_novel": dist.get("Novel", ""),
        "justification": (verdict.justification or "").replace("\t", " "),
        **spent,
    }


def run_knowledge(todo_pairs):
    path = OUT / "knowledge.tsv"
    rows = load_tsv(path)
    have = {(r["gene"], r["phenotype"]) for r in rows
            if r.get("verdict") in ("Novel", "Hypothesized", "Existing",
                                    "Established")}
    todo = [p for p in todo_pairs if (p["gene"], p["phenotype"]) not in have]
    print(f"prior   {len(todo_pairs) - len(todo)} done, {len(todo)} to run",
          flush=True)
    n_new = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = {pool.submit(knowledge_one, p["gene"], p["phenotype"]): p
                for p in todo}
        for fut in as_completed(futs):
            try:
                one = fut.result()
            except Exception as exc:
                p = futs[fut]
                one = {"gene": p["gene"], "phenotype": p["phenotype"],
                       "verdict": "", "lit_score": "", "error": str(exc)}
                print(f"  FAIL {p['gene']}: {exc}", flush=True)
            with LOCK:
                rows.append(one)
                save_tsv(path, rows, KNOW_FIELDS)
                n_new += 1
                usd = one.get("usd") or 0
                print(f"  prior {n_new}/{len(todo)}  {one['gene']}  "
                      f"{one.get('verdict')}  ${float(usd or 0):.3f}",
                      flush=True)
    return rows


def run_current(todo_pairs):
    path = OUT / "pesto.tsv"
    cost_path = str(OUT / "pesto.cost.tsv")
    runs_dir = str(OUT / "runs" / "flow")
    rows = load_tsv(path)
    have = {(r["gene"], r["phenotype"]) for r in rows if not r.get("error")}
    todo = [p for p in todo_pairs if (p["gene"], p["phenotype"]) not in have]
    arm = arms.get("current")
    print(f"current {len(todo_pairs) - len(todo)} done, {len(todo)} to run  "
          f"[{arm.fingerprint()}]", flush=True)
    cost = load_table(cost_path)

    def work(p):
        query = Query(p["gene"], p["phenotype"])
        try:
            res = flow_run.run(arm, query, cache=True, runs_dir=runs_dir)
            dist = (res.verdict.distribution.as_dict()
                    if res.verdict.distribution else {})
            spent = with_total(from_trace(p["gene"], p["phenotype"],
                                          res.trace, res.cached),
                               p["gene"], p["phenotype"], res.cached)
            return {
                "gene": p["gene"], "phenotype": p["phenotype"],
                "verdict": res.verdict.call, "lit_score": lit_score(dist),
                "found": res.counts.get("found", ""),
                "read": res.counts.get("read", ""),
                "cached": res.cached, "error": "",
                "usd": f"{step_usd(spent, 'total'):.6f}",
                "bands_usd": f"{step_usd(spent, 'Bands'):.6f}",
            }, spent
        except Exception as exc:
            return {
                "gene": p["gene"], "phenotype": p["phenotype"],
                "verdict": "", "lit_score": "",
                "found": "", "read": "", "cached": False,
                "error": f"{type(exc).__name__}: {exc}",
                "usd": "", "bands_usd": "",
            }, None

    n_new = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = [pool.submit(work, p) for p in todo]
        for fut in as_completed(futs):
            rec, spent = fut.result()
            with LOCK:
                rows.append(rec)
                if spent is not None:
                    cost[(rec["gene"], rec["phenotype"])] = spent
                save_tsv(path, rows, FLOW_FIELDS)
                write_table(cost_path, cost)
                n_new += 1
                print(f"  current {n_new}/{len(todo)}  {rec['gene']}  "
                      f"{rec.get('verdict') or rec.get('error', '')[:40]}  "
                      f"${rec.get('usd') or 0}", flush=True)
    return rows


def main():
    prior, current = arms.get("prior"), arms.get("current")
    if prior.fingerprint() != PRIOR_FP:
        print(f"ERROR: prior fingerprint {prior.fingerprint()}", file=sys.stderr)
        return 1
    if current.fingerprint() != CURRENT_FP:
        print(f"ERROR: current fingerprint {current.fingerprint()}",
              file=sys.stderr)
        return 1
    todo = pairs()
    print(f"n={len(todo)}  {WORKERS} workers  → {OUT}", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    know = run_knowledge(todo)
    pesto = run_current(todo)
    n_err = sum(1 for r in pesto if r.get("error"))
    n_k = sum(1 for r in know if r.get("verdict") in
              ("Novel", "Hypothesized", "Existing", "Established"))
    print(f"prior ok {n_k}/{len(todo)}  pesto errors {n_err}/{len(todo)}",
          flush=True)
    return 1 if n_err or n_k < len(todo) else 0


if __name__ == "__main__":
    sys.exit(main())
