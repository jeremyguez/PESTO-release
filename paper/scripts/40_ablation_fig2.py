"""The two figure-2 ablations on every pair in panels a and b.

  knowledge  one Opus call, pair only, no articles
  titles     current, abstracts never fetched

Resumes a partial file. Reuses the ten-pair pilot cache.
"""
from __future__ import annotations

import csv
import os
import shutil
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("PESTO_PROJECT_ROOT", str(ROOT))

from pesto.cost import capturing, from_trace, row, with_total, write_table
from pesto.flow import arms, run as flow_run
from pesto.flow.steps import Argmax
from pesto.flow.types import Query
from pesto.harness import MODELS
from pesto.services.llm_service import call_llm_with_usage
from pesto.utils.helpers import load_prompt

OUT = ROOT / "results" / "ablation_fig2"
PILOT = ROOT / "results" / "ablation_pilot10"
N_KNOWLEDGE = 1
WORKERS = 20
PROMPT = load_prompt("knowledge_probabilities")
LOCK = threading.Lock()


def lit_score(dist):
    if not dist:
        return None
    return (4 * dist.get("Established", 0) + 3 * dist.get("Existing", 0)
            + 2 * dist.get("Hypothesized", 0) + 1 * dist.get("Novel", 0)) / 100


def write_pairs():
    rows = []
    with (ROOT / "benchmark" / "gencc_clingen_g2p" / "pesto.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("gencc_class") in {
                "Definitive", "Strong", "Moderate", "Limited", "Absent"}:
                rows.append({"gene": r["gene"], "phenotype": r["phenotype"],
                             "set": f"gencc_{r['gencc_class'].lower()}"})
    curated = []
    with (ROOT / "benchmark" / "asd" / "pesto.tsv").open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            curated.append(r)
            rows.append({"gene": r["gene"], "phenotype": r["phenotype"],
                         "set": "asd_curated"})
    have = {r["gene"].upper() for r in curated}
    extra101 = []
    with (ROOT / "benchmark" / "asd_extra101" / "pesto.tsv").open() as fh:
        extra101 = [r for r in csv.DictReader(fh, delimiter="\t")
                    if r.get("set") == "random"]
    extra101.sort(key=lambda r: r["gene"].upper())
    extras = extra101[:50]
    with (ROOT / "benchmark" / "extra50_random" / "pesto.tsv").open() as fh:
        extras += [r for r in csv.DictReader(fh, delimiter="\t")
                   if r.get("phenotype") == "Autism Spectrum Disorder"]
    for r in extras:
        if r["gene"].upper() in have:
            continue
        have.add(r["gene"].upper())
        rows.append({"gene": r["gene"], "phenotype": r["phenotype"],
                     "set": "not_in_asc"})
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "pairs.tsv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, ["gene", "phenotype", "set"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return rows


def seed_from_pilot(titles):
    dest = OUT / "runs" / "titles"
    dest.mkdir(parents=True, exist_ok=True)
    src = PILOT / "runs" / "titles"
    if src.is_dir():
        for f in src.glob(f"*__{titles.fingerprint()}.json"):
            target = dest / f.name
            if not target.exists():
                shutil.copy2(f, target)
    know_src = PILOT / "knowledge.tsv"
    know_dst = OUT / "knowledge.tsv"
    if know_src.exists() and not know_dst.exists():
        shutil.copy2(know_src, know_dst)
    titles_src = PILOT / "titles.tsv"
    titles_dst = OUT / "titles.tsv"
    if titles_src.exists() and not titles_dst.exists():
        shutil.copy2(titles_src, titles_dst)


def load_tsv(path):
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def save_tsv(path, rows, fields):
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


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
        "gene": gene, "phenotype": phenotype, "verdict": verdict.call,
        "lit_score": lit_score(dist),
        "p_established": dist.get("Established", ""),
        "p_existing": dist.get("Existing", ""),
        "p_hypothesized": dist.get("Hypothesized", ""),
        "p_novel": dist.get("Novel", ""),
        "justification": (verdict.justification or "").replace("\t", " "),
        **spent,
    }


KNOW_FIELDS = [
    "gene", "phenotype", "verdict", "lit_score", "p_established",
    "p_existing", "p_hypothesized", "p_novel", "justification",
    "input_tokens", "output_tokens", "usd", "calls", "model", "draw",
]


def run_knowledge(pairs):
    path = OUT / "knowledge.tsv"
    rows = load_tsv(path)
    done = {}
    for r in rows:
        done.setdefault((r["gene"], r["phenotype"]), set()).add(int(r["draw"]))
    todo = [p for p in pairs
            if len(done.get((p["gene"], p["phenotype"]), ())) < N_KNOWLEDGE]
    print(f"knowledge  {len(pairs) - len(todo)} done, {len(todo)} to run",
          flush=True)

    def work(p):
        gene, pheno = p["gene"], p["phenotype"]
        have = done.get((gene, pheno), set())
        out = []
        for k in range(1, N_KNOWLEDGE + 1):
            if k in have:
                continue
            one = knowledge_one(gene, pheno)
            one["draw"] = k
            out.append(one)
        return out

    n_new = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futs = [pool.submit(work, p) for p in todo]
        for fut in as_completed(futs):
            added = fut.result()
            with LOCK:
                rows.extend(added)
                save_tsv(path, rows, KNOW_FIELDS)
                n_new += len(added)
                print(f"  knowledge {n_new}/{len(todo) * N_KNOWLEDGE} new  "
                      f"{added[-1]['gene'] if added else ''}", flush=True)
    return rows


def run_titles(pairs, titles):
    path = OUT / "titles.tsv"
    runs_dir = str(OUT / "runs" / "titles")
    rows = load_tsv(path)
    have = {(r["gene"], r["phenotype"]) for r in rows if not r.get("error")}
    todo = [p for p in pairs if (p["gene"], p["phenotype"]) not in have]
    print(f"titles     {len(pairs) - len(todo)} done, {len(todo)} to run",
          flush=True)
    fields = ["gene", "phenotype", "set", "verdict", "lit_score",
              "found", "read", "cached", "error"]
    cost = {}

    def work(p):
        gene, pheno = p["gene"], p["phenotype"]
        query = Query(gene, pheno)
        try:
            res = flow_run.run(titles, query, cache=True, runs_dir=runs_dir)
            dist = (res.verdict.distribution.as_dict()
                    if res.verdict.distribution else {})
            return {
                "gene": gene, "phenotype": pheno, "set": p.get("set", ""),
                "verdict": res.verdict.call, "lit_score": lit_score(dist),
                "found": res.counts.get("found", ""),
                "read": res.counts.get("read", ""),
                "cached": res.cached, "error": "",
            }, with_total(from_trace(gene, pheno, res.trace, res.cached),
                          gene, pheno, res.cached)
        except Exception as exc:
            return {
                "gene": gene, "phenotype": pheno, "set": p.get("set", ""),
                "verdict": "", "lit_score": "",
                "found": "", "read": "", "cached": False,
                "error": f"{type(exc).__name__}: {exc}",
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
                save_tsv(path, rows, fields)
                n_new += 1
                print(f"  titles {n_new}/{len(todo)}  {rec['gene']}  "
                      f"{rec.get('verdict') or rec.get('error', '')[:40]}",
                      flush=True)
    if cost:
        write_table(str(OUT / "titles.cost.tsv"), cost)
    return rows


def summarize(know, titles_rows):
    by = {}
    for r in know:
        by.setdefault((r["gene"], r["phenotype"]), []).append(r)
    k_usd = sum(float(r.get("usd") or 0) for r in know)
    print()
    print(f"knowledge  {len(by)} pairs, {len(know)} calls, ${k_usd:.2f}")
    print(f"titles     {sum(1 for r in titles_rows if not r.get('error'))} ok, "
          f"{sum(1 for r in titles_rows if r.get('error'))} errors")


def main():
    current = arms.get("current")
    titles = arms.get("titles")
    print(f"current [{current.fingerprint()}]", flush=True)
    print(f"titles  [{titles.fingerprint()}]", flush=True)
    if current.fingerprint() != "ff136250f372":
        print("ERROR: current fingerprint moved; aborting", file=sys.stderr)
        return 1
    pairs = write_pairs()
    print(f"{len(pairs)} pairs", flush=True)
    seed_from_pilot(titles)
    know = run_knowledge(pairs)
    titles_rows = run_titles(pairs, titles)
    summarize(know, titles_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
