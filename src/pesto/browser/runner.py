"""Running a run from the browser: the same pipeline as `pesto --bench`.

Each pair goes through cli.assess, the function the command line uses, so a run
made here and one made with `pesto --bench` are the same answers in the same
tables. Pairs run in a thread pool; the tables and run.json are rewritten after
each one, so the page can follow, and a stopped or interrupted run resumes from
what it has.
"""
from __future__ import annotations

import contextlib
import io
import threading
import types
from concurrent.futures import ThreadPoolExecutor
import os

from .. import cli
from ..cost import load_table, write_table as write_cost
from ..flow import arms
from . import runs

ACTIVE = {}
_LOCK = threading.Lock()


class Job:
    def __init__(self, path):
        self.path = path
        self.stop = threading.Event()


def start(path, workers=8):
    with _LOCK:
        if path in ACTIVE:
            return ACTIVE[path]
        job = Job(path)
        ACTIVE[path] = job
    threading.Thread(target=_execute, args=(job, workers), daemon=True).start()
    return job


def stop(path):
    job = ACTIVE.get(path)
    if job:
        job.stop.set()


def _execute(job, workers):
    path = job.path
    meta = runs.read_meta(path)
    meta.pop("path", None)
    write_lock = threading.Lock()
    try:
        pairs = runs.read_input(path)
        family = "fulltext" if meta.get("mode") == "fulltext" else "abstracts"
        data_dir = os.path.join(path, "runs")
        flow_dir = os.path.join(data_dir, "flow")
        os.makedirs(flow_dir, exist_ok=True)
        args = types.SimpleNamespace(
            no_cache=False, data_dir=data_dir, ot_encoder=None, ot_top=None,
            no_ot_gate=False, ot_reasons=False, fulltext=family == "fulltext")
        extra = [k for k in pairs[0] if k not in ("gene", "phenotype")]
        fields = ["gene", "phenotype"] + extra + cli.ANSWER_FIELDS
        out = os.path.join(path, "pesto.tsv")
        costs_out = os.path.join(path, "cost.tsv")
        allowed = {arms.get(n).fingerprint() for n in arms.FAMILIES[family]}
        done = cli.resumable(out, allowed, fields)
        costs = load_table(costs_out)
        order = [(p["gene"], p["phenotype"]) for p in pairs]

        def save_meta(**changes):
            meta.update(changes)
            meta["done"] = len(done)
            meta["usd"] = round(sum(float(r["usd"]) for rows in costs.values()
                                    for r in rows if r.get("step") == "total"), 4)
            runs.write_meta(path, meta)

        todo = [p for p in pairs if (p["gene"], p["phenotype"]) not in done]
        save_meta(status="running", error="")
        if todo:
            # One load, before the pool, so twenty workers do not each fetch
            # the embedding model at once.
            from ..ot_shortlist import encoder, resolve_encoder
            if resolve_encoder(None) != "none":
                save_meta(status="loading the embedding model")
                encoder(None)
                save_meta(status="running")

        def work(p):
            if job.stop.is_set():
                return
            gene, pheno = p["gene"], p["phenotype"]
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    arm = arms.resolve(arms.AUTO, pheno, family=family)
                row = cli.assess(arm, gene, pheno, args, None, flow_dir)
            except Exception as exc:
                row = {k: "" for k in fields}
                row.update({"gene": gene, "phenotype": pheno,
                            "error": f"{type(exc).__name__}: {exc}"})
            for k in extra:
                row[k] = p.get(k, "")
            with write_lock:
                done[(gene, pheno)] = row
                if row.get("_cost"):
                    costs[(gene, pheno)] = row["_cost"]
                cli.write_table(out, fields, [done[k] for k in order if k in done])
                write_cost(costs_out, costs)
                save_meta()

        with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            list(pool.map(work, todo))
        finished = len(done) >= len(pairs)
        errors = sum(1 for r in done.values() if (r.get("error") or "").strip())
        save_meta(status="done" if finished else "stopped",
                  error=f"{errors} pair(s) failed" if errors else "")
    except Exception as exc:
        meta.update(status="error", error=f"{type(exc).__name__}: {exc}")
        runs.write_meta(path, meta)
    finally:
        with _LOCK:
            ACTIVE.pop(path, None)
