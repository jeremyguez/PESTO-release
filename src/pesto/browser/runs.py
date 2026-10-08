"""Runs on disk, as the browser sees them.

A run is a folder. It holds the same files `pesto --bench` writes, so a folder
made on the command line opens in the browser and the other way round:

  input.tsv   the pairs asked about (gene, phenotype, and any other column)
  pesto.tsv   one row of answers per pair, rewritten as pairs finish
  cost.tsv    what each pair cost, step by step
  runs/       the evidence: Open Targets records, and runs/flow/ with one JSON
              per pair holding the corpus, the prompts' output and the verdict
  run.json    the browser's own note: name, date, mode, progress

Runs made in the browser go under the runs folder (settings); a folder opened
from elsewhere is remembered in the settings file and read where it is.
"""
from __future__ import annotations

import csv
import datetime
import json
import os
import re

from .. import config
from ..flow.store import slug
from ..flow.types import Query

_CONFIG_HOME = (os.environ.get("XDG_CONFIG_HOME")
                or os.path.join(os.path.expanduser("~"), ".config"))
SETTINGS = os.path.join(_CONFIG_HOME, "pesto", "browser.json")
DEFAULT_RUNS_DIR = os.path.join(config.PROJECT_ROOT, "runs")
BANDS = {4: "strong human", 3: "weak human", 2: "non-human", 1: "no bearing"}


# ------------------------------------------------------------------ settings

def load_settings():
    try:
        with open(SETTINGS, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def save_settings(settings):
    os.makedirs(os.path.dirname(SETTINGS), exist_ok=True)
    with open(SETTINGS, "w", encoding="utf-8") as fh:
        json.dump(settings, fh, indent=1)


def runs_dir():
    return load_settings().get("runs_dir") or DEFAULT_RUNS_DIR


def set_runs_dir(path):
    path = os.path.abspath(os.path.expanduser(path))
    os.makedirs(path, exist_ok=True)
    s = load_settings()
    s["runs_dir"] = path
    save_settings(s)
    return path


def remember(path):
    s = load_settings()
    opened = [p for p in s.get("opened", []) if p != path]
    s["opened"] = [path] + opened
    save_settings(s)


def forget(path):
    s = load_settings()
    s["opened"] = [p for p in s.get("opened", []) if p != path]
    save_settings(s)


# ---------------------------------------------------------------- run folders

def is_run(path):
    return os.path.isdir(path) and any(
        os.path.exists(os.path.join(path, f))
        for f in ("run.json", "pesto.tsv", "input.tsv"))


def allowed(path):
    """Only runs under the runs folder, or opened by the user, are served."""
    real = os.path.realpath(path)
    root = os.path.realpath(runs_dir())
    if real.startswith(root + os.sep):
        return True
    return real in {os.path.realpath(p) for p in load_settings().get("opened", [])}


def _tsv(path):
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def create(name, pairs, mode):
    """A new run folder with its input written, ready to be run."""
    stamp = datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    clean = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")[:50] or "run"
    path = os.path.join(runs_dir(), f"{stamp}_{clean}")
    os.makedirs(os.path.join(path, "runs", "flow"), exist_ok=True)
    cols = ["gene", "phenotype"] + [k for k in pairs[0] if k not in ("gene", "phenotype")]
    with open(os.path.join(path, "input.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(pairs)
    write_meta(path, {"name": name, "created": datetime.datetime.now().isoformat(
        timespec="seconds"), "pairs": len(pairs), "mode": mode, "status": "queued",
        "done": 0, "usd": 0.0, "version": _version()})
    return path


def _version():
    try:
        from importlib.metadata import version
        return version("pesto-genetics")
    except Exception:
        return ""


def read_input(path):
    rows = _tsv(os.path.join(path, "input.tsv"))
    if rows:
        return rows
    return [{"gene": r["gene"], "phenotype": r["phenotype"]}
            for r in _tsv(os.path.join(path, "pesto.tsv"))]


def write_meta(path, meta):
    tmp = os.path.join(path, "run.json.part")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=1)
    os.replace(tmp, os.path.join(path, "run.json"))


def read_meta(path):
    try:
        with open(os.path.join(path, "run.json"), encoding="utf-8") as fh:
            meta = json.load(fh)
    except (OSError, ValueError):
        # A folder written by `pesto --bench`: everything else is read off it.
        rows = _tsv(os.path.join(path, "pesto.tsv"))
        meta = {"name": os.path.basename(os.path.normpath(path)),
                "created": datetime.datetime.fromtimestamp(
                    os.path.getmtime(path)).isoformat(timespec="seconds"),
                "pairs": len(read_input(path)), "mode": "",
                "status": "done", "done": len(rows), "usd": total_cost(path)}
    meta["path"] = os.path.abspath(path)
    return meta


def total_cost(path):
    return round(sum(float(r.get("usd") or 0)
                     for r in _tsv(os.path.join(path, "cost.tsv"))
                     if r.get("step") == "total"), 4)


def list_runs():
    out, seen = [], set()
    root = runs_dir()
    paths = []
    if os.path.isdir(root):
        paths += [os.path.join(root, d) for d in os.listdir(root)]
    paths += load_settings().get("opened", [])
    for p in paths:
        real = os.path.realpath(p)
        if real in seen or not is_run(p):
            continue
        seen.add(real)
        out.append(read_meta(p))
    return sorted(out, key=lambda m: m.get("created", ""), reverse=True)


def rows(path):
    """The input pairs, each with its answer if it has one."""
    answered = {(r["gene"], r["phenotype"]): r
                for r in _tsv(os.path.join(path, "pesto.tsv"))}
    costs = {}
    for r in _tsv(os.path.join(path, "cost.tsv")):
        if r.get("step") == "total":
            costs[(r["gene"], r["phenotype"])] = float(r.get("usd") or 0)
    out = []
    for p in read_input(path):
        key = (p["gene"], p["phenotype"])
        r = dict(answered.get(key) or p)
        r["answered"] = key in answered
        r["usd"] = costs.get(key)
        r.pop("justification", None)
        out.append(r)
    return out


def pair(path, gene, phenotype):
    """Everything saved about one pair: the answer row, the literature run,
    the Open Targets record and the cost."""
    row = next((r for r in _tsv(os.path.join(path, "pesto.tsv"))
                if r["gene"] == gene and r["phenotype"] == phenotype), None)
    out = {"row": row, "flow": None, "ot": None, "cost": [],
           "bands": BANDS}
    if row and row.get("fingerprint"):
        f = os.path.join(path, "runs", "flow",
                         f"{slug(Query(gene, phenotype))}__{row['fingerprint']}.json")
        if os.path.exists(f):
            with open(f, encoding="utf-8") as fh:
                flow = json.load(fh)
            flow.pop("raw", None)
            flow.pop("trace", None)
            labels = flow.get("labels") or {}
            for a in flow.get("articles") or []:
                lab = labels.get(str(a.get("pmid"))) or {}
                a["band"] = lab.get("relevance")
            flow["corpus_size"] = len(labels)
            out["flow"] = flow
    from ..open_targets import _ot_match_sidecar_path
    side = _ot_match_sidecar_path(os.path.join(path, "runs"), gene, phenotype)
    if os.path.exists(side):
        with open(side, encoding="utf-8") as fh:
            rec = json.load(fh)
        out["ot"] = {"ensembl_id": rec.get("ensembl_id"),
                     "of": rec.get("of") or len(rec.get("traits") or []),
                     "graded": rec.get("graded") or [],
                     "dropped": rec.get("gate_dropped_ids") or [],
                     "traits": rec.get("traits") or []}
    out["cost"] = [r for r in _tsv(os.path.join(path, "cost.tsv"))
                   if r["gene"] == gene and r["phenotype"] == phenotype]
    return out
