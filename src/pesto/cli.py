#!/usr/bin/env python3
"""Ask whether a gene and a phenotype are already linked in the literature.

    pesto --gene DDX41 --phenotype "myelodysplastic syndrome"
    pesto --bench pairs.tsv
    pesto browser

Two answers come back, from two sources that are kept apart on purpose. The
literature verdict is read off the abstracts the pipeline retrieves and weighs.
The Open Targets verdict is read off that database's own evidence for the gene.
They disagree often, and the disagreement is informative, so neither is folded
into the other. A table of pairs is the same two answers, in parallel, one
process, so the embedding model is loaded once.
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import logging
import os
import sys
import textwrap
import threading
import time
from concurrent.futures import ThreadPoolExecutor

from . import config
from .credentials import anthropic_api_key

CREDENTIALS_SOURCE = ("https://github.com/jeremyguez/PESTO-release/blob/main/"
                      "src/pesto/credentials.py")

KEY_HELP = """\
No Anthropic API key found. PESTO reads the literature with Claude, so it needs
one, from console.anthropic.com:

  export ANTHROPIC_API_KEY=sk-ant-...

or the same line, without `export`, in a .env file you create here:
  {env_path}

The key stays on this machine. PESTO never writes it to disk, never prints it
and never stores it in a saved run; it is handed only to the official Anthropic
client, which sends it to api.anthropic.com. The code that reads it is short:
  {source}
"""


def pick_models(reader, worker):
    """What to read with, where None means: whatever the arm already says.

    An arm names the models its published numbers were produced by. Passing
    --model or --model-fast swaps them, and the fingerprint records that it did.
    """
    if not anthropic_api_key():
        raise SystemExit(KEY_HELP.format(
            env_path=os.path.join(config.PROJECT_ROOT, ".env"),
            source=CREDENTIALS_SOURCE))
    return reader or None, worker or None


def open_targets(gene, phenotype, synonyms, model, runs_dir, args=None):
    """That database's own answer, kept apart from the literature's.

    It is a lookup and a match, not a reading, so it does not belong in the arm.
    Where the two disagree is informative, which is the reason neither is folded
    into the other.
    """
    from .open_targets import fetch_and_save_open_targets_traits
    from .scoring import open_targets_verdict
    traits = fetch_and_save_open_targets_traits(
        gene, runs_dir, phenotype, list(synonyms or ()), model=model,
        matcher=getattr(args, "ot_matcher", None),
        encoder=getattr(args, "ot_encoder", None),
        top_k=getattr(args, "ot_top", None),
        # The flags read the other way round from the parameters, so that both
        # read naturally where they are written: `--no-ot-gate` on the command
        # line, `gate=False` in the call.
        gate=False if getattr(args, "no_ot_gate", False) else None,
        reasons=True if getattr(args, "ot_reasons", False) else None)
    match = (traits or {}).get("match") or {}
    ot_tokens = match.get("tokens") or match.get("usage") or {}
    return open_targets_verdict(traits) + (ot_tokens, match)


def render(result, ot_verdict, arm, show_reasoning):
    """The two verdicts, and what they were read from."""
    q = result.query
    out = [f"{q.gene}  /  {q.phenotype}", ""]
    out.append(f"  literature      {result.verdict.call}")
    out.append(f"  open targets    {ot_verdict or 'unknown'}")
    out.append("")

    read = result.counts.get("read", 0)
    found = result.counts.get("found", 0)
    source = ("reused an identical earlier run" if result.cached else
              f"{found} articles found, {read} read")
    out.append(f"  {source}, pipeline {arm.name} [{result.fingerprint}]")

    if show_reasoning and result.verdict.justification:
        out.append("")
        out.extend("  " + line for line in
                   textwrap.wrap(result.verdict.justification, width=76))
    if show_reasoning and result.verdict.pmids:
        out.append("")
        out.append("  rests on  " + ", ".join(result.verdict.pmids))
    return "\n".join(out)


WEIGHT = {"Established": 4, "Existing": 3, "Hypothesized": 2, "Novel": 1}
ORDER = ("Novel", "Hypothesized", "Existing", "Established")


def mean_call(dist):
    """The category the hundred points average to, as figure 2a already plots."""
    if not dist:
        return ""
    total = sum(dist.get(k, 0) for k in WEIGHT) or 1
    rank = sum(WEIGHT[k] * dist.get(k, 0) for k in WEIGHT) / total
    return ORDER[max(1, min(4, round(rank))) - 1]


def read_bench(path):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit(f"pesto: {path} is empty")
    if "gene" not in rows[0] or "phenotype" not in rows[0]:
        raise SystemExit(f"pesto: {path} needs gene and phenotype columns")
    return rows


def assess(arm, gene, phenotype, args, worker, flow_dir):
    """One pair: literature, then Open Targets. Same path as a single call."""
    from .flow import run as flow_run
    from .flow.types import Query
    query = Query(gene, phenotype)
    chatter = io.StringIO()
    with contextlib.redirect_stdout(chatter):
        result = flow_run.run(arm, query, cache=not args.no_cache,
                              runs_dir=flow_dir)
        ot_verdict, ot_score, ot_tokens, ot_match = open_targets(
            gene, phenotype, result.terms.get("synonyms"),
            worker or config.DEFAULT_MODEL, args.data_dir, args)
    dist = (result.verdict.distribution.as_dict()
            if result.verdict.distribution else None)
    from .cost import from_ot_tokens, from_trace, with_total
    return {
        "gene": gene, "phenotype": phenotype,
        "arm": arm.name, "fingerprint": result.fingerprint,
        "lit_argmax": result.verdict.call,
        "lit_mean": mean_call(dist or {}),
        "p_established": (dist or {}).get("Established", ""),
        "p_existing": (dist or {}).get("Existing", ""),
        "p_hypothesized": (dist or {}).get("Hypothesized", ""),
        "p_novel": (dist or {}).get("Novel", ""),
        "open_targets_verdict": ot_verdict,
        "open_targets_max_score": "" if ot_score is None else ot_score,
        "ot_channel": ot_match.get("channel") or "",
        "ot_tag": ot_match.get("deciding_tag") or "",
        "ot_trait": ot_match.get("basis") or "",
        "ot_cap": ot_match.get("cap") or "",
        "found": result.counts.get("found", 0),
        "read": result.counts.get("read", 0),
        "cached": result.cached,
        "justification": (result.verdict.justification or "").replace("\t", " "),
        "error": "",
        "_cost": with_total(
            from_trace(gene, phenotype, result.trace, result.cached)
            + from_ot_tokens(gene, phenotype, ot_tokens,
                             cached=bool(ot_match.get("cached"))),
            gene, phenotype, result.cached),
    }


def resumable(out, fingerprints, fields):
    """Pairs an earlier attempt already answered, keyed like the work queue.

    The criterion is the arm's fingerprint, exactly as in the flow cache: a row
    answers this question only if this pipeline produced it. Two rows are
    refused on top of that. One carrying an error is not an answer, and one
    without a literature verdict came from a run that stopped in the middle --
    which is precisely what an interrupted batch leaves behind.

    An empty Open Targets verdict is *not* refused: a gene that branch cannot
    resolve legitimately has none, and asking again would never fill it.

    Unlike the flow cache this reuses the summary row rather than reopening the
    stored run, so it also skips the Open Targets call, which has no cache of
    its own and is otherwise repaid on every attempt.
    """
    if not os.path.exists(out):
        return {}
    keep = {}
    if hasattr(fingerprints, "fingerprint"):
        allowed = {fingerprints.fingerprint()}
    else:
        allowed = set(fingerprints)
    with open(out, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("fingerprint") not in allowed:
                continue
            if (row.get("error") or "").strip():
                continue
            if not (row.get("lit_mean") or "").strip():
                continue
            keep[(row["gene"], row["phenotype"])] = {k: row.get(k, "")
                                                     for k in fields}
    return keep


def write_table(out, fields, rows):
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# The columns a table of answers carries after the input's own.
ANSWER_FIELDS = ["lit_mean", "lit_argmax", "p_established", "p_existing",
                 "p_hypothesized", "p_novel", "open_targets_verdict",
                 "open_targets_max_score", "ot_channel", "ot_tag",
                 "ot_trait", "ot_cap", "found", "read", "cached",
                 "arm", "fingerprint", "justification", "error"]


def run_bench(args, arm, worker, reader=None):
    """Every row of a TSV, in parallel, one process."""
    from .flow import arms
    from .ot_shortlist import encoder, resolve_encoder
    rows = read_bench(args.bench)
    # Traits live under --data-dir; literature runs sit beside them so a
    # benchmark folder holds both and does not mix with data/runs/.
    flow_dir = os.path.join(args.data_dir, "flow")
    os.makedirs(args.data_dir, exist_ok=True)
    os.makedirs(flow_dir, exist_ok=True)
    out = args.out or os.path.join(os.path.dirname(os.path.abspath(args.bench)),
                                   "pesto.tsv")
    from .cost import load_table as load_cost, path_for as cost_path, write_table as write_cost
    costs_out = cost_path(out)
    extra = [k for k in rows[0] if k not in ("gene", "phenotype")]
    fields = ["gene", "phenotype"] + extra + ANSWER_FIELDS

    def fitted(name):
        return arms.get(name).using(reader=reader, worker=worker)

    family = family_of(args)
    if arm is None:
        allowed = {fitted(n).fingerprint() for n in arms.FAMILIES[family]}
        label = "auto: " + " or ".join(arms.FAMILIES[family])
    else:
        allowed = {arm.fingerprint()}
        label = f"{arm.name} [{arm.fingerprint()}]"

    done, lock = {}, threading.Lock()
    costs = load_cost(costs_out) if getattr(args, "resume", False) else {}
    if getattr(args, "resume", False):
        done.update(resumable(out, allowed, fields))
        rows = [r for r in rows if (r["gene"], r["phenotype"]) not in done]
        print(f"resumed {len(done)} rows from {out}", file=sys.stderr)
        if not rows:
            write_table(out, fields, done.values())
            write_cost(costs_out, costs)
            print(f"nothing left to do; wrote {out}", file=sys.stderr)
            return 0
    resumed = len(done)

    # One load, on this thread, before anyone asks for a ranking. After the
    # resume check, so a batch with nothing left to do does not pay for it.
    if resolve_encoder(args.ot_encoder) != "none":
        encoder(args.ot_encoder)

    started = time.time()
    print(f"{len(rows)} pairs, {args.workers} workers, {label}", file=sys.stderr)

    def work(r):
        gene, pheno = r["gene"], r["phenotype"]
        chosen = arm
        if chosen is None:
            with contextlib.redirect_stdout(io.StringIO()):
                chosen = arms.resolve(arms.AUTO, pheno, family=family).using(
                    reader=reader, worker=worker)
        try:
            row = assess(chosen, gene, pheno, args, worker, flow_dir)
        except Exception as exc:
            row = {k: "" for k in fields}
            row.update({"gene": gene, "phenotype": pheno,
                        "error": f"{type(exc).__name__}: {exc}"})
        for k in extra:
            row[k] = r.get(k, "")
        with lock:
            done[(gene, pheno)] = row
            if row.get("_cost"):
                costs[(gene, pheno)] = row["_cost"]
            write_table(out, fields, done.values())
            write_cost(costs_out, costs)
            n = len(done) - resumed
            left = (len(rows) - n) * (time.time() - started) / max(n, 1)
            print(f"  {n:3d}/{len(rows)}  {gene:10s} "
                  f"{row.get('lit_mean') or row.get('error', '')[:20]:14s} "
                  f"{row.get('open_targets_verdict') or '':13s} "
                  f"~{left / 60:.0f} min left", file=sys.stderr, flush=True)
        return row

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(work, rows))
    print(f"wrote {out}", file=sys.stderr)
    print(f"wrote {costs_out}", file=sys.stderr)
    return 0


def family_of(args):
    """Which family `auto` picks a disease or trait arm from."""
    return "fulltext" if getattr(args, "fulltext", False) else "abstracts"


def download_models():
    """Fetch the embedding models now, and report which were already here."""
    from .ot_shortlist import download
    try:
        got = download()
    except Exception as exc:
        print(f"pesto: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    for key, repo, had in got:
        print(f"  {key:<9} {repo:<48} "
              f"{'already cached' if had else 'downloaded'}")
    print("\nWeights live in ~/.cache/huggingface and are shared between "
          "installs.\nSet HF_HUB_OFFLINE=1 to forbid any further download.")
    return 0


def build_parser():
    p = argparse.ArgumentParser(
        prog="pesto",
        description="Is this gene-phenotype association already known?",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Verdicts, weakest to strongest:
              Novel          nothing in the literature links the two
              Hypothesized   suggested, or shown for a closely related phenotype
              Existing       reported, in one study or in a form that falls short
              Established    reported repeatedly, in humans, with genetic evidence

            `pesto browser` opens the same pipeline in a web page.
            """))
    p.add_argument("--gene", metavar="SYMBOL",
                   help="HGNC gene symbol, e.g. DDX41")
    p.add_argument("--phenotype", metavar="NAME",
                   help="disease or trait, spelled as you would search for it")
    p.add_argument("--bench", metavar="TSV",
                   help="table of pairs (columns gene, phenotype). Runs them "
                        "in parallel instead of --gene / --phenotype")
    p.add_argument("--workers", type=int, default=config.DEFAULT_WORKERS,
                   metavar="N",
                   help="parallel pairs when --bench is set "
                        "(default: %(default)s)")
    p.add_argument("--out", metavar="TSV",
                   help="where the table of answers is written (--bench only; "
                        "default: pesto.tsv beside the input)")
    # Not a `choices` list: that would need the arms imported before --help can
    # be printed, and the declaration is the only place their names should live.
    # An unknown name is caught below, by the module that knows them.
    p.add_argument("--arm", default="auto", metavar="NAME",
                   help="which pipeline reads the literature. `auto` asks Opus "
                        "whether the name is a disease and then runs "
                        "`abstracts` or `abstracts-trait` (PESTO as in the "
                        "paper: twenty abstracts weighed in one call). "
                        "`fulltext` and `fulltext-trait` add a Europe PMC "
                        "search of article bodies; `titles`, `knowledge`, "
                        "`abstracts-bare`, `abstracts-open` and "
                        "`abstracts-score` are the paper's comparisons "
                        "(default: %(default)s)")
    p.add_argument("--fulltext", action="store_true",
                   help="with --arm auto, choose between `fulltext` and "
                        "`fulltext-trait` instead")
    p.add_argument("--model", metavar="ID",
                   help="Anthropic model that reads the corpus (default: the "
                        "arm's own, claude-opus-5)")
    p.add_argument("--model-fast", metavar="ID",
                   help="Anthropic model for the cheap steps: expanding the "
                        "phenotype, sieving and banding (default: the arm's "
                        "own, claude-haiku-4-5 and claude-opus-5 for synonyms)")
    p.add_argument("--no-cache", action="store_true",
                   help="read the literature again even if the same pipeline "
                        "has already answered this pair")
    p.add_argument("--resume", action="store_true",
                   help="keep the rows an earlier attempt already answered "
                        "(--bench only) and ask only for the rest. Matches on "
                        "the arm's fingerprint, so it never hands back an "
                        "answer from a different pipeline. Unlike the run "
                        "cache this also skips the Open Targets call, which "
                        "has no cache and is otherwise repaid every time")
    # The Open Targets side. Its default reads the twenty traits nearest the
    # question, sieved by a cheap model and graded one by one by an expensive
    # one, which is about half the price of reading the whole answer and grades
    # the tie rather than only finding it.
    p.add_argument("--ot-encoder", metavar="NAME",
                   help="which model ranks the answer for closeness: `biolord` "
                        "finds ties that run through a disease rather than "
                        "through a word, `sapbert` is a little faster, `none` "
                        "reads the whole answer (default: %s; `none` when "
                        "PyTorch is not installed)" % config.OT_ENCODER)
    p.add_argument("--ot-top", type=int, metavar="N",
                   help="how many of the nearest traits are read; 0 reads them "
                        "all (default: %d)" % config.OT_TOP_K)
    p.add_argument("--no-ot-gate", action="store_true",
                   help="grade every shortlisted trait, without first letting "
                        "a cheap model remove the plainly unrelated ones")
    p.add_argument("--ot-reasons", action="store_true",
                   help="make the grader state why each trait got its tag. "
                        "Costs more and holds the grading up, so it is worth "
                        "it when the grades are being read and not just the "
                        "verdict")
    p.add_argument("--show-pipeline", action="store_true",
                   help="print the arm's blocks and their fingerprints, then stop")
    p.add_argument("--download-models", action="store_true",
                   help="fetch the two embedding models, 439 MB each, and stop. "
                        "Optional: the first assessment that needs one fetches "
                        "it anyway. This only moves that wait somewhere you "
                        "chose, and prepares a machine that will later run "
                        "without a network")
    p.add_argument("--data-dir", metavar="DIR", default=config.NOVEL_RUNS_DIR,
                   help="where a run writes its evidence. A single pair "
                        "defaults to data/runs/novel; --bench defaults to "
                        "runs/ beside the table, so a benchmark does not mix "
                        "with earlier answers")
    p.add_argument("--json", action="store_true",
                   help="print the whole result as JSON instead")
    p.add_argument("--quiet", action="store_true",
                   help="print the two verdicts and nothing else")
    p.add_argument("--verbose", action="store_true",
                   help="report each step on stderr while it runs")
    return p


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    # Read before parsing, because fetching weights is the one thing here that
    # is not about a pair.
    if "--download-models" in argv:
        return download_models()
    if argv[:1] == ["browser"]:
        from .browser.server import main as browser_main
        return browser_main(argv[1:])

    args = build_parser().parse_args(argv)
    if (not config.TORCH_INSTALLED and not args.ot_encoder
            and not os.environ.get("OT_ENCODER") and not args.show_pipeline):
        print("pesto: PyTorch is not installed, so the Open Targets branch grades "
              "every trait associated with the gene instead of the 20 nearest, "
              "which costs up to about $0.50 a pair for well-studied genes. See "
              "'Without PyTorch' in the README.", file=sys.stderr)
    # Describing a pipeline calls nothing, so it needs no key.
    if args.show_pipeline:
        reader, worker = args.model, args.model_fast or args.model
    else:
        reader, worker = pick_models(args.model, args.model_fast or args.model)

    logging.basicConfig(
        stream=sys.stderr, format="%(message)s",
        level=logging.INFO if args.verbose else logging.WARNING)

    from .flow import arms, run as flow_run
    from .flow.types import Query

    def fitted(name):
        try:
            return arms.get(name).using(reader=reader, worker=worker)
        except KeyError as exc:
            raise SystemExit(str(exc).strip('"'))

    auto = args.arm == arms.AUTO
    arm = None if auto else fitted(args.arm)
    if args.show_pipeline:
        if auto:
            disease, trait = arms.FAMILIES[family_of(args)]
            print(fitted(disease).describe())
            print()
            print(fitted(trait).describe())
        else:
            print(arm.describe())
        return 0

    if args.bench:
        if args.gene or args.phenotype:
            raise SystemExit("pesto: --bench cannot be combined with "
                             "--gene / --phenotype")
        # A table of pairs writes its own tree, so it does not mix with
        # whatever a single-pair run left under data/runs/novel.
        if args.data_dir == config.NOVEL_RUNS_DIR:
            args.data_dir = os.path.join(
                os.path.dirname(os.path.abspath(args.bench)), "runs")
        return run_bench(args, arm, worker, reader)

    if not args.gene or not args.phenotype:
        raise SystemExit("pesto: --gene and --phenotype are required, "
                         "or pass --bench")

    query = Query(args.gene, args.phenotype)
    # The services report progress and failures on standard output. That is the
    # channel the answer goes out on, and --json makes it a parsed one, so it is
    # held aside and released to stderr where it belongs.
    chatter = io.StringIO()
    try:
        with contextlib.redirect_stdout(chatter):
            if auto:
                arm = arms.resolve(arms.AUTO, args.phenotype,
                                   family=family_of(args)).using(
                    reader=reader, worker=worker)
            result = flow_run.run(arm, query, cache=not args.no_cache)
            ot_verdict, ot_score, ot_tokens, ot_match = open_targets(
                args.gene, args.phenotype, result.terms.get("synonyms"),
                worker or config.DEFAULT_MODEL, args.data_dir, args)
    except Exception as exc:
        sys.stderr.write(chatter.getvalue())
        if type(exc).__name__ == "AuthenticationError":
            print("pesto: the Anthropic API refused the key in ANTHROPIC_API_KEY "
                  "(401). Check it at console.anthropic.com.", file=sys.stderr)
        else:
            print(f"pesto: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1
    if args.verbose:
        sys.stderr.write(chatter.getvalue())

    from .cost import from_ot_tokens, from_trace, with_total
    cost_rows = with_total(
        from_trace(args.gene, args.phenotype, result.trace, result.cached)
        + from_ot_tokens(args.gene, args.phenotype, ot_tokens,
                         cached=bool(ot_match.get("cached"))),
        args.gene, args.phenotype, result.cached)
    total = next((r for r in cost_rows if r["step"] == "total"), None)

    if args.json:
        print(json.dumps({
            "gene": args.gene, "phenotype": args.phenotype,
            "arm": arm.name, "fingerprint": result.fingerprint,
            "verdict": result.verdict.call,
            "distribution": (result.verdict.distribution.as_dict()
                             if result.verdict.distribution else None),
            "justification": result.verdict.justification,
            "supporting_pmids": list(result.verdict.pmids),
            "open_targets_verdict": ot_verdict,
            "open_targets_max_score": ot_score,
            "counts": result.counts, "cached": result.cached,
            "detail": result.verdict.detail,
            "cost": cost_rows,
        }, indent=2, default=str))
    elif args.quiet:
        print(f"{result.verdict.call}\t{ot_verdict}")
    else:
        print(render(result, ot_verdict, arm, show_reasoning=True))
        if total:
            print(f"\n  cost  ${float(total['usd']):.4f}  "
                  f"({total['input_tokens']} in / {total['output_tokens']} out)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
