"""What must not change, checked without calling a model.

Half of this pipeline is deterministic: which articles a rule selects, what text
a prompt ends up being, how a saved model answer is parsed, which traits Open
Targets counts as a match, and the fingerprint of every arm. Those are frozen
here. Corpora and saved runs come from tests/fixtures/, nothing is paid for,
nothing is written outside this directory, and the whole thing runs in a second.

  python3 tests/replay.py            check against the recorded digests
  python3 tests/replay.py --update   record them again, after a deliberate change

The fingerprint hashes a block's parameters and prompts, not the body of the
function in pesto.harness it calls. A change there can alter what an arm reads
while every fingerprint stays put; the select and render facts below are what
catch it.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(HERE, "fixtures")
os.environ.setdefault("PESTO_PROJECT_ROOT", FIXTURES)

from pesto import config, scoring  # noqa: E402
from pesto.config import APP_ROOT  # noqa: E402
from pesto.flow import arms, steps  # noqa: E402
from pesto.flow.render import render  # noqa: E402
from pesto.flow.types import Corpus, FIELDS, Labels  # noqa: E402
from pesto.open_targets import _derive_open_targets_verdict  # noqa: E402
from pesto.parsers.llm_parsers import parse_novelty_assessment_response  # noqa: E402
from pesto.utils.data_loader import data_loader  # noqa: E402
from pesto.utils.helpers import load_prompt  # noqa: E402

GOLDEN = os.path.join(HERE, "replay.json")
RUNS = os.path.join(FIXTURES, "runs")
CORPORA = os.path.join(FIXTURES, "corpora")
ALIAS_GENES = ["DMD", "APC", "DDX41", "ERBB3", "NLRP1", "TTN"]


def digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()[:16]


def corpora():
    """Three corpora: one large enough for batching to matter, one whose
    disease name no author writes, one whose gene symbol is an English word."""
    return {os.path.basename(p)[:-5]: json.load(open(p, encoding="utf-8"))
            for p in sorted(glob.glob(os.path.join(CORPORA, "*.json")))}


def bands_for(corpus):
    """Bands that do not depend on a model, spread over all four, so the
    selection rules can be checked alone."""
    labels = {}
    for i, a in enumerate(corpus):
        b = 4 - i % 4
        labels[a.pmid] = (b, steps.harness.BANDS[b])
    return labels


def facts():
    out = {"arms": {}, "prompts": {}, "render": {}, "select": {}, "judge": {},
           "ot_verdicts": {}, "novelty_parses": {}, "aliases": {}, "constants": {}}

    for name, arm in sorted(arms.ARMS.items()):
        out["arms"][name] = {"fingerprint": arm.fingerprint(),
                             "blocks": [b.fingerprint() for b in arm.blocks()]}
    out["arms"]["+aliases"] = digest(arms.ALIASES)

    for path in sorted(glob.glob(os.path.join(APP_ROOT, "prompts", "*.txt"))):
        name = os.path.basename(path)[:-4]
        out["prompts"][name] = digest(load_prompt(name))
    # The reading prompts as assembled, not as stored: the probabilities
    # template carries the closing rule of the assessment prompt appended.
    for cls in (steps.ProbabilityRead, steps.ProbabilityReadGeneral,
                steps.BareProbabilityRead, steps.OpenProbabilityRead,
                steps.NoveltyScoreRead, steps.KnowledgeRead):
        out["prompts"][f"+{cls.__name__}"] = digest(cls().template())

    for name, saved in corpora().items():
        corpus = Corpus.of(saved["articles"], FIELDS)
        labels = Labels.from_tiers(bands_for(corpus))
        out["render"][name] = {
            "with_abstracts": digest(render(corpus, steps.ProbabilityRead.needs)),
            "titles_only": digest(render(corpus, steps.TitleKick.needs)),
            "n": len(corpus),
        }
        out["select"][name] = {
            "strength_20_20": list(steps.ByStrength(20, 20).run(corpus, labels).pmids),
            "strength_14_14": list(steps.ByStrength(14, 14).run(corpus, labels).pmids),
            "keep_all": len(steps.KeepAll().run(corpus, labels)),
        }

    q = None
    for text in ("<established>70</established><existing>20</existing>"
                 "<hypothesized>10</hypothesized><novel>0</novel>"
                 "<justification>two families</justification>",
                 "<novelty>0</novelty>", "<novelty>50</novelty>",
                 "<novelty>83</novelty>", "<novelty>100</novelty>", "nothing"):
        out["judge"][text[:40]] = [steps.Argmax().run(q, text).call,
                                   steps.ScoreRule().run(q, text).call]

    # Read off saved runs without asking anything of a model: how the Open
    # Targets payload becomes a verdict, and how a recorded answer is parsed.
    for path in sorted(glob.glob(os.path.join(RUNS, "*.json"))):
        d = json.load(open(path, encoding="utf-8"))
        key = os.path.basename(path)
        traits = (d.get("api_response") or {}).get("open_targets_traits")
        if traits:
            out["ot_verdicts"][key] = digest(_derive_open_targets_verdict(traits))
        raw, pmids = d.get("llm_raw_response"), d.get("pubmed_pmids") or []
        if raw and pmids:
            out["novelty_parses"][key] = digest(
                parse_novelty_assessment_response(raw, pmids))

    # Sorted: the loader builds them from a set.
    for gene in ALIAS_GENES:
        out["aliases"][gene] = digest(sorted(data_loader.get_aliases_for_gene(gene)))

    out["constants"] = {
        "config": digest({"verdict_order": config.VERDICT_ORDER,
                          "default_temperature": config.DEFAULT_TEMPERATURE,
                          "ot_match_version": config.OT_MATCH_VERSION}),
        "scoring": digest({"hypothesized_below": scoring.OT_HYPOTHESIZED_BELOW,
                           "existing_upto": scoring.OT_EXISTING_UPTO}),
        "reading": digest({"budget": steps.harness.READING_BUDGET,
                           "quotas": [steps.harness.SPECIFIC_QUOTA,
                                      steps.harness.BROAD_QUOTA,
                                      steps.harness.SECONDARY_QUOTA,
                                      steps.harness.OLDEST_QUOTA,
                                      steps.harness.VARIANT_QUOTA,
                                      steps.harness.DECOMPOSED_QUOTA,
                                      steps.harness.GENETICS_QUOTA,
                                      steps.harness.NAME_HEAD_QUOTA],
                           "genetics_terms": steps.harness.GENETICS_TERMS,
                           "models": steps.harness.MODELS}),
    }
    return out


def compare(was, now):
    """Changed, gone, and new, which are three different pieces of news."""
    moved, gone, added = [], [], []
    for section in sorted(set(was) | set(now)):
        old, new = was.get(section, {}), now.get(section, {})
        for key in sorted(set(old) | set(new)):
            name = f"{section}/{key}"
            if key not in new:
                gone.append(name)
            elif key not in old:
                added.append(name)
            elif old[key] != new[key]:
                moved.append((name, old[key], new[key]))
    return moved, gone, added


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--update", action="store_true")
    args = ap.parse_args()

    was = (json.load(open(GOLDEN, encoding="utf-8"))
           if os.path.exists(GOLDEN) else {})
    now = facts()
    total = sum(len(v) for v in now.values())
    if args.update or not was:
        json.dump(now, open(GOLDEN, "w", encoding="utf-8"), indent=1, sort_keys=True)
        print(f"recorded {total} facts to {os.path.relpath(GOLDEN)}")
        return 0

    moved, gone, added = compare(was, now)
    for name, old, new in moved:
        print(f"MOVED  {name}\n  was {old}\n  now {new}")
    for name in gone:
        print(f"GONE   {name}, which the reference expects")
    if added:
        print(f"new since the reference, not a failure: {len(added)} "
              f"({', '.join(added[:3])}{', ...' if len(added) > 3 else ''})")
    print(f"{total - len(moved) - len(added)}/{total - len(added)} unchanged")
    return 1 if (moved or gone) else 0


if __name__ == "__main__":
    sys.exit(main())
