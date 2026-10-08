#!/usr/bin/env python3
"""Figure 2a's benchmark: gene-disease pairs that no submitter contests.

The panel's first version pooled every GenCC submitter, kept each pair at the
strongest class anyone had filed, and then repaired the result: because
submitters work at different granularities,
`scripts/14_gencc_parent_promotion.py` had to read a pair at its MONDO parent's
class whenever the parent carried a stronger one. Nothing is promoted here, and
nothing needs to be, because a pair whose submitters disagree is not drawn at
all.

Eligibility is one sentence: **every submission for the pair carries the same
class, and that class is on the ordered scale.** Unanimity among several
platforms and a lone platform nobody has echoed both satisfy it -- there is
nothing to reconcile in either case. Orphanet's `Supportive`, which is its only
positive label and sits outside the scale, therefore excludes a pair as soon as
Orphanet has filed it alongside anyone. A curated pair also needs at least one
PMID from the platforms that vote: a Limited with no paper is a leftover label,
not a question the literature can answer.

Pairs are keyed on `(gene, MONDO term)` rather than on the disease title, since
two submitters naming one term differently would otherwise look like two pairs
and their disagreement would go unseen.

One further cut, and it is about the question rather than the answer: no *other*
phenotype of the same gene, close enough to be the same clinical question, may
sit at a different class. A gene curated Definitive for one cardiomyopathy and
Limited for a neighbouring one cannot be asked about either, since the pipeline
reading the literature will meet both. "Close enough" is two tests, either
sufficient: BioLORD-2023 cosine at or above the threshold, or a MONDO path
making one term the other's ancestor, descendant or self. BioLORD rather than
SapBERT because it was measured to be the better of the two at relations that
run through a disease rather than through a shared word -- see
`src/pesto/ot_shortlist.py`.

The fifth column is a negative control, not a weaker rung. A gene is drawn at
random from every gene GenCC records and a phenotype at random from every
disease it records, independently of what the four curated columns happen to
hold, and the pair is kept only where GenCC holds no such association, the
phenotype names no other gene, and the phenotype is neither near nor
ontologically related to any disease its new gene is curated for. Drawing from
the whole database rather than recombining the benchmark's own rows is what
makes the column a sample of the questions nobody has answered, rather than a
shuffle of the questions somebody has. Both pools are real GenCC entries, so a
control gene has a literature and a control phenotype has a name the search can
use: the column is an answer, not an absence of one.

The same cosine threshold serves both cuts. At 0.60 it is deliberately stricter
than the 0.80 used by earlier draws: the pool is large enough to pay for it, and
0.80 let two controls through that were not absent at all -- MED27 against
`complex neurodevelopmental disorder` at 0.79, NDUFS1 against `syndromic complex
neurodevelopmental disorder` at 0.61. Neither is caught by any ontology: the
real MONDO, checked through OLS4, asserts none of those links, so the cosine is
the only instrument that sees them.

Nothing here calls a model or the network, so it is free and repeatable. The
sixth guard `scripts/draw_clingen_mini40.py` applies -- expanding each control
phenotype through Opus and rejecting the ones whose synonyms land on a curated
disease -- costs money and is left out until the table below is agreed.

Usage:
  python3 scripts/draw_figure2_benchmark.py --report-only
  python3 scripts/draw_figure2_benchmark.py --n 30 --absent 30
"""
from __future__ import annotations

import argparse
import collections
import csv
import os
import random
import re
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("PESTO_PROJECT_ROOT", ROOT)
sys.path.insert(0, os.path.join(ROOT, "efo"))

GENCC = os.path.join(ROOT, "data", "raw", "gencc-submissions.tsv")
OUT = os.path.join(ROOT, "results", "bench_figure2_uncontested.tsv")
FUNNEL = os.path.join(ROOT, "results", "bench_figure2_uncontested_funnel.tsv")

# The ordered scale the panel plots. Negative classes are excluded by design: an
# axis running from weak to strong support cannot hold a class meaning the
# association was withdrawn. They still count as "a different class" in the
# ambiguity cut, where they are the loudest kind of difference.
CLASSES = ["Definitive", "Strong", "Moderate", "Limited"]


def norm(s):
    return re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip()


ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--n", type=int, default=30, help="pairs per curated class")
ap.add_argument("--absent", type=int, default=30, help="negative controls")
ap.add_argument("--threshold", type=float, default=0.60,
                help="cosine at which two disease labels are one question")
ap.add_argument("--max-descendants", type=int, default=100,
                help="a term with this many descendants is a category, not a "
                     "question; 0 disables the cut")
ap.add_argument("--encoder", default="biolord", help="biolord or sapbert")
ap.add_argument("--seed", type=int, default=2026)
ap.add_argument("--submitters", default="",
                help="comma-separated platforms that vote on class; empty "
                     "means every GenCC submitter. Absence is always tested "
                     "against the whole database.")
ap.add_argument("--report-only", action="store_true",
                help="print the funnel and stop, writing nothing")
ap.add_argument("--out", default=OUT)
ap.add_argument("--funnel", default="",
                help="funnel table; default is the --out path with _funnel")
ap.add_argument("--keep", default="",
                help="existing pairs TSV; pin rows that still pass the filters")
args = ap.parse_args()
if not args.funnel:
    stem, ext = os.path.splitext(args.out)
    args.funnel = f"{stem}_funnel{ext or '.tsv'}"

rng = random.Random(args.seed)

gencc = pd.read_csv(GENCC, sep="\t", dtype=str, low_memory=False)
gencc = gencc.dropna(subset=["gene_symbol", "disease_title", "classification_title",
                             "disease_curie"])

# ------------------------------------------------------------- eligibility

voters = {s.strip() for s in args.submitters.split(",") if s.strip()}

classes_of = collections.defaultdict(set)
voter_classes_of = collections.defaultdict(set)
submitters_of = collections.defaultdict(set)
voter_submitters_of = collections.defaultdict(set)
titles_of = collections.defaultdict(collections.Counter)
filed_as = collections.defaultdict(set)
cited_of = collections.defaultdict(set)
for gene, title, curie, cls, sub, filed, pmids in zip(
        gencc.gene_symbol, gencc.disease_title, gencc.disease_curie,
        gencc.classification_title, gencc.submitter_title,
        gencc.submitted_as_hgnc_symbol, gencc.submitted_as_pmids):
    key = (gene, curie)
    who = str(sub)
    classes_of[key].add(cls)
    submitters_of[key].add(who)
    titles_of[key][str(title)] += 1
    if isinstance(filed, str) and filed.strip():
        filed_as[key].add(filed.strip())
    if not voters or who in voters:
        voter_classes_of[key].add(cls)
        voter_submitters_of[key].add(who)
        if isinstance(pmids, str) and pmids.strip() and pmids.strip().upper() != "NULL":
            cited_of[key].update(re.findall(r"\d{6,9}", pmids))


def title_of(key):
    """The wording most submitters used, ties broken alphabetically."""
    return min(titles_of[key].items(), key=lambda kv: (-kv[1], kv[0]))[0]


# Class is what the chosen platforms say. Everyone else is ignored for
# eligibility, so an Invitae Strong next to a ClinGen Definitive does not
# throw the pair out. Absence, below, still reads the whole file.
source = voter_classes_of if voters else classes_of
agreed = {k: next(iter(v)) for k, v in source.items()
          if len(v) == 1 and next(iter(v)) in CLASSES}
who = ", ".join(sorted(voters)) if voters else "every submitter"
print(f"{len(classes_of)} (gene, MONDO term) pairs in GenCC, every submitter")
print(f"{len(source)} seen by {who}")
print(f"{len(agreed)} carry one class among those platforms, on the scale")

# --------------------------------------------- the gene's other phenotypes

by_gene = collections.defaultdict(list)
for (gene, curie), cs in classes_of.items():
    by_gene[gene].append((curie, cs))

live = {k for k in agreed if len(by_gene[k[0]]) > 1}
texts = sorted({title_of((gene, curie)) for gene, _ in ((k[0], 0) for k in live)
                for curie, _ in by_gene[gene]})
print(f"{len(live)} of them sit on a gene carrying another phenotype, "
      f"{len(texts)} labels to embed")

from pesto import ot_shortlist  # noqa: E402

enc = ot_shortlist.encoder(args.encoder)
if enc is None:
    sys.exit(f"encoder {args.encoder!r} resolves to none; this draw needs one")
vec = enc.encode(texts, max_length=64)
idx = {t: i for i, t in enumerate(texts)}

import relate  # noqa: E402

graph = relate.Graph()
RELATED = (relate.SAME, relate.DESCENDANT, relate.ANCESTOR)


def canon(curie):
    return graph.canonical(str(curie).replace(":", "_"))


def cosine(a, b):
    if a in idx and b in idx:
        return float(vec[idx[a]] @ vec[idx[b]])
    return 0.0


def sibling_clash(key, cls, threshold):
    """Another phenotype of this gene, at another class, but one question."""
    gene, curie = key
    mine, q = title_of(key), canon(curie)
    for other, other_classes in sorted(by_gene[gene]):
        if other == curie or other_classes == {cls}:
            continue
        label = title_of((gene, other))
        cos = cosine(mine, label)
        if cos >= threshold:
            return (f"{'/'.join(sorted(other_classes))} for {label!r} "
                    f"at cos {cos:.2f}")
        rel = graph.relate(q, canon(other))
        if rel in RELATED:
            return (f"{'/'.join(sorted(other_classes))} for {label!r}, "
                    f"{rel} in MONDO")
    return None


SYMBOLS = set(gencc.gene_symbol.dropna().astype(str))
from pesto.config import APP_ROOT as _PESTO  # noqa: E402
ALIASES = os.path.join(_PESTO, "aliases.tsv")
approved_name = {}
alias_of = {}
gene_aliases = collections.defaultdict(set)
with open(ALIASES, encoding="utf-8") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        approved = (row.get("Approved symbol") or "").strip()
        if not approved:
            continue
        approved_name[approved.upper()] = (row.get("Approved name") or "").strip()
        alias_of[approved.upper()] = approved.upper()
        for cell in (row.get("Previous symbols") or "", row.get("Alias symbols") or ""):
            for tok in re.split(r"[,;]", cell):
                tok = tok.strip()
                if tok:
                    alias_of.setdefault(tok.upper(), approved.upper())
                    gene_aliases[approved.upper()].add(tok)

# Definitive is the top of the scale. enumerate(CLASSES) would invert it.
RANK = {c: i for i, c in enumerate(reversed(CLASSES))}
# ClinGen (and sometimes G2P) files a heading when the gene already has a
# named disease at a higher class. Asking the heading then asks the wrong
# question. A heading versus another heading is left alone.
COARSE = {
    "amyotrophic lateral sclerosis",
    "arthrogryposis syndrome",
    "ciliopathy",
    "complex neurodevelopmental disorder",
    "congenital heart disease",
    "congenital nervous system disorder",
    "dilated cardiomyopathy",
    "familial thoracic aortic aneurysm and aortic dissection",
    "hereditary nonpolyposis colon cancer",
    "intellectual disability",
    "intestinal cancer",
    "Leigh syndrome",
    "mitochondrial disease",
    "neurodevelopmental disorder",
    "pulmonary arterial hypertension",
    "retinitis pigmentosa",
}
ENZYME = (
    "dehydrogenase", "synthase", "synthetase", "phosphatase", "bisphosphatase",
    "transferase", "reductase", "oxidase", "kinase", "lyase", "hydratase",
    "carboxylase", "acetyltransferase",
)


def names_a_gene(phenotype, gene=""):
    """Whether the label carries a gene symbol, its own or anyone's.

    `ATF6-related retinopathy` asked about ATF6 is not a question: the search
    finds the gene in the name of the disease. GenCC names diseases this way
    once an association is beyond doubt, so the wording tracks the class --
    18% of Definitive pairs against under 3% everywhere else -- and leaving them
    in would make the strongest column the easiest one for a reason that has
    nothing to do with the strength of the evidence.

    Symbols that are also ordinary words are excluded by the case test, and
    labels like `C3 glomerulonephritis` are lost with them. That is the
    conservative direction, and the pool can afford it. The pair's own HGNC
    name is also checked, so `recombinase activating gene 2 deficiency` is
    seen as RAG2 even though the symbol never appears.
    """
    pheno_n = norm(phenotype)
    if "deficien" in pheno_n and any(e in pheno_n for e in ENZYME):
        return True
    tokens = re.split(r"[^A-Za-z0-9]+", str(phenotype))

    def token_is_a_gene(t):
        if t.islower():
            return False
        if t in SYMBOLS:
            return True
        # Aliases like DK1, not single letters (T) or short English (AS).
        return t in alias_of and (any(c.isdigit() for c in t) or len(t) >= 4)

    if any(token_is_a_gene(t) for t in tokens):
        return True
    if gene:
        own = {gene} | gene_aliases.get(gene.upper(), set())
        if any(t in own and not t.islower() for t in tokens):
            return True
    if not gene:
        return False
    name = approved_name.get(gene.upper(), "")
    if not name:
        return False
    stop = {"gene", "protein", "the", "of", "and", "a"}
    need = [w for w in norm(name).split() if w not in stop and not w.isdigit()]
    have = [w for w in norm(phenotype).split() if w not in stop]
    if not need:
        return False

    def hit(word):
        if word in have:
            return True
        if len(word) < 6:
            return False
        return any(h.startswith(word[:6]) and word.startswith(h[:6])
                   for h in have if len(h) >= 6)

    return all(hit(w) for w in need)


def wrong_gene(key):
    """The row was filed under a different gene, then remapped.

    SLC4A2 / spherocytosis type 4 is the TRPV6 case again: G2P's
    `submitted_as_hgnc_symbol` is SLC4A1, which is the OMIM gene. Renames
    (CCDC103 to DNAAF19) resolve to the same approved symbol and stay.
    """
    gene = alias_of.get(key[0].upper(), key[0].upper())
    for filed in filed_as.get(key, ()):
        other = alias_of.get(filed.upper(), filed.upper())
        if other != gene:
            return filed
    return None


def weaker_bin(key, cls):
    """A coarse label, while the same gene has a named disease at a higher class."""
    title = title_of(key)
    if title not in COARSE:
        return None
    gene, curie = key
    mine = RANK[cls]
    for other, other_classes in by_gene[gene]:
        if other == curie:
            continue
        label = title_of((gene, other))
        if label in COARSE:
            continue
        higher = [c for c in other_classes if c in RANK and RANK[c] > mine]
        if higher:
            return f"{max(higher, key=RANK.get)} for {label!r}"
    return None


def descendants(curie):
    tid = canon(curie)
    return len(graph.descendants(tid)) if tid in graph.terms else 0


def is_category(curie):
    """A term with a large progeny is a heading, not a question.

    MONDO numbers one subtype per gene found, so this counts genetic maturity as
    much as vagueness: Noonan syndrome has ten descendants and is a perfectly
    concrete question. The cut is therefore loose, and it is the cosine that
    does the fine work -- `syndromic complex neurodevelopmental disorder` is an
    umbrella with two descendants and no ontology sees it.
    """
    return bool(args.max_descendants) and descendants(curie) >= args.max_descendants


HEADINGS = {"congenital nervous system disorder", "thrombotic disease"}

clean, ambiguity_of, categories = [], {}, []
wrong, bins, named, uncited = [], [], [], []
for key, cls in agreed.items():
    title = title_of(key)
    why = sibling_clash(key, cls, args.threshold) if key in live else None
    if why:
        ambiguity_of[key] = why
    elif is_category(key[1]) or title in HEADINGS:
        categories.append(key)
    elif names_a_gene(title, key[0]):
        named.append(key)
    elif wrong_gene(key):
        wrong.append(key)
    elif weaker_bin(key, cls):
        bins.append(key)
    elif not cited_of[key]:
        uncited.append(key)
    else:
        clean.append((key, cls))

# ------------------------------------------------------------------ funnel

rows = []
for cls in CLASSES:
    eligible = [k for k, c in agreed.items() if c == cls]
    kept = [k for k, c in clean if c == cls]
    unamb = [k for k in eligible if k not in ambiguity_of]
    specific = [k for k in unamb if not is_category(k[1])]
    unnamed = [k for k in specific if not names_a_gene(title_of(k), k[0])]
    shown = voter_submitters_of if voters else submitters_of
    rows.append({"class": cls, "uncontested": len(eligible),
                 "unambiguous": len(unamb),
                 "specific": len(specific),
                 "unnamed": len(unnamed),
                 "one_platform": sum(1 for k in kept if len(shown[k]) == 1),
                 "several_agree": sum(1 for k in kept if len(shown[k]) > 1),
                 "wanted": args.n})
print(f"\n{'class':12s} {'uncontested':>12s} {'unambiguous':>12s} {'specific':>9s} "
      f"{'unnamed':>8s} {'1 platform':>11s} {'2+ agree':>9s} {'wanted':>7s}")
for r in rows:
    print(f"{r['class']:12s} {r['uncontested']:12d} {r['unambiguous']:12d} "
          f"{r['specific']:9d} {r['unnamed']:8d} {r['one_platform']:11d} "
          f"{r['several_agree']:9d} {r['wanted']:7d}"
          + ("" if r["unnamed"] >= args.n else "   SHORT"))
print(f"\n{len(ambiguity_of)} pairs dropped for a near phenotype at another class:")
for key, why in list(ambiguity_of.items())[:4]:
    print(f"  {key[0]:9s} {title_of(key)[:42]:42s} vs {why}")
print(f"{len(categories)} dropped as headings rather than questions, the widest:")
for key in sorted(categories, key=lambda k: -descendants(k[1]))[:3]:
    print(f"  {descendants(key[1]):4d} descendants  {title_of(key)[:48]}")
print(f"{len(named)} dropped because the phenotype names a gene, "
      f"{len(wrong)} filed under the wrong gene, "
      f"{len(bins)} a coarse bin below a named disease, "
      f"{len(uncited)} with no paper cited by the voting platforms")
short = [r["class"] for r in rows if r["unnamed"] < args.n]
if not args.report_only:
    with open(args.funnel, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
if args.report_only:
    print("\nreport only, nothing written")
    sys.exit(0)
if short:
    sys.exit(f"\n{', '.join(short)} cannot reach {args.n}; decide the fallback")

# -------------------------------------------------------------------- draw

# Scarcest class first, because a gene is spent once and the classes are far
# from equally stocked.
by_class = collections.defaultdict(list)
for key, cls in clean:
    by_class[cls].append(key)
clean_keys = {key: cls for key, cls in clean}

drawn, spent = [], set()
keep_absent = []
if args.keep:
    with open(args.keep, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["gencc_class"] == "Absent":
                keep_absent.append(row)
                continue
            key = (row["gene"], row["disease_curie"])
            cls = row["gencc_class"]
            if clean_keys.get(key) != cls:
                print(f"  drop kept {row['gene']} / {row['phenotype'][:48]} "
                      f"({cls})")
                continue
            if key[0] in spent:
                continue
            spent.add(key[0])
            drawn.append((key, cls))
    print(f"kept {len(drawn)} curated pairs from {args.keep}")

for cls in sorted(CLASSES, key=lambda c: len(by_class[c])):
    have = sum(1 for _, c in drawn if c == cls)
    candidates = sorted(by_class[cls])
    rng.shuffle(candidates)
    taken = have
    for key in candidates:
        if taken == args.n:
            break
        if key[0] in spent:
            continue
        spent.add(key[0])
        drawn.append((key, cls))
        taken += 1
    if taken < args.n:
        sys.exit(f"{cls} ran out at {taken} once scarcer classes had spent genes")
print(f"\ndrew {len(drawn)} curated pairs over {len(spent)} distinct genes")

# ---------------------------------------------------------------- controls

# Absence is tested against every submitter and every class, negative ones
# included: a control must be an association nobody has recorded at all.
curated = collections.defaultdict(set)
for col in ("disease_title", "disease_original_title"):
    for gene, disease in zip(gencc.gene_symbol, gencc[col]):
        if isinstance(disease, str):
            curated[norm(gene)].add(disease)
known = {(g, norm(d)) for g, ds in curated.items() for d in ds}

curies_of_gene = collections.defaultdict(set)
for gene, curie in zip(gencc.gene_symbol, gencc.disease_curie):
    curies_of_gene[norm(gene)].add(canon(curie))

# Both pools are the whole database, not the benchmark's own rows. Phenotypes
# pass the same two label tests the curated columns pass: a control asking
# whether a gene causes `mitochondrial disease` is no more answerable than a
# curated pair saying it does, and one naming a gene answers itself.
gene_pool = sorted(SYMBOLS)
pheno_pool = sorted({title_of(k): k[1] for k in classes_of
                     if not is_category(k[1])
                     and not names_a_gene(title_of(k))}.items())
print(f"\ncontrols drawn from {len(gene_pool)} genes and {len(pheno_pool)} "
      f"phenotypes, the whole of GenCC")

# Cheap tests first, on many more candidates than are needed, so that the model
# is asked to embed once for a batch rather than once per candidate.
absent, seen_g, seen_p = [], set(), set()
for row in keep_absent:
    gene, pheno = row["gene"], row["phenotype"]
    if gene in spent or gene in seen_g or names_a_gene(pheno):
        print(f"  drop kept control {gene} / {pheno[:48]}")
        continue
    absent.append((gene, pheno))
    seen_g.add(gene)
    seen_p.add(pheno)
    spent.add(gene)

candidates = []
while len(candidates) < args.absent * 8 and len(seen_p) < len(pheno_pool):
    gene = rng.choice(gene_pool)
    pheno, curie = rng.choice(pheno_pool)
    if gene in seen_g or pheno in seen_p or gene in spent:
        continue
    if (norm(gene), norm(pheno)) in known or (gene, curie) in classes_of:
        continue
    if any(graph.relate(canon(curie), c) in RELATED
           for c in sorted(curies_of_gene[norm(gene)])):
        continue
    candidates.append((gene, pheno, curie))
    seen_g.add(gene)
    seen_p.add(pheno)

need = [(g, p) for g, p, _ in candidates] + absent
ctexts = sorted({p for _, p in need}
                | {d for g, _ in need for d in curated[norm(g)]})
cvec = enc.encode(ctexts, max_length=64)
cidx = {t: i for i, t in enumerate(ctexts)}


def nearest(gene, pheno):
    sims = [(float(cvec[cidx[pheno]] @ cvec[cidx[d]]), d)
            for d in sorted(curated[norm(gene)]) if d in cidx]
    return max(sims) if sims else (0.0, "")


refused = 0
for gene, pheno, _ in candidates:
    if len(absent) == args.absent:
        break
    if nearest(gene, pheno)[0] >= args.threshold:
        refused += 1
        continue
    absent.append((gene, pheno))
if len(absent) < args.absent:
    sys.exit(f"only {len(absent)} controls survive at cosine {args.threshold}")
print(f"drew {len(absent)} controls, {refused} of {len(candidates)} candidates "
      f"rejected on cosine after the exact and ontology tests")

# ------------------------------------------------------------------- write

fields = ["gene", "phenotype", "gencc_class", "disease_curie", "submitters",
          "n_submitters", "max_cosine", "nearest_curated_disease"]
with open(args.out, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for key, cls in drawn:
        w.writerow({"gene": key[0], "phenotype": title_of(key), "gencc_class": cls,
                    "disease_curie": key[1],
                    "submitters": "; ".join(sorted(
                        (voter_submitters_of if voters else submitters_of)[key])),
                    "n_submitters": len(
                        (voter_submitters_of if voters else submitters_of)[key]),
                    "max_cosine": "", "nearest_curated_disease": ""})
    for gene, pheno in absent:
        cos, d = nearest(gene, pheno)
        w.writerow({"gene": gene, "phenotype": pheno, "gencc_class": "Absent",
                    "disease_curie": "", "submitters": "", "n_submitters": 0,
                    "max_cosine": f"{cos:.3f}", "nearest_curated_disease": d})

repeats = collections.Counter(title_of(k) for k, _ in drawn)
print(f"\n{len(repeats)} distinct phenotypes over {len(drawn)} curated pairs")
for p, n in repeats.most_common(5):
    if n > 1:
        print(f"  {n:3d} x  {p[:58]}")
shown = voter_submitters_of if voters else submitters_of
solo = sum(1 for k, _ in drawn if len(shown[k]) == 1)
print(f"{solo} of {len(drawn)} curated pairs rest on a single platform")
worst = max(absent, key=lambda gp: nearest(*gp)[0])
cos, d = nearest(*worst)
print(f"closest control kept: {worst[0]} / {worst[1]} at cos {cos:.2f} of {d!r}")
print(f"\nwrote {os.path.relpath(args.out, ROOT)} and "
      f"{os.path.relpath(args.funnel, ROOT)}")
