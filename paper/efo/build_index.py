"""Parse the EFO OTAR slim OBO-graph release into a flat term table.

The release file is the exact ontology Open Targets uses, so the identifier
space here matches the one the Open Targets API returns for associated
diseases. Everything downstream is a pure function of the release pinned in
data/, which is what makes the matcher reproducible.
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SRC = os.path.join(DATA, "efo_otar_slim.json")
OUT = os.path.join(DATA, "efo_terms.json")

SYN_PREDS = {
    "hasExactSynonym": "exact",
    "hasRelatedSynonym": "related",
    "hasNarrowSynonym": "narrow",
    "hasBroadSynonym": "broad",
}


def short_form(iri):
    return iri.rsplit("/", 1)[-1].rsplit("#", 1)[-1]


def build():
    with open(SRC) as fh:
        graph = json.load(fh)["graphs"][0]

    terms = {}
    n_obsolete = 0
    for node in graph["nodes"]:
        if node.get("type") != "CLASS":
            continue
        iri = node.get("id", "")
        if "://" not in iri:
            continue
        tid = short_form(iri)
        label = node.get("lbl")
        meta = node.get("meta") or {}

        # A term can be deprecated either by the explicit flag or, in practice,
        # by an "obsolete_" label prefix. Both must be excluded: an exact match
        # on a dead term looks perfect and is silently wrong.
        deprecated = bool(meta.get("deprecated"))
        if label and label.lower().startswith("obsolete"):
            deprecated = True
        if deprecated:
            n_obsolete += 1

        syns = {"exact": [], "related": [], "narrow": [], "broad": []}
        for s in meta.get("synonyms") or []:
            kind = SYN_PREDS.get(short_form(s.get("pred", "")))
            if kind and s.get("val"):
                syns[kind].append(s["val"])

        gwas = False
        replaced_by = None
        for b in meta.get("basicPropertyValues") or []:
            pred = b.get("pred", "")
            if pred.endswith("gwas_trait") and b.get("val") == "true":
                gwas = True
            # Open Targets lags the ontology by a release or two, so it still
            # returns identifiers that EFO has since retired. Keeping the
            # replacement lets both sides be compared on the live term.
            elif pred.endswith("IAO_0100001") and b.get("val"):
                replaced_by = short_form(b["val"])

        terms[tid] = {
            "id": tid,
            "label": label,
            "definition": (meta.get("definition") or {}).get("val"),
            "synonyms": syns,
            "deprecated": deprecated,
            "replaced_by": replaced_by,
            "gwas_trait": gwas,
            "parents": [],
            "children": [],
        }

    n_edges = 0
    for e in graph.get("edges", []):
        if e.get("pred") not in ("is_a", "subClassOf"):
            continue
        sub, obj = short_form(e["sub"]), short_form(e["obj"])
        if sub in terms and obj in terms:
            terms[sub]["parents"].append(obj)
            terms[obj]["children"].append(sub)
            n_edges += 1

    with open(OUT, "w") as fh:
        json.dump({"version": graph.get("meta", {}).get("version"), "terms": terms}, fh)

    live = sum(1 for t in terms.values() if not t["deprecated"])
    n_syn = sum(len(t["synonyms"]["exact"]) for t in terms.values())
    print(f"terms          {len(terms)}")
    print(f"  live         {live}")
    print(f"  deprecated   {n_obsolete}")
    print(f"is_a edges     {n_edges}")
    print(f"exact synonyms {n_syn}")
    print(f"gwas_trait     {sum(1 for t in terms.values() if t['gwas_trait'])}")
    print(f"replaced_by    {sum(1 for t in terms.values() if t['replaced_by'])}")
    print(f"version        {graph.get('meta', {}).get('version')}")
    print(f"wrote          {OUT}")


if __name__ == "__main__":
    if not os.path.exists(SRC):
        sys.exit(f"missing {SRC}")
    build()
