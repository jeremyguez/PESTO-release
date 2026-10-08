"""Classify one EFO term against another using the is_a graph.

This is what replaces substring containment in the Open Targets matcher: once
the query phenotype is resolved to a term, every gene-associated disease that
Open Targets returns already carries an identifier in the same space, so the
relationship between the two is a graph question rather than a string question.
"""

import json
import os
from collections import deque

HERE = os.path.dirname(os.path.abspath(__file__))
TERMS_FILE = os.path.join(HERE, "data", "efo_terms.json")

SAME = "same"
DESCENDANT = "descendant"   # candidate is more specific than the query
ANCESTOR = "ancestor"       # candidate is broader than the query
SIBLING = "sibling"         # shares a parent but is a different trait
OTHER = "unrelated"
UNKNOWN = "not_in_efo"

# Only "same" counts as evidence that the gene has been linked to this exact
# trait. Descendants are informative but narrower; ancestors and siblings are
# the two failure modes the old substring rule silently accepted.
EVIDENCE = {SAME: 1.0, DESCENDANT: 0.6, ANCESTOR: 0.0, SIBLING: 0.0,
            OTHER: 0.0, UNKNOWN: 0.0}


class Graph:
    def __init__(self):
        with open(TERMS_FILE) as fh:
            self.terms = json.load(fh)["terms"]

    def _walk(self, tid, key, cap=100000):
        seen, q = set(), deque([tid])
        while q:
            cur = q.popleft()
            for nxt in self.terms.get(cur, {}).get(key, []):
                if nxt not in seen:
                    seen.add(nxt)
                    q.append(nxt)
                    if len(seen) >= cap:
                        return seen
        return seen

    def ancestors(self, tid):
        return self._walk(tid, "parents")

    def descendants(self, tid):
        return self._walk(tid, "children")

    def canonical(self, tid, _depth=0):
        """Follow retirement links to the term that is current.

        Open Targets is built on an older ontology release than the one pinned
        here, so it still returns identifiers that have since been retired.
        Without this step a retired identifier looks like an unknown trait and
        the association is wrongly reported as unseen.
        """
        t = self.terms.get(tid)
        if not t or _depth > 5:
            return tid
        nxt = t.get("replaced_by")
        if nxt and nxt in self.terms and nxt != tid:
            return self.canonical(nxt, _depth + 1)
        return tid

    def relate(self, query_id, cand_id):
        query_id, cand_id = self.canonical(query_id), self.canonical(cand_id)
        if cand_id not in self.terms or query_id not in self.terms:
            return UNKNOWN
        if cand_id == query_id:
            return SAME
        if cand_id in self.descendants(query_id):
            return DESCENDANT
        if cand_id in self.ancestors(query_id):
            return ANCESTOR
        qp = set(self.terms[query_id]["parents"])
        if qp & set(self.terms[cand_id]["parents"]):
            return SIBLING
        return OTHER

    def label(self, tid):
        t = self.terms.get(tid)
        return t["label"] if t else "(not in EFO slim)"
