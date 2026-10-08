"""One run as one HTML file, readable by double-clicking it.

The page is the same one `pesto browser` serves, with the run's data written
into it. Opened from disk it has no server behind it, so it reads and does not
run: the pairs, their summaries, the articles and the links all work.
"""
from __future__ import annotations

import json
import os

from . import runs

STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")
SLOT = '<script id="pesto-export" type="application/json"></script>'


def html(path):
    meta = runs.read_meta(path)
    meta.pop("path", None)
    rows = runs.rows(path)
    details = {}
    for r in rows:
        if r.get("answered"):
            details[f"{r['gene']}\t{r['phenotype']}"] = runs.pair(
                path, r["gene"], r["phenotype"])
    data = json.dumps({"meta": meta, "rows": rows, "details": details})
    # A closing tag inside the JSON would end the script element early.
    data = data.replace("</", "<\\/")
    with open(os.path.join(STATIC, "index.html"), encoding="utf-8") as fh:
        page = fh.read()
    if SLOT not in page:
        raise RuntimeError("the page has no slot for exported data")
    return page.replace(SLOT, SLOT.replace("></script>", f">{data}</script>"))
