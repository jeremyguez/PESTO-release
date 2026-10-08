#!/usr/bin/env python3
"""Check that the paper's scripts can import the package on disk.

Running the scripts to find out would cost money and hours, so this runs their
import statements alone, each in a fresh interpreter, and reports the ones that
no longer resolve. It says nothing about whether a script still produces the
right answer; replay.py next to it does that.

It then checks build/ for a copy of the package left behind by an earlier
install. Both directories are ignored by git, so nothing else in the repository
can see them go stale, and setuptools reuses build/lib rather than rebuilding
it: a tree from before flow/ existed will happily ship a package with no arms
in it. That copy is also why a first count of dead code here came back at zero,
every function having been found twice.
"""
from __future__ import annotations

import ast
import glob
import os
import re
import subprocess
import sys
import textwrap

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Errors that mean the import block is broken, as opposed to the ones that mean
# the fragment we extracted is not a runnable program on its own.
BROKEN = ("ModuleNotFoundError", "ImportError", "SyntaxError")


def preamble(path):
    """The script's imports and the few assignments they depend on."""
    text = open(path, encoding="utf-8").read()
    lines = text.splitlines()
    future, body = [], []
    for node in ast.parse(text, path).body:
        chunk = "\n".join(lines[node.lineno - 1:node.end_lineno])
        if isinstance(node, ast.ImportFrom) and node.module == "__future__":
            future.append(chunk)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            body.append(chunk)
        elif isinstance(node, (ast.Assign, ast.Expr)):
            # ROOT and the project-root pin are set between the imports and
            # come before them in effect, so they have to come along.
            if "PESTO_PROJECT_ROOT" in chunk or "sys.path" in chunk or \
                    chunk.startswith(("ROOT", "_HERE", "HERE")):
                body.append(chunk)
    return "\n".join(future + ["import os, sys", f"__file__ = {path!r}"] +
                     [textwrap.dedent(c) for c in body])


def leftovers():
    """Copies of the package that an install would prefer to the source."""
    said = []
    src, built = os.path.join(ROOT, "src", "pesto"), os.path.join(ROOT, "build", "lib", "pesto")
    if os.path.isdir(built):
        def modules(root):
            return {os.path.relpath(p, root)
                    for p in glob.glob(os.path.join(root, "**", "*.py"), recursive=True)}
        missing = modules(src) - modules(built)
        stale = [m for m in modules(src) & modules(built)
                 if os.path.getmtime(os.path.join(src, m))
                 > os.path.getmtime(os.path.join(built, m))]
        if missing or stale:
            said.append(f"build/lib/pesto is out of date: {len(missing)} modules "
                        f"absent from it, {len(stale)} older than the source"
                        + (f" (missing {', '.join(sorted(missing)[:3])}"
                           f"{', ...' if len(missing) > 3 else ''})" if missing else ""))
    version = re.search(r'^version = "(.+)"', open(os.path.join(ROOT, "pyproject.toml"),
                                                   encoding="utf-8").read(), re.M)
    for wheel in glob.glob(os.path.join(ROOT, "dist", "*.whl")):
        if version and f"-{version.group(1)}-" not in os.path.basename(wheel):
            said.append(f"dist/{os.path.basename(wheel)} is not the declared "
                        f"version {version.group(1)}")
    return said


def main():
    paper = os.path.join(ROOT, "paper")
    scripts = sorted(glob.glob(os.path.join(paper, "scripts", "*.py")))
    bad = []
    for path in scripts:
        r = subprocess.run([sys.executable, "-c", preamble(path)],
                           capture_output=True, text=True, cwd=paper)
        if r.returncode:
            last = (r.stderr.strip().splitlines() or ["?"])[-1]
            if last.split(":")[0] in BROKEN:
                bad.append((os.path.relpath(path, ROOT), last))
    print(f"{len(scripts)} scripts, {len(bad)} whose imports no longer resolve")
    for name, err in bad:
        print(f"  {name}: {err}")
    left = leftovers()
    for line in left:
        print(f"  {line}; remove build/ and dist/ before installing")
    return 1 if (bad or left) else 0


if __name__ == "__main__":
    raise SystemExit(main())
