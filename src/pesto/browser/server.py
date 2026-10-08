"""`pesto browser`: the pipeline behind a web page, on this computer only.

The server listens on 127.0.0.1 and nowhere else. Every API call must carry the
token printed in the address it opens, and a request whose Host header is not
local is refused, so another web page open in the same browser cannot start a
paid run.

The Anthropic key is read from the environment, as on the command line. When it
is not there the page asks for it; the key then lives in this process's memory
until the server stops. It is never written to disk and never sent anywhere but
api.anthropic.com (see pesto/credentials.py). The optional NCBI key can be
pasted in the settings and is kept the same way, sent only to NCBI.

  pesto browser                     serve and open the page
  pesto browser --runs-dir DIR      keep new runs in DIR (remembered)
  pesto browser --port 8765 --no-open
"""
from __future__ import annotations

import argparse
import json
import os
import secrets
import sys
import threading
import urllib.parse
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import requests

from .. import config
from ..credentials import anthropic_api_key
from ..services import pubmed_service
from . import export, runner, runs

STATIC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")
TOKEN = secrets.token_urlsafe(18)
_KEY_FROM_PAGE = {"set": False}
_KEY_LOCK = threading.Lock()


def _check_key(key):
    """True if Anthropic accepts the key. Lists models: costs nothing."""
    try:
        r = requests.get("https://api.anthropic.com/v1/models",
                         headers={"x-api-key": key,
                                  "anthropic-version": "2023-06-01"},
                         timeout=20)
        return r.status_code == 200
    except requests.RequestException:
        return False


def _check_ncbi_key(key):
    """True if NCBI accepts the key. An empty PubMed search: costs nothing."""
    try:
        r = requests.get(f"{config.BASE_URL_NCBI}esearch.fcgi",
                         params={"db": "pubmed", "term": "pesto", "retmax": 0,
                                 "retmode": "json", "api_key": key},
                         timeout=20)
        return r.status_code == 200 and "error" not in r.json()
    except (requests.RequestException, ValueError):
        return False


def _use_key(key):
    from ..services import llm_service
    with _KEY_LOCK:
        os.environ["ANTHROPIC_API_KEY"] = key
        llm_service._client = None
        _KEY_FROM_PAGE["set"] = True


class Handler(BaseHTTPRequestHandler):
    server_version = "pesto"

    def log_message(self, fmt, *args):
        pass

    # ---------------------------------------------------------------- helpers
    def _local(self):
        host = (self.headers.get("Host") or "").rsplit(":", 1)[0]
        return host in ("127.0.0.1", "localhost")

    def _authorised(self, query):
        return (self.headers.get("X-Pesto-Token") == TOKEN
                or (query.get("token") or [""])[0] == TOKEN)

    def _send(self, code, body, ctype="application/json", extra=None):
        data = body if isinstance(body, bytes) else (
            json.dumps(body).encode() if ctype == "application/json" else body.encode())
        self.send_response(code)
        self.send_header("Content-Type", f"{ctype}; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(data)

    def _body(self):
        n = int(self.headers.get("Content-Length") or 0)
        if n > 50_000_000:
            raise ValueError("request too large")
        return json.loads(self.rfile.read(n) or b"{}")

    def _run_path(self, value):
        path = os.path.abspath(os.path.expanduser(value or ""))
        if not value or not runs.is_run(path) or not runs.allowed(path):
            raise LookupError("no such run")
        return path

    # ------------------------------------------------------------------ verbs
    def do_GET(self):
        url = urllib.parse.urlparse(self.path)
        q = urllib.parse.parse_qs(url.query)
        if not self._local():
            return self._send(403, {"error": "local requests only"})
        if url.path in ("/", "/index.html"):
            if not self._authorised(q):
                return self._send(403, "Open PESTO from the address printed by "
                                       "`pesto browser`.", "text/plain")
            with open(os.path.join(STATIC, "index.html"), encoding="utf-8") as fh:
                return self._send(200, fh.read(), "text/html")
        if not url.path.startswith("/api/") or not self._authorised(q):
            return self._send(403, {"error": "forbidden"})
        try:
            if url.path == "/api/status":
                return self._send(200, {
                    "version": runs._version(),
                    "key": bool(anthropic_api_key()),
                    "key_from_page": _KEY_FROM_PAGE["set"],
                    "ncbi_key": bool(pubmed_service.NCBI_API_KEY),
                    "runs_dir": runs.runs_dir(),
                    "torch": config.TORCH_INSTALLED})
            if url.path == "/api/runs":
                listed = runs.list_runs()
                for m in listed:
                    m["active"] = m["path"] in runner.ACTIVE
                return self._send(200, {"runs": listed})
            if url.path == "/api/run":
                path = self._run_path(q.get("path", [""])[0])
                meta = runs.read_meta(path)
                meta["active"] = path in runner.ACTIVE
                return self._send(200, {"meta": meta, "rows": runs.rows(path)})
            if url.path == "/api/pair":
                path = self._run_path(q.get("path", [""])[0])
                return self._send(200, runs.pair(path, q.get("gene", [""])[0],
                                                 q.get("phenotype", [""])[0]))
            if url.path == "/api/export":
                path = self._run_path(q.get("path", [""])[0])
                name = os.path.basename(os.path.normpath(path)) + ".html"
                return self._send(200, export.html(path), "text/html", {
                    "Content-Disposition": f'attachment; filename="{name}"'})
        except LookupError as exc:
            return self._send(404, {"error": str(exc)})
        except Exception as exc:
            return self._send(500, {"error": f"{type(exc).__name__}: {exc}"})
        return self._send(404, {"error": "not found"})

    def do_POST(self):
        url = urllib.parse.urlparse(self.path)
        if not self._local() or self.headers.get("X-Pesto-Token") != TOKEN:
            return self._send(403, {"error": "forbidden"})
        try:
            body = self._body()
            if url.path == "/api/key":
                key = (body.get("key") or "").strip()
                if not key or not _check_key(key):
                    return self._send(400, {"error": "Anthropic did not accept "
                                                     "this key."})
                _use_key(key)
                return self._send(200, {"ok": True})
            if url.path == "/api/ncbi-key":
                key = (body.get("key") or "").strip()
                if key and not _check_ncbi_key(key):
                    return self._send(400, {"error": "NCBI did not accept "
                                                     "this key."})
                pubmed_service.use_key(key)
                return self._send(200, {"ok": True})
            if url.path == "/api/settings":
                return self._send(200, {"runs_dir": runs.set_runs_dir(
                    body.get("runs_dir") or runs.DEFAULT_RUNS_DIR)})
            if url.path == "/api/runs":
                if not anthropic_api_key():
                    return self._send(400, {"error": "no API key"})
                pairs = [{k: str(v).strip() for k, v in p.items()}
                         for p in body.get("pairs") or []]
                pairs = [p for p in pairs if p.get("gene") and p.get("phenotype")]
                if not pairs:
                    return self._send(400, {"error": "no pairs"})
                mode = "fulltext" if body.get("mode") == "fulltext" else "abstracts"
                path = runs.create(body.get("name") or pairs[0]["gene"], pairs, mode)
                runner.start(path, int(body.get("workers") or 8))
                return self._send(200, {"path": path})
            if url.path == "/api/run/stop":
                runner.stop(self._run_path(body.get("path")))
                return self._send(200, {"ok": True})
            if url.path == "/api/run/resume":
                if not anthropic_api_key():
                    return self._send(400, {"error": "no API key"})
                path = self._run_path(body.get("path"))
                runner.start(path, int(body.get("workers") or 8))
                return self._send(200, {"ok": True})
            if url.path == "/api/open":
                path = os.path.abspath(os.path.expanduser(body.get("path") or ""))
                if not runs.is_run(path):
                    return self._send(400, {"error": "This folder holds no PESTO "
                                                     "run (no pesto.tsv)."})
                runs.remember(path)
                return self._send(200, {"path": path})
            if url.path == "/api/forget":
                runs.forget(os.path.abspath(body.get("path") or ""))
                return self._send(200, {"ok": True})
        except LookupError as exc:
            return self._send(404, {"error": str(exc)})
        except Exception as exc:
            return self._send(500, {"error": f"{type(exc).__name__}: {exc}"})
        return self._send(404, {"error": "not found"})


def _server(port):
    try:
        return ThreadingHTTPServer(("127.0.0.1", port), Handler)
    except OSError:
        return ThreadingHTTPServer(("127.0.0.1", 0), Handler)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="pesto browser", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--runs-dir", metavar="DIR",
                    help="where new runs are kept (remembered for next time; "
                         f"default {runs.DEFAULT_RUNS_DIR})")
    ap.add_argument("--no-open", action="store_true",
                    help="print the address instead of opening a browser")
    args = ap.parse_args(argv)
    if args.runs_dir:
        runs.set_runs_dir(args.runs_dir)
    os.makedirs(runs.runs_dir(), exist_ok=True)
    httpd = _server(args.port)
    url = f"http://127.0.0.1:{httpd.server_address[1]}/?token={TOKEN}"
    print(f"PESTO is running at\n  {url}\nRuns are kept in {runs.runs_dir()}\n"
          f"Press Ctrl+C to stop.", flush=True)
    if not args.no_open:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped; unfinished runs resume from the page.", file=sys.stderr)
    finally:
        httpd.server_close()
    return 0
