#!/usr/bin/env bash
# Fetch the data archives deposited on Zenodo (doi:10.5281/zenodo.23288837) and
# extract them into paper/. Only the upstream scripts need them: every figure
# and table of the manuscript is redrawn from what the repository already holds.
#
#   bash paper/fetch_data.sh          both archives (about 105 MB)
#   bash paper/fetch_data.sh runs     the saved model answers only
#   bash paper/fetch_data.sh inputs   the upstream inputs only
set -euo pipefail

RECORD="https://zenodo.org/records/23288837/files"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(dirname "$HERE")"

declare -A MD5=(
  [pesto_paper_runs.tar.gz]=ba1a7d997ca942456007efab910b0c82
  [pesto_paper_inputs.tar.gz]=84ebeaa428f766eea6c90039c4901725
)

want=("pesto_paper_runs.tar.gz" "pesto_paper_inputs.tar.gz")
case "${1:-all}" in
  runs)   want=("pesto_paper_runs.tar.gz") ;;
  inputs) want=("pesto_paper_inputs.tar.gz") ;;
  all)    ;;
  *) echo "usage: $0 [runs|inputs|all]" >&2; exit 2 ;;
esac

mkdir -p "$HERE/downloads"
for name in "${want[@]}"; do
  out="$HERE/downloads/$name"
  if [ ! -f "$out" ]; then
    echo "downloading $name"
    curl -L --fail -o "$out" "$RECORD/$name?download=1"
  fi
  got=$(md5sum "$out" | cut -d' ' -f1)
  if [ "$got" != "${MD5[$name]}" ]; then
    echo "$name: checksum $got, expected ${MD5[$name]}" >&2
    exit 1
  fi
  echo "extracting $name into paper/"
  tar -xzf "$out" -C "$REPO"
done
