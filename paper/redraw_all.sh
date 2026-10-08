#!/usr/bin/env bash
# Redraw every figure and table of the manuscript from paper/results/.
# Calls no model and needs no API key. Outputs go to paper/figures/ and
# paper/results/.
set -euo pipefail
cd "$(dirname "$0")"
export PESTO_PROJECT_ROOT="$PWD"
PY="${PYTHON:-python3}"
mkdir -p figures

echo "== Figure 1, Extended Data Fig. 1"
Rscript scripts/20_figure1.R
Rscript scripts/25_figure_ed1.R
"$PY" scripts/142_fig1_text_numbers.py > results/fig1_text_numbers.txt

echo "== Extended Data Table 1, Supplementary Table 1"
"$PY" scripts/19_build_ed1_v2_table.py
Rscript scripts/22_figure_ed1_v3.R
"$PY" scripts/145_supplementary_table1.py

echo "== Figure 2, Extended Data Fig. 2"
Rscript scripts/21_figure2.R
"$PY" scripts/143_ed2_titles_data.py
Rscript scripts/144_figure_ed2_titles.R

echo "== Figure 3"
Rscript scripts/117_figure3_chd_top5.R --loeuf-decile --brain-decile

echo "== Supplementary Figures 1 and 2"
Rscript scripts/122_figure2_current_vs_score.R
"$PY" scripts/127_fig2_roc_candidate.py
Rscript scripts/128_figure2_roc_candidate.R

echo "== Supplementary Figure 3 (body-only benchmark)"
"$PY" scripts/150_bodyonly_data.py
Rscript scripts/151_figure_supp_bodyonly.R
