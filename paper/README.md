# Reproducing the paper

Everything needed to redraw the figures and tables of *Large language model
classifies prior evidence in gene–phenotype associations*. The scripts compute
the project root as the parent of `scripts/`, so they read and write under this
directory.

## Redraw everything

```bash
pip install pesto-genetics scipy openpyxl    # or `pip install -e ..` from a clone
# R 4.x with: ggplot2 dplyr tidyr readr patchwork cowplot ggforce ggpubr
#             gridExtra gtable jsonlite
bash redraw_all.sh
```

This calls no model and needs no API key. It takes about two minutes.

| Item | Script(s) | Output |
|---|---|---|
| Figure 1 | `20_figure1.R` (+ `pipeline_schema.R`) | `figures/figure1_v7.png` |
| Numbers quoted for Figure 1 | `142_fig1_text_numbers.py` | `results/fig1_text_numbers.txt` |
| Extended Data Fig. 1 | `25_figure_ed1.R` | `figures/figure_ed1.png` |
| Extended Data Table 1 | `19_build_ed1_v2_table.py`, `22_figure_ed1_v3.R` | `figures/table_ed1_v3.png` |
| Supplementary Table 1 | `145_supplementary_table1.py` | `results/supplementary_table1.xlsx` |
| Figure 2 | `21_figure2.R` | `figures/figure2.png` |
| Extended Data Fig. 2 | `143_ed2_titles_data.py`, `144_figure_ed2_titles.R` | `figures/figure_ed2_titles.png` |
| Figure 3 | `117_figure3_chd_top5.R --loeuf-decile --brain-decile` (+ `_panel_dendrogram.R`) | `figures/figure3.png` |
| Supplementary Figure 1 | `122_figure2_current_vs_score.R` | `figures/figure2_current_vs_score.png` |
| Supplementary Figure 2 | `127_fig2_roc_candidate.py`, `128_figure2_roc_candidate.R` | `figures/figure2_roc_candidate.png` |
| Supplementary Figure 3 | `150_bodyonly_data.py`, `151_figure_supp_bodyonly.R` | `figures/figure_supp_bodyonly.png` |

Redrawn here, Figures 2 and 3, Extended Data Fig. 2, Extended Data Table 1 and
both tables are identical to the published ones pixel for pixel or byte for byte.
Figure 1 and Extended Data Fig. 1 differ only in where the jittered points fall,
which was random and is now seeded. The Supplementary Figures differ only in
their labels, which carry the arms' current names.

## Names in the tables

The result tables record the names the arms had when they were run. The package
declares them under new names, and accepts the old ones as aliases:

| In the tables | In the package |
|---|---|
| `current` | `abstracts` (PESTO in the manuscript) |
| `current-general` | `abstracts-trait` |
| `prior` | `knowledge` (model internal knowledge) |
| `titles` | `titles` (PESTO-titles) |
| `current-bare`, `current-open`, `score` | `abstracts-bare`, `abstracts-open`, `abstracts-score` |
| `currentv2-cheap` | `fulltext` |

## How the tables were produced

The other scripts in `scripts/` drew the benchmark pairs and asked the models.
They are kept as the record of how every table in `results/` was made. Running
them again costs money and, because Claude Opus does not accept a temperature,
will not reproduce each answer exactly.

| Step | Scripts |
|---|---|
| One arm over a table of pairs | `run_arm.py` (or `pesto --bench`) |
| Figure 1 cohort and negative controls | `11_build_all_runs_auto.py`, `12_draw_cohort_absent.py`, `13_draw_random_absent.py`, `14_replay_ot.py`, `16_rematch_ot_cohort.py`, `17_replay_ot_prompt.py`, `00_aggregate_runs.py` (rebuilds `results/all_runs.tsv` from earlier runs that are not included; the table itself is) |
| Figure 2 benchmarks | `draw_figure2_benchmark.py`, `score_figure2_uncontested.py`, `16_consolidate_benchmark.py`, `123_fig2a_extra30.py`, `85_`/`94_`/`100_draw_fig2a_absent*.py`, `86_`/`95_`/`101_score_fig2a_absent*.py`, `80_prior_year_sample.py`, `102_prior_year10_extra10bin.py`, `103_prior_year10_fill3year.py`, `40_ablation_fig2.py` |
| Extended Data Fig. 2 | `104_score_fig2a_absent70_titles.py`, `78_fig2_ed5_data.py` |
| Figure 3 gene universe and enrichments | `82_fig3_four_phenos.py`, `92_fig3_random50.py`, `104_fig3_alzheimer.py`, `106_fig3_random50b.py`, `115_fig3_chd.py`, `108_fig3_go_nocc.py`, `110_fig3_nocc_meancall.py`, `113_fig3_top5_parent.py`, `116_fig3_chd_top5.py`, `86_go_matrix_5pheno.py`, `42_fig3_cd_bars.py` |
| Supplementary Figures 1 and 2 | `121_score_fig2ab.py`, `125_bare_fig2ab.py`, `126_open_fig2ab.py` |
| Supplementary Figure 3 benchmark | `139_build_bench_bodyonly.py` (first 22 pairs), `153_extend_bench_bodyonly.py` (18 more, same criteria), `138_score_bodyonly.py` |

Their inputs and the saved model answers (one JSON per pair, with the corpus,
the prompts and the response) are on Zenodo,
[doi:10.5281/zenodo.23230222](https://doi.org/10.5281/zenodo.23230222):

```bash
bash fetch_data.sh        # about 105 MB, extracted into paper/
```

Not redistributed: the PEPPER gene features used to draw the random genes of
Figure 3 (`92_fig3_random50.py` says where to put them), the GWAS Catalog used to
build the body-only benchmark, and the embedding weights, which download from
HuggingFace on first use.
