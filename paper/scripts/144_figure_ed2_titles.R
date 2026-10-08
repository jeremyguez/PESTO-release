#!/usr/bin/env Rscript
# Draft Extended Data Figure 2 - what reading only titles costs.
#
#   a  average precision on the 2023-2025 bar of Fig. 2d (40 Limited vs the
#      100 Absent of Fig. 2a): model internal knowledge, PESTO-titles, PESTO
#   b  USD a pair on the Fig. 2a,b pairs
#
# Reads results/ed2_titles_*.tsv from scripts/143_ed2_titles_data.py.
# Writes figures/figure_ed2_titles.{png,pdf}; figure_ed2.png is left alone.
#
# Usage: Rscript scripts/144_figure_ed2_titles.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(readr)
})

root <- normalizePath(file.path(dirname(sub("--file=", "",
         grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
if (is.na(root) || !dir.exists(root)) root <- normalizePath(".")

RESULTS <- file.path(root, "results")
FIGURES <- file.path(root, "figures")

METHODS <- c("Model knowledge", "PESTO-titles", "PESTO")
MCOL <- c(`Model knowledge` = "#7A4E6A", `PESTO-titles` = "#8A6B4A",
          PESTO = "#2E5F4F")
CHANCE <- 40 / 140

base_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    plot.subtitle = element_text(size = 12, colour = "grey30"),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

fmt_p <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)

ap <- read_tsv(file.path(RESULTS, "ed2_titles_ap.tsv"), show_col_types = FALSE) %>%
  mutate(method = factor(method, levels = METHODS))
tests <- read_tsv(file.path(RESULTS, "ed2_titles_tests.tsv"), show_col_types = FALSE)
cost <- read_tsv(file.path(RESULTS, "ed2_titles_cost.tsv"), show_col_types = FALSE) %>%
  mutate(method = factor(method, levels = METHODS))

p_of <- function(a, b) tests$p_one_sided[tests$better == a & tests$worse == b]
x_of <- function(m) match(m, METHODS)
# The two short brackets share a height, so each stops short of the middle bar.
brackets <- tibble::tibble(
  x1 = c(x_of("PESTO-titles") + 0.06, x_of("Model knowledge"), x_of("Model knowledge")),
  x2 = c(x_of("PESTO"), x_of("PESTO-titles") - 0.06, x_of("PESTO")),
  y = c(1.04, 1.04, 1.13),
  label = c(fmt_p(p_of("PESTO", "PESTO-titles")),
            fmt_p(p_of("PESTO-titles", "Model knowledge")),
            fmt_p(p_of("PESTO", "Model knowledge"))))

pA <- ggplot(ap, aes(x = method, y = ap, fill = method)) +
  geom_col(width = 0.66) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.16, linewidth = 0.5) +
  geom_hline(yintercept = CHANCE, linetype = "dashed", colour = "grey55",
             linewidth = 0.4) +
  geom_segment(data = brackets, aes(x = x1, xend = x2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey30") +
  geom_segment(data = brackets, aes(x = x1, xend = x1, y = y, yend = y - 0.02),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey30") +
  geom_segment(data = brackets, aes(x = x2, xend = x2, y = y, yend = y - 0.02),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey30") +
  geom_text(data = brackets, aes(x = (x1 + x2) / 2, y = y + 0.025, label = label),
            inherit.aes = FALSE, size = 3.6, colour = "grey20") +
  scale_fill_manual(values = MCOL, guide = "none") +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "a",
       x = NULL, y = "Average precision") +
  base_theme

pB <- ggplot(cost, aes(x = method, y = mean_usd, fill = method)) +
  geom_col(width = 0.66) +
  geom_errorbar(aes(ymin = mean_usd - sd_usd, ymax = mean_usd + sd_usd),
                width = 0.16, linewidth = 0.5) +
  geom_text(aes(y = mean_usd + sd_usd + 0.006, label = sprintf("$%.3f", mean_usd)),
            size = 3.8, colour = "grey20") +
  scale_fill_manual(values = MCOL, guide = "none") +
  scale_y_continuous(labels = function(v) sprintf("$%.2f", v),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(title = "b",
       x = NULL, y = "USD per pair (mean ± s.d.)") +
  base_theme

fig <- pA | pB

ggsave(file.path(FIGURES, "figure_ed2_titles.png"), fig,
       width = 9.5, height = 5.2, dpi = 300, bg = "white")
ggsave(file.path(FIGURES, "figure_ed2_titles.pdf"), fig,
       width = 9.5, height = 5.2, bg = "white", device = cairo_pdf)
cat("Wrote", file.path(FIGURES, "figure_ed2_titles.png"), "and .pdf\n")
