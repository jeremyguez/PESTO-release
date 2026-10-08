#!/usr/bin/env Rscript
# Supplementary Figure 3. abstracts and fulltext on the body-only benchmark.
#
#   a  Returned label for the GWAS associations whose gene is named only in
#      the body of the article and for the negative controls, per arm.
#   b  Key PMIDs (the GWAS articles naming the gene) cited in the answer.
#   c  Literature-branch cost per pair, mean and s.d.
#
# Run scripts/150_bodyonly_data.py first.
#
# Usage: Rscript scripts/151_figure_supp_bodyonly.R

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
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)

VERDICTS <- c("Novel", "Hypothesized", "Existing", "Established")
VCOL <- c(Novel = "#A9573F", Hypothesized = "#BFA45C",
          Existing = "#5B7F95", Established = "#2E5F4F")
ARMS <- c("abstracts", "fulltext")
ACOL <- c(abstracts = "#000000", fulltext = "#0072B2")
GROUPS <- c("GWAS, gene in body only", "negative control")

base_theme <- theme_classic(base_size = 16) +
  theme(axis.text = element_text(colour = "black", size = 14),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5))

summ <- read_tsv(file.path(RESULTS, "bodyonly_summary.tsv"),
                 show_col_types = FALSE) %>%
  mutate(arm = factor(arm, levels = ARMS))
verd <- read_tsv(file.path(RESULTS, "bodyonly_verdicts.tsv"),
                 show_col_types = FALSE) %>%
  mutate(arm = factor(arm, levels = ARMS),
         group = factor(group, levels = GROUPS,
                        labels = c(sprintf("GWAS association,\ngene in body only (n = %d)",
                                           summ$positives[1]),
                                   sprintf("negative\ncontrols (n = %d)",
                                           summ$negatives[1]))),
         verdict = factor(verdict, levels = VERDICTS))

panel_a <- ggplot(verd, aes(arm, n, fill = verdict)) +
  geom_col(width = 0.7, position = "fill") +
  facet_wrap(~ group) +
  scale_fill_manual(values = VCOL, drop = FALSE, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(100 * x),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "% of pairs", title = "Returned label") +
  base_theme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 13),
        legend.position = "bottom", legend.text = element_text(size = 13))

panel_b <- ggplot(summ, aes(arm, 100 * pmids_cited / pmids_total, fill = arm)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d/%d", pmids_cited, pmids_total)),
            vjust = -0.4, size = 4.6, colour = "grey25") +
  scale_fill_manual(values = ACOL, guide = "none") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.04))) +
  labs(x = NULL, y = "% of key PMIDs cited", title = "Deciding articles cited") +
  base_theme

panel_c <- ggplot(summ, aes(arm, usd_mean, fill = arm)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = usd_mean - usd_sd, ymax = usd_mean + usd_sd),
                width = 0.18, linewidth = 0.5, colour = "grey25") +
  geom_text(aes(y = usd_mean + usd_sd, label = sprintf("$%.3f", usd_mean)),
            vjust = -0.6, size = 4.6, colour = "grey25") +
  scale_fill_manual(values = ACOL, guide = "none") +
  scale_y_continuous(labels = function(x) sprintf("$%.2f", x),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "USD per pair (mean ± s.d.)",
       title = "Literature branch cost") +
  base_theme

fig <- (panel_a | panel_b | panel_c) +
  plot_layout(widths = c(2, 1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 20))

out <- file.path(FIGURES, "figure_supp_bodyonly")
ggsave(paste0(out, ".png"), fig, width = 14, height = 5.6, dpi = 600, bg = "white")
ggsave(paste0(out, ".pdf"), fig, width = 14, height = 5.6, bg = "white")
cat("Wrote", paste0(out, ".png"), "and .pdf\n")
