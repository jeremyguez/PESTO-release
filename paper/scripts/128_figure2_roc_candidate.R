#!/usr/bin/env Rscript
# Supplementary Figure 2. What the four-colour bars of figure 2 cannot show,
# in two panels.
#
#   a  ROC of the returned label over the four GenCC classes a curator
#      called strong (Definitive, Strong) against the two they did not
#      (Moderate, Limited), one curve per reading prompt, n = 400 pairs.
#      AUCs compared by DeLong on the same pairs.
#   b  How often each prompt calls Established one of the 56 ASC candidate
#      genes, which by the curators' own classification cannot be. Exact
#      binomial intervals, paired mid-p McNemar against abstracts.
#
# Both read the label, not the 100-point distribution, because the label is
# what the pipeline returns. Run scripts/127_fig2_roc_candidate.py first.
#
# Usage: Rscript scripts/128_figure2_roc_candidate.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(readr)
  library(grid)
})

root <- normalizePath(file.path(dirname(sub("--file=", "",
         grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
if (is.na(root) || !dir.exists(root)) root <- normalizePath(".")

RESULTS <- file.path(root, "results")
FIGURES <- file.path(root, "figures")
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)

ARMS <- c("current", "current-bare", "current-open", "score")
# current in black because it is the claim; score in grey because it is the
# floor. The two ablations take the colourblind-safe blue and vermillion.
ACOL <- c(current = "#000000", `current-bare` = "#0072B2",
          `current-open` = "#D55E00", score = "#8C8C8C")
# The tables keep the names the arms had when they were run; the figure shows
# the names the package declares them under.
SHOWN <- c(current = "abstracts", `current-bare` = "abstracts-bare",
           `current-open` = "abstracts-open", score = "abstracts-score")

base_theme <- theme_classic(base_size = 16) +
  theme(
    axis.text  = element_text(colour = "black", size = 14),
    axis.title = element_text(size = 15),
    plot.tag   = element_text(face = "bold", size = 20)
  )

fmt_p <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 1e-4, sprintf("P = %.0e", p), sprintf("P = %.2g", p)))
}

auc <- read_tsv(file.path(RESULTS, "fig2_roc_auc.tsv"), show_col_types = FALSE) %>%
  mutate(arm = factor(arm, levels = ARMS))
roc <- read_tsv(file.path(RESULTS, "fig2_roc_points.tsv"), show_col_types = FALSE) %>%
  mutate(arm = factor(arm, levels = ARMS)) %>%
  arrange(arm, fpr, tpr)
cand <- read_tsv(file.path(RESULTS, "fig2_asd_candidate.tsv"),
                 show_col_types = FALSE) %>%
  mutate(arm = factor(arm, levels = ARMS))

stopifnot(setequal(levels(auc$arm), ARMS), nrow(cand) == 4L)

# ------------------------------------------------------------------ panel a
lab <- auc %>%
  arrange(arm) %>%
  mutate(label = sprintf("%s   AUC %.3f", SHOWN[as.character(arm)], auc))

panel_roc <- ggplot(roc, aes(fpr, tpr, colour = arm)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 0.35, colour = "grey60") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = ACOL, breaks = ARMS, labels = lab$label,
                      guide = guide_legend(ncol = 1)) +
  scale_x_continuous(expand = expansion(mult = 0.01)) +
  scale_y_continuous(expand = expansion(mult = 0.01)) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "false positive rate (Moderate, Limited)",
       y = "true positive rate (Definitive, Strong)",
       title = "Label separates curated strength") +
  base_theme +
  theme(legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "white", colour = "grey80"),
        legend.key.height = unit(0.5, "cm"),
        plot.title = element_text(size = 15, hjust = 0.5))

# ------------------------------------------------------------------ panel b
# One bracket per comparison against current, stacked clear of the tallest
# interval each one spans.
brackets <- cand %>%
  filter(arm != "current") %>%
  transmute(to = as.integer(arm), p = mcnemar_midp) %>%
  arrange(to) %>%
  mutate(from = 1L, y = c(34, 44, 70), tick = 1.8,
         label = fmt_p(p))

panel_cand <- ggplot(cand, aes(arm, rate, fill = arm)) +
  geom_col(width = 0.66) +
  geom_errorbar(aes(ymin = rate_lo, ymax = rate_hi), width = 0.18,
                linewidth = 0.5, colour = "grey25") +
  geom_text(aes(y = -3.2, label = sprintf("%d/%d", established, n)),
            size = 4.2, colour = "grey30", vjust = 1) +
  geom_segment(data = brackets, aes(x = from, xend = to, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey25") +
  geom_segment(data = brackets, aes(x = from, xend = from,
                                    y = y, yend = y - tick),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey25") +
  geom_segment(data = brackets, aes(x = to, xend = to,
                                    y = y, yend = y - tick),
               inherit.aes = FALSE, linewidth = 0.4, colour = "grey25") +
  geom_text(data = brackets, aes(x = (from + to) / 2, y = y + 1,
                                 label = label),
            inherit.aes = FALSE, size = 4.1, colour = "grey25", vjust = 0) +
  scale_fill_manual(values = ACOL, guide = "none") +
  scale_x_discrete(labels = SHOWN) +
  scale_y_continuous(limits = c(-8, 80), breaks = seq(0, 70, 10),
                     expand = expansion(mult = 0)) +
  labs(x = NULL, y = "% called Established",
       title = "ASC candidate genes (n = 56)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        aspect.ratio = 1,
        plot.title = element_text(size = 15, hjust = 0.5))

fig <- (panel_roc | panel_cand) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 20),
        plot.tag.position = c(0, 1))

ggsave(file.path(FIGURES, "figure2_roc_candidate.png"), fig,
       width = 12.4, height = 6.2, dpi = 600, bg = "white")
ggsave(file.path(FIGURES, "figure2_roc_candidate.pdf"), fig,
       width = 12.4, height = 6.2, bg = "white")
cat("Wrote", file.path(FIGURES, "figure2_roc_candidate.png"), "and .pdf\n")
