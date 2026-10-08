#!/usr/bin/env Rscript
# Supplementary Figure 1. Figure 2a and 2b under abstracts, the 0–100
# novelty-score prompt (abstracts-score), abstracts-bare (category definitions
# only), and abstracts-open (intensity wording only). The result tables keep
# the names these arms had when they were run: current, score, current-bare
# and current-open.
#
#   a abstracts       GenCC    b abstracts       ASD
#   c abstracts-score GenCC    d abstracts-score ASD
#   e abstracts-bare  GenCC    f abstracts-bare  ASD
#   g abstracts-open  GenCC    h abstracts-open  ASD
#
# GenCC columns are n = 100 per class (original 30 plus extra30 and
# extra40). Score panels map novelty to the same four labels as abstracts:
# 1 + 3 * (100 - novelty) / 100, rounded, clipped. Bare and open use the
# same mean_call as abstracts on the four-bin distribution.
#
# Usage: Rscript scripts/122_figure2_current_vs_score.R

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

VERDICTS <- c("Novel", "Hypothesized", "Existing", "Established")
VCOL <- c(Novel = "#A9573F", Hypothesized = "#BFA45C",
          Existing = "#5B7F95", Established = "#2E5F4F")
LEV_A <- c("Definitive", "Strong", "Moderate", "Limited", "Absent")
ASD_CATS <- c("ASD known", "ASD candidate", "novel candidate")
B_CATS <- c(ASD_CATS, "not in ASC")

base_theme <- theme_classic(base_size = 16) +
  theme(
    axis.text     = element_text(colour = "black", size = 14),
    axis.title    = element_text(size = 15),
    legend.title  = element_blank(),
    legend.text   = element_text(size = 14),
    legend.key.size = unit(0.65, "cm"),
    plot.tag      = element_text(face = "bold", size = 20)
  )

as_verdict <- function(x) factor(x, levels = VERDICTS)

mean_call <- function(score) {
  i <- pmax(1L, pmin(4L, as.integer(round(as.numeric(score)))))
  VERDICTS[i]
}

score_call <- function(novelty) {
  mean_call(1 + 3 * (100 - as.numeric(novelty)) / 100)
}

asd_group <- function(classification) {
  dplyr::recode(classification,
    `autism known` = "ASD known",
    `ID known` = "ASD candidate",
    `autism candidate` = "ASD candidate",
    `novel candidate` = "novel candidate",
    .default = NA_character_)
}

verdict_fill <- function() {
  scale_fill_manual(
    values = VCOL, drop = FALSE, breaks = VERDICTS,
    guide = guide_legend(nrow = 1, byrow = TRUE))
}

bar_panel <- function(d, x, xlab, ylab = "% of pairs") {
  d %>% count(.data[[x]], verdict) %>%
    group_by(.data[[x]]) %>% mutate(pct = 100 * n / sum(n), cls_n = sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = .data[[x]], y = pct, fill = verdict)) +
    geom_col(width = 0.72) +
    geom_text(aes(x = .data[[x]], y = 101, label = paste0("n=", cls_n)),
              data = ~ distinct(.x, .data[[x]], cls_n), inherit.aes = FALSE,
              size = 4.2, colour = "grey35", vjust = 0) +
    verdict_fill() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.07))) +
    labs(x = xlab, y = ylab) +
    base_theme +
    theme(legend.position = "bottom",
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 25, hjust = 1))
}

current_2a <- function() {
  lev <- LEV_A
  d <- read_tsv(file.path(root, "benchmark", "gencc_clingen_g2p",
                          "results", "scored.tsv"),
                show_col_types = FALSE) %>%
    filter(lit_mean %in% VERDICTS, gencc_class %in% lev) %>%
    transmute(gencc = factor(gencc_class, levels = lev),
              verdict = as_verdict(lit_mean))
  extra30 <- read_tsv(file.path(RESULTS, "fig2a_extra30", "current_calls.tsv"),
                      show_col_types = FALSE) %>%
    transmute(gencc = factor(gencc_class, levels = lev),
              verdict = as_verdict(lit_mean))
  extra40 <- read_tsv(file.path(RESULTS, "fig2a_extra40", "current_calls.tsv"),
                      show_col_types = FALSE) %>%
    transmute(gencc = factor(gencc_class, levels = lev),
              verdict = as_verdict(lit_mean))
  extra <- bind_rows(lapply(
    c("fig2a_absent_extra20", "fig2a_absent_extra20b",
      "fig2a_absent_extra30"),
    function(dir) {
      read_tsv(file.path(RESULTS, dir, "pesto.tsv"),
               show_col_types = FALSE) %>%
        filter(is.na(error) | error == "") %>%
        transmute(gencc = factor("Absent", levels = lev),
                  verdict = as_verdict(mean_call(lit_score)))
    }))
  d <- bind_rows(d, extra30, extra40, extra)
  stopifnot(sum(d$gencc == "Absent") == 100L)
  stopifnot(all(sapply(c("Definitive", "Strong", "Moderate", "Limited"),
                       function(x) sum(d$gencc == x) == 100L)))
  bar_panel(d, "gencc", "GenCC classification") +
    geom_vline(xintercept = 4.5, linetype = "dashed", linewidth = 0.35,
               colour = "grey55")
}

not_in_asc_current <- function(exclude) {
  parts <- list()
  f101 <- file.path(root, "benchmark", "asd_extra101", "pesto.tsv")
  parts[[1]] <- read_tsv(f101, show_col_types = FALSE) %>%
    filter(set == "random", lit_mean %in% VERDICTS) %>%
    arrange(gene) %>%
    slice_head(n = 50)
  f50 <- file.path(root, "benchmark", "extra50_random", "pesto.tsv")
  parts[[2]] <- read_tsv(f50, show_col_types = FALSE) %>%
    filter(phenotype == "Autism Spectrum Disorder", lit_mean %in% VERDICTS)
  bind_rows(parts) %>%
    mutate(gene = toupper(gene)) %>%
    filter(!gene %in% exclude) %>%
    distinct(gene, .keep_all = TRUE) %>%
    transmute(group = factor("not in ASC", levels = B_CATS),
              verdict = as_verdict(lit_mean))
}

current_2b <- function() {
  curated <- read_tsv(file.path(RESULTS, "bench_asd_current.tsv"),
                      show_col_types = FALSE) %>%
    mutate(gene = toupper(gene),
           group = factor(asd_group(classification), levels = B_CATS)) %>%
    filter(!is.na(group), lit_mean %in% VERDICTS) %>%
    transmute(gene, group, verdict = as_verdict(lit_mean))
  extra <- not_in_asc_current(curated$gene)
  d <- bind_rows(curated %>% select(group, verdict), extra)
  bar_panel(d, "group", NULL, "% of genes")
}

score <- read_tsv(file.path(RESULTS, "arm_current-score_fig2ab.tsv"),
                  show_col_types = FALSE) %>%
  filter(!is.na(novelty), is.na(error) | error == "") %>%
  mutate(verdict = as_verdict(score_call(novelty)))

score_2a <- function() {
  extra_calls <- bind_rows(lapply(
    c("fig2a_extra30", "fig2a_extra40"),
    function(dir) {
      read_tsv(file.path(RESULTS, dir, "current_calls.tsv"),
               show_col_types = FALSE) %>%
        transmute(gencc = factor(gencc_class, levels = LEV_A),
                  verdict = as_verdict(score_mean))
    }))
  d <- bind_rows(
    score %>%
      filter(panel == "2a") %>%
      transmute(gencc = factor(group, levels = LEV_A), verdict),
    extra_calls)
  stopifnot(n_distinct(d$gencc) == 5L)
  stopifnot(all(sapply(c("Definitive", "Strong", "Moderate", "Limited"),
                       function(x) sum(d$gencc == x) == 100L)))
  bar_panel(d, "gencc", "GenCC classification") +
    geom_vline(xintercept = 4.5, linetype = "dashed", linewidth = 0.35,
               colour = "grey55")
}

score_2b <- function() {
  d <- score %>%
    filter(panel == "2b") %>%
    transmute(group = factor(group, levels = B_CATS), verdict)
  bar_panel(d, "group", NULL, "% of genes")
}

bare <- read_tsv(file.path(RESULTS, "arm_current-bare_fig2ab.tsv"),
                 show_col_types = FALSE) %>%
  filter(lit_mean %in% VERDICTS, is.na(error) | error == "") %>%
  mutate(verdict = as_verdict(lit_mean))

bare_2a <- function() {
  d <- bare %>%
    filter(panel == "2a") %>%
    transmute(gencc = factor(group, levels = LEV_A), verdict)
  stopifnot(n_distinct(d$gencc) == 5L)
  stopifnot(all(sapply(c("Definitive", "Strong", "Moderate", "Limited"),
                       function(x) sum(d$gencc == x) == 100L)))
  stopifnot(sum(d$gencc == "Absent") == 100L)
  bar_panel(d, "gencc", "GenCC classification") +
    geom_vline(xintercept = 4.5, linetype = "dashed", linewidth = 0.35,
               colour = "grey55")
}

bare_2b <- function() {
  d <- bare %>%
    filter(panel == "2b") %>%
    transmute(group = factor(group, levels = B_CATS), verdict)
  bar_panel(d, "group", NULL, "% of genes")
}

open <- read_tsv(file.path(RESULTS, "arm_current-open_fig2ab.tsv"),
                 show_col_types = FALSE) %>%
  filter(lit_mean %in% VERDICTS, is.na(error) | error == "") %>%
  mutate(verdict = as_verdict(lit_mean))

open_2a <- function() {
  d <- open %>%
    filter(panel == "2a") %>%
    transmute(gencc = factor(group, levels = LEV_A), verdict)
  stopifnot(n_distinct(d$gencc) == 5L)
  stopifnot(all(sapply(c("Definitive", "Strong", "Moderate", "Limited"),
                       function(x) sum(d$gencc == x) == 100L)))
  stopifnot(sum(d$gencc == "Absent") == 100L)
  bar_panel(d, "gencc", "GenCC classification") +
    geom_vline(xintercept = 4.5, linetype = "dashed", linewidth = 0.35,
               colour = "grey55")
}

open_2b <- function() {
  d <- open %>%
    filter(panel == "2b") %>%
    transmute(group = factor(group, levels = B_CATS), verdict)
  bar_panel(d, "group", NULL, "% of genes")
}

fig <- ((current_2a() + labs(title = "abstracts") |
         current_2b() + labs(title = "abstracts")) /
        (score_2a() + labs(title = "abstracts-score") |
         score_2b() + labs(title = "abstracts-score")) /
        (bare_2a() + labs(title = "abstracts-bare") |
         bare_2b() + labs(title = "abstracts-bare")) /
        (open_2a() + labs(title = "abstracts-open") |
         open_2b() + labs(title = "abstracts-open"))) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        legend.box.margin = margin(t = -6, b = 0),
        plot.title = element_text(size = 15, hjust = 0.5),
        plot.tag = element_text(face = "bold", size = 20),
        plot.tag.position = c(0, 1))

ggsave(file.path(FIGURES, "figure2_current_vs_score.png"), fig,
       width = 11, height = 22.4, dpi = 600, bg = "white")
ggsave(file.path(FIGURES, "figure2_current_vs_score.pdf"), fig,
       width = 11, height = 22.4, bg = "white")
cat("Wrote", file.path(FIGURES, "figure2_current_vs_score.png"), "and .pdf\n")
