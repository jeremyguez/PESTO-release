#!/usr/bin/env Rscript
# Figure 2 - validation against two external references, bars then ranking.
#
# Layout:
#   A B
#   C D
#
#   A  ClinGen+G2P literature, one hundred pairs in each curated class (30 of
#      the original ladder plus the extra30 and extra40 draws) and 100 random
#      Absent recombinations (30 of the ladder plus 70 extra)
#   B  consensus autism curation: ASD known, ASD candidate (incl. ID known),
#      novel candidate, plus genes drawn at random and absent from ASC
#   C  mean score of Limited pairs by 3-year first-report bin (2010 aside)
#   D  average precision of those bins vs 100 random Absents (Fig. 2a-style)
#
# Former C and D (40 Limited 2023–2025 vs 40 adjacent Absent) are
# Supplementary Figure 1.
#
# The permutation control moved to Figure 1. The three-model precision-recall
# curves on the weak end of the GenCC ladder stay Extended Data 4.
#
# Every panel is annotated with the model that produced it, because the runs
# were not all made with the same one (see paper/QUESTIONS.md).
#
# Usage: Rscript scripts/21_figure2.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(readr)
  library(jsonlite)
  library(grid)
})

root <- normalizePath(file.path(dirname(sub("--file=", "",
         grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
if (is.na(root) || !dir.exists(root)) root <- normalizePath(".")

RESULTS <- file.path(root, "results")
FIGURES <- file.path(root, "figures")
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)
# opus (default) or haiku: same layout, Haiku pesto + prior. Titles and
# found_* stay the Opus ablation; they are not the comparison.
FIG2_MODEL <- Sys.getenv("PESTO_FIG2_MODEL", "opus")
HAIKU <- file.path(RESULTS, "haiku_fig2")
IS_HAIKU <- identical(FIG2_MODEL, "haiku")

VERDICTS <- c("Novel", "Hypothesized", "Existing", "Established")

# Panel a is the literature reading of the ClinGen+G2P draw
# (benchmark/gencc_clingen_g2p), thirty pairs in each class. The earlier
# all-submitter draw and the opus/claims switch that used to live here sit in
# results/bench_gencc_opus5_absent.tsv and results/bench_kick_vs_opus.tsv.
# Desaturated earth and slate tones. The hues still run warm at the unsupported
# end and cool at the settled one, and their lightness is spread far enough that
# the four levels stay apart in greyscale and under red-green colour blindness,
# where the distance between Novel and Established carries most of the meaning.
VCOL <- c(Novel = "#A9573F", Hypothesized = "#BFA45C",
          Existing = "#5B7F95", Established = "#2E5F4F")

base_theme <- theme_classic(base_size = 16) +
  theme(
    axis.text     = element_text(colour = "black", size = 14),
    axis.title    = element_text(size = 15),
    legend.title  = element_blank(),
    legend.text   = element_text(size = 14),
    legend.key.size = unit(0.65, "cm"),
    plot.tag      = element_text(face = "bold", size = 20)
  )

verdict_fill <- function() {
  scale_fill_manual(
    values = VCOL, drop = FALSE, breaks = VERDICTS,
    guide = guide_legend(nrow = 1, byrow = TRUE))
}

as_verdict <- function(x) factor(x, levels = VERDICTS)

# Same rule as pesto.cli.mean_call and as panels a and b: the 1–4 mean of
# the hundred points, rounded, clipped to the four labels.
mean_call <- function(score) {
  i <- pmax(1L, pmin(4L, as.integer(round(as.numeric(score)))))
  VERDICTS[i]
}

# The permutation control that stood here now lives in scripts/20_figure1.R,
# beside the counts it controls.


# --------------------------------------------------------------------------
# B - agreement with human curators
# --------------------------------------------------------------------------
# Plot labels. The sheet says autism known / autism candidate / ID known;
# ID known is folded into ASD candidate at draw time, the table keeps both.
ASD_CATS <- c("ASD known", "ASD candidate", "Novel candidate")
B_CATS <- c(ASD_CATS, "Not in ASD gene list")
asd_group <- function(classification) {
  dplyr::recode(classification,
    `autism known` = "ASD known",
    `ID known` = "ASD candidate",
    `autism candidate` = "ASD candidate",
    `novel candidate` = "Novel candidate",
    .default = NA_character_)
}

# Random genes scored on ASD and absent from the ASC list: fifty from the
# extra101 random draw (that list has 51 after a replacement; keep the
# first fifty in gene order) plus the later fifty-gene extra50 draw.
# Hardnegatives stay out; they were chosen for transcript count, not at random.
not_in_asc_pesto <- function(exclude) {
  parts <- list()
  f101 <- file.path(root, "benchmark", "asd_extra101", "pesto.tsv")
  if (file.exists(f101)) {
    parts[[length(parts) + 1]] <- read_tsv(f101, show_col_types = FALSE) %>%
      filter(set == "random", lit_mean %in% VERDICTS) %>%
      arrange(gene) %>%
      slice_head(n = 50)
  }
  f50 <- file.path(root, "benchmark", "extra50_random", "pesto.tsv")
  if (file.exists(f50)) {
    parts[[length(parts) + 1]] <- read_tsv(f50, show_col_types = FALSE) %>%
      filter(phenotype == "Autism Spectrum Disorder", lit_mean %in% VERDICTS)
  }
  if (!length(parts)) return(NULL)
  bind_rows(parts) %>%
    mutate(gene = toupper(gene)) %>%
    filter(!gene %in% exclude) %>%
    distinct(gene, .keep_all = TRUE)
}

not_in_asc <- function(exclude) {
  rows <- not_in_asc_pesto(exclude)
  if (is.null(rows) || !nrow(rows)) return(NULL)
  rows %>% transmute(group = factor("Not in ASD gene list", levels = B_CATS),
                     verdict = as_verdict(lit_mean))
}

panel_B <- function() {
  if (IS_HAIKU) {
    f <- file.path(HAIKU, "asd_bars.tsv")
    if (!file.exists(f)) return(patchwork::plot_spacer())
    d <- read_tsv(f, show_col_types = FALSE) %>%
      mutate(group = factor(dplyr::recode(group,
                              `novel candidate` = "Novel candidate",
                              `not in ASC` = "Not in ASD gene list"),
                            levels = B_CATS),
             verdict = as_verdict(lit_mean)) %>%
      filter(!is.na(group), lit_mean %in% VERDICTS) %>%
      select(group, verdict)
  } else {
  f <- file.path(RESULTS, "bench_asd_current.tsv")
  if (!file.exists(f)) return(patchwork::plot_spacer())
  curated <- read_tsv(f, show_col_types = FALSE) %>%
    mutate(gene = toupper(gene),
           group = factor(asd_group(classification), levels = B_CATS)) %>%
    filter(!is.na(group), lit_mean %in% VERDICTS) %>%
    transmute(gene, group, verdict = as_verdict(lit_mean))
  extra <- not_in_asc(curated$gene)
  d <- bind_rows(curated %>% select(group, verdict), extra)
  }

  d %>% count(group, verdict) %>%
    group_by(group) %>% mutate(pct = 100 * n / sum(n), cls_n = sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = group, y = pct, fill = verdict)) +
    geom_col(width = 0.72) +
    geom_text(aes(x = group, y = 101, label = paste0("n=", cls_n)),
              data = ~ distinct(.x, group, cls_n), inherit.aes = FALSE,
              size = 4.2, colour = "grey35", vjust = 0) +
    verdict_fill() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.07))) +
    labs(x = "ASD gene curation", y = "% of genes") +
    base_theme +
    theme(legend.position = "bottom",
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 25, hjust = 1))
}

# --------------------------------------------------------------------------
# C - GenCC positive control
# --------------------------------------------------------------------------
panel_C <- function() {
  # ClinGen and G2P only, no disagreement between them, one hundred pairs in
  # each curated class: the original thirty, plus the extra30 and extra40 draws
  # from the same pool. Absent is 100 random recombinations: the original thirty
  # plus seventy later draws of the same kind (the negatives of panel d).
  lev <- c("Definitive", "Strong", "Moderate", "Limited", "Absent")
  f <- if (IS_HAIKU) file.path(HAIKU, "gencc_scored.tsv") else
    file.path(root, "benchmark", "gencc_clingen_g2p", "results", "scored.tsv")
  if (!file.exists(f)) return(patchwork::plot_spacer())
  d <- read_tsv(f, show_col_types = FALSE) %>%
    filter(lit_mean %in% VERDICTS, gencc_class %in% lev) %>%
    transmute(gencc = factor(gencc_class, levels = lev),
              call = as_verdict(lit_mean))
  if (!IS_HAIKU) {
    # The original 120 keep their verdict from scored.tsv; combined_n100's
    # `current` column reproduces it pair for pair, so only the 280 later draws
    # are read from there.
    curated_extra <- read_tsv(file.path(RESULTS, "fig2a_extra40",
                                        "combined_n100.tsv"),
                              show_col_types = FALSE) %>%
      filter(source != "original", gencc_class %in% lev) %>%
      transmute(gencc = factor(gencc_class, levels = lev),
                call = as_verdict(mean_call(current)))
    extra <- bind_rows(lapply(
      c("fig2a_absent_extra20", "fig2a_absent_extra20b",
        "fig2a_absent_extra30"),
      function(dir) {
        read_tsv(file.path(RESULTS, dir, "pesto.tsv"),
                 show_col_types = FALSE) %>%
          filter(is.na(error) | error == "") %>%
          transmute(gencc = factor("Absent", levels = lev),
                    call = as_verdict(mean_call(lit_score)))
      }))
    d <- bind_rows(d, curated_extra, extra)
    stopifnot(all(table(d$gencc) == 100L))
  }

  d %>% count(gencc, call) %>%
    group_by(gencc) %>% mutate(pct = 100 * n / sum(n), cls_n = sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = gencc, y = pct, fill = call)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 4.5, linetype = "dashed", linewidth = 0.35,
               colour = "grey55") +
    geom_text(aes(x = gencc, y = 101, label = paste0("n=", cls_n)),
              data = ~ distinct(.x, gencc, cls_n), inherit.aes = FALSE,
              size = 4.2, colour = "grey35", vjust = 0) +
    verdict_fill() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.07))) +
    labs(x = "GenCC classification", y = "% of pairs") +
    base_theme +
    theme(legend.position = "bottom",
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 25, hjust = 1))
}

# --------------------------------------------------------------------------
# D - OMIM positive control
# --------------------------------------------------------------------------
MCOL <- c(`claude-haiku-4-5` = "#8C8C8C",
          `gemini-3.1-flash-lite` = "#B07A4B",
          `claude-opus-5` = "#35586B")

# Ranked by score, precision and recall are read at the end of each group of
# equal scores rather than at every pair. Ties are the rule here, not the
# exception: a model that puts all hundred points on one category gives exactly
# 1 or exactly 4, and splitting such a group in an arbitrary order would let the
# average precision depend on the order the rows happen to sit in.
pr_curve <- function(score, positive) {
  o <- order(score, decreasing = TRUE)
  s <- score[o]
  y <- positive[o]
  ends <- which(c(s[-1] != s[-length(s)], TRUE))
  tp <- cumsum(y)[ends]
  precision <- tp / ends
  recall <- tp / sum(y)
  list(curve = tibble::tibble(recall, precision, threshold = s[ends]),
       ap = sum(diff(c(0, recall)) * precision))
}

roc_curve <- function(score, positive) {
  o <- order(score, decreasing = TRUE)
  s <- score[o]
  y <- positive[o]
  ends <- which(c(s[-1] != s[-length(s)], TRUE))
  tp <- cumsum(y)[ends]
  fp <- ends - tp
  tpr <- tp / sum(y)
  fpr <- fp / sum(y == 0)
  fpr <- c(0, fpr)
  tpr <- c(0, tpr)
  list(curve = tibble::tibble(fpr, tpr, threshold = c(Inf, s[ends])),
       auc = sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2))
}

# Panels c and d: the two summary bars of gencc_pr_found_vs_litscore
# (average precision and ROC AUC) in one tagged cell.
SCOL <- c(lit_score = "#2E5F4F",
          titles = "#8A6B4A",
          knowledge = "#7A4E6A",
          found_specific = "#C47A3A",
          found_general = "#5B7F95")
BAR_LAB <- c(lit_score = "PESTO",
             titles = "PESTO-titles",
             knowledge = "Model knowledge",
             found_specific = "n papers (pair)",
             found_general = "n papers (gene)")
BOOT_B <- 2000
PERM_B <- 20000

fmt_p <- function(p) {
  if (p >= 0.001) return(list(label = sprintf("p = %.3f", p), parse = FALSE))
  e <- floor(log10(p))
  m <- p / 10^e
  list(label = sprintf("p == %.2f %%*%% 10^{%d}", m, e), parse = TRUE)
}

lit_score_of <- function(d) {
  (4 * d$p_established + 3 * d$p_existing +
   2 * d$p_hypothesized + 1 * d$p_novel) / 100
}

# First no-article Opus call, and the titles-only arm.
attach_ablations <- function(d) {
  if (IS_HAIKU) {
    know <- read_tsv(file.path(HAIKU, "knowledge_all.tsv"),
                     show_col_types = FALSE) %>%
      transmute(gene = toupper(gene), phenotype,
                knowledge = as.numeric(lit_score))
  } else {
    know <- read_tsv(file.path(RESULTS, "ablation_fig2", "knowledge.tsv"),
                     show_col_types = FALSE) %>%
      mutate(gene = toupper(gene), draw = as.integer(draw),
             lit_score = as.numeric(lit_score)) %>%
      filter(draw == 1) %>%
      transmute(gene, phenotype, knowledge = lit_score)
  }
  titles <- read_tsv(file.path(RESULTS, "ablation_fig2", "titles.tsv"),
                     show_col_types = FALSE) %>%
    filter(error == "" | is.na(error)) %>%
    mutate(gene = toupper(gene), titles = as.numeric(lit_score)) %>%
    select(gene, phenotype, titles)
  d %>%
    mutate(gene = toupper(gene)) %>%
    left_join(know, by = c("gene", "phenotype")) %>%
    left_join(titles, by = c("gene", "phenotype"))
}

haiku_lit <- function() {
  read_tsv(file.path(HAIKU, "pesto_all.tsv"), show_col_types = FALSE) %>%
    transmute(gene = toupper(gene), phenotype,
              lit_score = as.numeric(lit_score))
}

found_auc_panel <- function(d, out_tsv, roc = TRUE, bar_lab = BAR_LAB) {
  scores <- intersect(names(SCOL), names(d))
  stopifnot(length(scores) >= 2, all(c("positive") %in% names(d)))
  aps <- setNames(sapply(scores, function(s) pr_curve(d[[s]], d$positive)$ap),
                  scores)
  aucs <- setNames(sapply(scores, function(s) roc_curve(d[[s]], d$positive)$auc),
                   scores)
  n <- nrow(d)
  boot_ap <- matrix(NA_real_, BOOT_B, length(scores),
                    dimnames = list(NULL, scores))
  boot_auc <- boot_ap
  set.seed(1)
  for (b in seq_len(BOOT_B)) {
    i <- sample.int(n, n, replace = TRUE)
    y <- d$positive[i]
    if (sum(y) == 0L || sum(y) == n) next
    for (s in scores) {
      boot_ap[b, s] <- pr_curve(d[[s]][i], y)$ap
      boot_auc[b, s] <- roc_curve(d[[s]][i], y)$auc
    }
  }
  ci <- function(mat, name) {
    tibble(score = name,
           lo = as.numeric(quantile(mat[, name], 0.025, na.rm = TRUE)),
           hi = as.numeric(quantile(mat[, name], 0.975, na.rm = TRUE)))
  }
  ap_only <- function(score, positive) {
    o <- order(score, decreasing = TRUE)
    s <- score[o]
    y <- positive[o]
    ends <- which(c(s[-1] != s[-length(s)], TRUE))
    tp <- cumsum(y)[ends]
    sum(diff(c(0, tp / sum(y))) * (tp / ends))
  }
  auc_only <- function(score, positive) {
    o <- order(score, decreasing = TRUE)
    s <- score[o]
    y <- positive[o]
    ends <- which(c(s[-1] != s[-length(s)], TRUE))
    tp <- cumsum(y)[ends]
    fp <- ends - tp
    tpr <- c(0, tp / sum(y))
    fpr <- c(0, fp / sum(y == 0))
    sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)
  }
  # One-sided paired permutation: swap the two methods on each pair.
  paired_perm_p <- function(a, b, y, fun) {
    obs <- fun(a, y) - fun(b, y)
    set.seed(1)
    ge <- 0L
    n <- length(y)
    for (i in seq_len(PERM_B)) {
      flip <- runif(n) < 0.5
      sa <- ifelse(flip, b, a)
      sb <- ifelse(flip, a, b)
      if ((fun(sa, y) - fun(sb, y)) >= obs) ge <- ge + 1L
    }
    (1 + ge) / (1 + PERM_B)
  }
  bars <- bind_rows(
    bind_rows(lapply(scores, function(s) {
      ci(boot_ap, s) %>% mutate(value = aps[[s]], metric = "ap")
    })),
    bind_rows(lapply(scores, function(s) {
      ci(boot_auc, s) %>% mutate(value = aucs[[s]], metric = "auc")
    }))
  ) %>% mutate(score = factor(score, levels = scores))
  y <- d$positive
  p_ap_prior <- paired_perm_p(d$lit_score, d$knowledge, y, ap_only)
  p_ap_titles <- paired_perm_p(d$lit_score, d$titles, y, ap_only)
  p_auc_prior <- paired_perm_p(d$lit_score, d$knowledge, y, auc_only)
  p_auc_titles <- paired_perm_p(d$lit_score, d$titles, y, auc_only)
  write_tsv(bars %>%
              mutate(p_pesto_vs_prior = ifelse(metric == "ap",
                                               p_ap_prior, p_auc_prior),
                     p_pesto_vs_titles = ifelse(metric == "ap",
                                                p_ap_titles, p_auc_titles)),
            out_tsv)

  chance <- mean(d$positive)
  x_pesto <- match("lit_score", scores)
  x_prior <- match("knowledge", scores)
  x_titles <- match("titles", scores)
  metric_col <- function(dat, ylab, title, chance_y, p_prior, p_titles, ylim) {
    tick <- 0.012 * diff(ylim)
    label_pad <- 0.055
    step <- 0.075 * diff(ylim)
    y1 <- max(dat$hi[dat$score %in% c("lit_score", "titles")]) + label_pad + step
    y2 <- max(dat$hi[dat$score %in% c("lit_score", "knowledge")]) + label_pad + 2 * step
    bracket <- function(x1, x2, y, pval) {
      list(
        annotate("segment", x = x1, xend = x2, y = y, yend = y,
                 colour = "grey25", linewidth = 0.4),
        annotate("segment", x = x1, xend = x1, y = y - tick, yend = y,
                 colour = "grey25", linewidth = 0.4),
        annotate("segment", x = x2, xend = x2, y = y - tick, yend = y,
                 colour = "grey25", linewidth = 0.4),
        annotate("text", x = (x1 + x2) / 2, y = y,
                 label = fmt_p(pval)$label, parse = fmt_p(pval)$parse,
                 vjust = -0.35, size = 4.0, colour = "grey25")
      )
    }
    ggplot(dat, aes(x = score, y = value, fill = score)) +
      geom_hline(yintercept = chance_y, linetype = "dashed", linewidth = 0.35,
                 colour = "grey60") +
      geom_col(width = 0.62, show.legend = FALSE) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15,
                    linewidth = 0.55, colour = "grey25") +
      geom_text(aes(y = pmin(hi, 1) + 0.018, label = sprintf("%.3f", value)),
                vjust = 0, size = 4.0, colour = "grey25") +
      bracket(x_pesto, x_titles, y1, p_titles) +
      bracket(x_pesto, x_prior, y2, p_prior) +
      scale_fill_manual(values = SCOL) +
      scale_x_discrete(labels = bar_lab) +
      scale_y_continuous(expand = expansion(mult = 0.01), name = ylab) +
      coord_cartesian(ylim = ylim) +
      labs(x = NULL, title = title) +
      base_theme +
      theme(plot.title = element_text(size = 15, hjust = 0.5),
            plot.tag = element_blank(),
            axis.text.x = element_text(angle = 25, hjust = 1))
  }
  ap_ylim <- c(min(0.70, min(filter(bars, metric == "ap")$lo) - 0.04), 1.32)
  auc_ylim <- c(min(0.35, min(filter(bars, metric == "auc")$lo) - 0.04), 1.36)
  ap_plot <- metric_col(filter(bars, metric == "ap"),
                        "Average precision", "Precision-recall",
                        chance, p_ap_prior, p_ap_titles, ap_ylim)
  if (!roc) return(ap_plot)
  inner <- ap_plot |
           metric_col(filter(bars, metric == "auc"),
                      "AUC", "ROC",
                      0.5, p_auc_prior, p_auc_titles, auc_ylim)
  patchwork::wrap_elements(full = inner)
}

gencc_found_data <- function() {
  hits <- read_tsv(file.path(RESULTS, "gencc_found_specific_general.tsv"),
                   show_col_types = FALSE) %>%
    select(gene, phenotype, found_specific, found_general)
  d <- read_tsv(file.path(root, "benchmark", "gencc_clingen_g2p", "pesto.tsv"),
           show_col_types = FALSE) %>%
    filter(gencc_class %in% c("Definitive", "Strong", "Moderate",
                              "Limited", "Absent")) %>%
    left_join(hits, by = c("gene", "phenotype")) %>%
    mutate(lit_score = lit_score_of(.),
           found_specific = as.numeric(found_specific),
           found_general = as.numeric(found_general),
           positive = as.integer(gencc_class != "Absent")) %>%
    attach_ablations()
  if (IS_HAIKU) {
    d <- d %>% select(-lit_score) %>%
      inner_join(haiku_lit(), by = c("gene", "phenotype"))
  }
  d %>%
    filter(!is.na(found_specific), !is.na(found_general),
           !is.na(lit_score), !is.na(knowledge), !is.na(titles),
           !is.na(positive))
}

asc_found_data <- function() {
  hits <- read_tsv(file.path(RESULTS, "asd_found_specific_general.tsv"),
                   show_col_types = FALSE) %>%
    mutate(gene = toupper(gene)) %>%
    select(gene, phenotype, found_specific, found_general)
  curated <- read_tsv(file.path(root, "benchmark", "asd", "pesto.tsv"),
                      show_col_types = FALSE) %>%
    mutate(gene = toupper(gene), lit_score = lit_score_of(.), positive = 1L)
  extra <- not_in_asc_pesto(curated$gene) %>%
    mutate(lit_score = lit_score_of(.), positive = 0L)
  d <- bind_rows(curated, extra) %>%
    left_join(hits, by = c("gene", "phenotype")) %>%
    mutate(found_specific = as.numeric(found_specific),
           found_general = as.numeric(found_general)) %>%
    attach_ablations()
  if (IS_HAIKU) {
    d <- d %>% mutate(gene = toupper(gene)) %>%
      select(-lit_score) %>%
      inner_join(haiku_lit(), by = c("gene", "phenotype"))
  }
  d %>%
    filter(!is.na(found_specific), !is.na(found_general),
           !is.na(lit_score), !is.na(knowledge), !is.na(titles),
           !is.na(positive))
}

panel_limited_absent_bars <- function() {
  # lit_mean, as in panels a and b. The TSV's *_verdict columns are argmax.
  f <- file.path(RESULTS, "limited30_vs_adjacent.tsv")
  if (!file.exists(f)) return(patchwork::plot_spacer())
  d <- read_tsv(f, show_col_types = FALSE)
  long <- bind_rows(
    d %>% transmute(method = "Closed-book",
                    class = factor(set, levels = c("Limited", "Absent")),
                    verdict = as_verdict(mean_call(prior_score))),
    d %>% transmute(method = "pesto-abstracts",
                    class = factor(set, levels = c("Limited", "Absent")),
                    verdict = as_verdict(mean_call(pesto_score)))
  ) %>% filter(!is.na(verdict))
  long$method <- factor(long$method, levels = c("Closed-book", "pesto-abstracts"))
  long %>% count(method, class, verdict) %>%
    group_by(method, class) %>% mutate(pct = 100 * n / sum(n), cls_n = sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = class, y = pct, fill = verdict)) +
    geom_col(width = 0.72) +
    geom_text(aes(x = class, y = 101, label = paste0("n=", cls_n)),
              data = ~ distinct(.x, method, class, cls_n), inherit.aes = FALSE,
              size = 4.2, colour = "grey35", vjust = 0) +
    facet_wrap(~ method, nrow = 1) +
    scale_fill_manual(values = VCOL, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.07))) +
    labs(x = NULL, y = "% of pairs") +
    base_theme +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(size = 15, face = "bold"))
}

limited_adjacent_found_data <- function() {
  read_tsv(file.path(RESULTS, "limited30_vs_adjacent.tsv"),
           show_col_types = FALSE) %>%
    mutate(lit_score = as.numeric(pesto_score),
           found_specific = as.numeric(found_specific),
           found_general = as.numeric(found_general),
           knowledge = as.numeric(prior_score),
           titles = as.numeric(titles_score),
           positive = as.integer(set == "Limited")) %>%
    filter(!is.na(found_specific), !is.na(found_general),
           !is.na(lit_score), !is.na(knowledge), !is.na(titles),
           !is.na(positive))
}

panel_limited_adjacent_found <- function() {
  f <- file.path(RESULTS, "limited30_vs_adjacent.tsv")
  if (!file.exists(f)) return(patchwork::plot_spacer())
  found_auc_panel(limited_adjacent_found_data(),
                  file.path(RESULTS, "limited30_pr_roc_bootstrap.tsv"),
                  roc = FALSE,
                  bar_lab = c(lit_score = "pesto-abstracts",
                              titles = "pesto-titles",
                              knowledge = "Closed-book",
                              found_specific = "n papers (pair)",
                              found_general = "n papers (gene)"))
}

# --------------------------------------------------------------------------
# C and D - Limited pairs dated by oldest cited report, vs Fig. 2a Absents
# --------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

YEAR_LEV <- c("2011-2013", "2014-2016", "2017-2019", "2020-2022", "2023-2025")
YEAR_LAB <- c("2011-2013" = "2011–13", "2014-2016" = "2014–16",
              "2017-2019" = "2017–19", "2020-2022" = "2020–22",
              "2023-2025" = "2023–25")
YEAR_COL <- c(`Model knowledge` = "#7A4E6A", PESTO = "#2E5F4F")

year_bin_of <- function(y) {
  if (y <= 2013) "2011-2013"
  else if (y <= 2016) "2014-2016"
  else if (y <= 2019) "2017-2019"
  else if (y <= 2022) "2020-2022"
  else "2023-2025"
}

lit_score_pts <- function(est, exi, hyp, nov) {
  (4 * as.numeric(est) + 3 * as.numeric(exi) +
   2 * as.numeric(hyp) + 1 * as.numeric(nov)) / 100
}

year_bin_prep <- function() {
  out <- file.path(RESULTS, "prior_year10")
  prior_lim <- read_tsv(file.path(out, "knowledge.tsv"),
                        show_col_types = FALSE) %>%
    mutate(year = as.integer(year),
           score = as.numeric(lit_score),
           method = "Model knowledge",
           positive = 1L) %>%
    filter(!is.na(year), !is.na(score), year >= 2011) %>%
    select(gene, phenotype, year, score, method, positive)

  pesto_raw <- fromJSON(file.path(out, "pesto.json"), simplifyDataFrame = FALSE)
  current_lim <- bind_rows(lapply(pesto_raw, function(r) {
    d <- r$distribution
    tibble(gene = r$gene, phenotype = r$phenotype,
           score = lit_score_pts(d[["Established"]] %||% 0,
                                 d[["Existing"]] %||% 0,
                                 d[["Hypothesized"]] %||% 0,
                                 d[["Novel"]] %||% 0),
           method = "PESTO", positive = 1L)
  })) %>%
    inner_join(prior_lim %>% select(gene, phenotype, year),
               by = c("gene", "phenotype"))

  abs_p30 <- read_tsv(file.path(root, "benchmark", "gencc_clingen_g2p",
                                "pesto.tsv"), show_col_types = FALSE) %>%
    filter(gencc_class == "Absent") %>%
    transmute(gene = toupper(gene), phenotype,
              score = lit_score_of(.),
              method = "PESTO", positive = 0L)
  abs_k30 <- read_tsv(file.path(RESULTS, "ablation_fig2", "knowledge.tsv"),
                      show_col_types = FALSE) %>%
    filter(as.integer(draw) == 1) %>%
    transmute(gene = toupper(gene), phenotype,
              score = as.numeric(lit_score),
              method = "Model knowledge", positive = 0L)
  abs_block <- function(dir) {
    p <- read_tsv(file.path(RESULTS, dir, "pesto.tsv"),
                  show_col_types = FALSE) %>%
      filter(is.na(error) | error == "") %>%
      transmute(gene = toupper(gene), phenotype,
                score = as.numeric(lit_score),
                method = "PESTO", positive = 0L)
    k <- read_tsv(file.path(RESULTS, dir, "knowledge.tsv"),
                  show_col_types = FALSE) %>%
      transmute(gene = toupper(gene), phenotype,
                score = as.numeric(lit_score),
                method = "Model knowledge", positive = 0L)
    list(p = p, k = k)
  }
  e20 <- abs_block("fig2a_absent_extra20")
  e20b <- abs_block("fig2a_absent_extra20b")
  e30 <- abs_block("fig2a_absent_extra30")

  join_abs <- function(k, p) {
    k %>% inner_join(p %>% select(gene, phenotype), by = c("gene", "phenotype"))
  }
  prior_abs <- bind_rows(
    join_abs(abs_k30, abs_p30), join_abs(e20$k, e20$p),
    join_abs(e20b$k, e20b$p), join_abs(e30$k, e30$p)) %>%
    filter(!is.na(score))
  pesto_abs <- bind_rows(
    join_abs(abs_p30, abs_k30), join_abs(e20$p, e20$k),
    join_abs(e20b$p, e20b$k), join_abs(e30$p, e30$k)) %>%
    filter(!is.na(score))
  stopifnot(nrow(prior_abs) == 100, nrow(pesto_abs) == 100)

  lim <- bind_rows(prior_lim, current_lim) %>%
    mutate(bin = factor(vapply(year, year_bin_of, ""), levels = YEAR_LEV),
           method = factor(method, levels = names(YEAR_COL)))
  abs_d <- bind_rows(prior_abs, pesto_abs) %>%
    mutate(method = factor(method, levels = names(YEAR_COL)))
  stopifnot(all(table(lim$bin, lim$method) == 40))

  sc_summ <- lim %>%
    group_by(bin, method) %>%
    summarise(mean = mean(score) - 1,
              se = sd(score) / sqrt(n()),
              .groups = "drop") %>%
    mutate(lo = mean - se, hi = mean + se)
  last2 <- c("2020-2022", "2023-2025")
  fmt_p_score <- function(p) {
    if (p >= 0.001) return(list(label = sprintf("P = %.2f", p), parse = FALSE))
    e <- floor(log10(p))
    m <- p / 10^e
    list(label = sprintf("P == %.1f %%*%% 10^{%d}", m, e), parse = TRUE)
  }
  br <- bind_rows(lapply(names(YEAR_COL), function(m) {
    a <- lim$score[lim$method == m & lim$bin == last2[1]]
    b <- lim$score[lim$method == m & lim$bin == last2[2]]
    fp <- fmt_p_score(t.test(a, b)$p.value)
    tibble(method = factor(m, levels = names(YEAR_COL)),
           x1 = match(last2[1], YEAR_LEV), x2 = match(last2[2], YEAR_LEV),
           lab = fp$label, use_plotmath = fp$parse)
  }))
  br$y <- max((sc_summ %>% filter(bin %in% last2))$hi) + 0.22

  set.seed(2026)
  boot <- bind_rows(lapply(YEAR_LEV, function(b) {
    bind_rows(lapply(names(YEAR_COL), function(m) {
      pos <- lim %>% filter(bin == b, method == m)
      neg <- abs_d %>% filter(method == m)
      sc <- c(pos$score, neg$score)
      yv <- c(pos$positive, neg$positive)
      n <- length(sc)
      tibble(bin = b, method = m,
             ap_pt = pr_curve(sc, yv)$ap,
             ap = replicate(BOOT_B, {
               i <- sample.int(n, n, replace = TRUE)
               pr_curve(sc[i], yv[i])$ap
             }))
    }))
  })) %>%
    mutate(bin = factor(bin, levels = YEAR_LEV),
           method = factor(method, levels = names(YEAR_COL)))
  ci <- boot %>%
    group_by(bin, method) %>%
    summarise(lo = as.numeric(quantile(ap, 0.025)),
              hi = as.numeric(quantile(ap, 0.975)),
              ap = ap_pt[1],
              .groups = "drop")

  last <- "2023-2025"
  pos_p <- lim %>% filter(bin == last, method == "PESTO") %>%
    arrange(gene, phenotype)
  pos_k <- lim %>% filter(bin == last, method == "Model knowledge") %>%
    arrange(gene, phenotype)
  neg_p <- abs_d %>% filter(method == "PESTO") %>%
    arrange(gene, phenotype)
  neg_k <- abs_d %>% filter(method == "Model knowledge") %>%
    arrange(gene, phenotype)
  stopifnot(identical(pos_p$gene, pos_k$gene),
            identical(pos_p$phenotype, pos_k$phenotype),
            identical(neg_p$gene, neg_k$gene),
            identical(neg_p$phenotype, neg_k$phenotype))
  paired_perm_p <- function(a, b, y) {
    obs <- pr_curve(a, y)$ap - pr_curve(b, y)$ap
    set.seed(1)
    ge <- 0L
    n <- length(y)
    for (i in seq_len(PERM_B)) {
      flip <- runif(n) < 0.5
      sa <- ifelse(flip, b, a)
      sb <- ifelse(flip, a, b)
      if ((pr_curve(sa, y)$ap - pr_curve(sb, y)$ap) >= obs) ge <- ge + 1L
    }
    (1 + ge) / (1 + PERM_B)
  }
  p_one <- paired_perm_p(c(pos_p$score, neg_p$score),
                         c(pos_k$score, neg_k$score),
                         c(pos_p$positive, neg_p$positive))
  chance <- 40 / (40 + 100)
  list(lim = lim, sc_summ = sc_summ, br = br, bars = ci,
       chance = chance, p_one = p_one, last = last)
}

panel_year_score <- function(yb) {
  bands <- tibble(
    ymin = c(0, 0.5, 1.5, 2.5),
    ymax = c(0.5, 1.5, 2.5, 3),
    odd  = c(TRUE, FALSE, TRUE, FALSE)
  )
  tick <- 0.09
  ggplot(yb$sc_summ, aes(x = bin, y = mean, fill = method)) +
    geom_rect(data = filter(bands, odd),
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = "#F3F3F3", colour = NA) +
    geom_rect(data = filter(bands, !odd),
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = "#E8E8E8", colour = NA) +
    geom_hline(yintercept = c(0.5, 1.5, 2.5),
               linewidth = 0.25, colour = "grey78") +
    geom_col(width = 0.72) +
    geom_errorbar(aes(ymin = lo, ymax = hi),
                  width = 0.18, linewidth = 0.45, colour = "black") +
    geom_segment(data = yb$br, aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, colour = "grey25", linewidth = 0.4) +
    geom_segment(data = yb$br, aes(x = x1, xend = x1, y = y - tick, yend = y),
                 inherit.aes = FALSE, colour = "grey25", linewidth = 0.4) +
    geom_segment(data = yb$br, aes(x = x2, xend = x2, y = y - tick, yend = y),
                 inherit.aes = FALSE, colour = "grey25", linewidth = 0.4) +
    geom_text(data = filter(yb$br, use_plotmath),
              aes(x = (x1 + x2) / 2, y = y, label = lab),
              inherit.aes = FALSE, parse = TRUE, vjust = -0.45,
              size = 4.0, colour = "grey25") +
    geom_text(data = filter(yb$br, !use_plotmath),
              aes(x = (x1 + x2) / 2, y = y, label = lab),
              inherit.aes = FALSE, parse = FALSE, vjust = -0.45,
              size = 4.0, colour = "grey25") +
    facet_wrap(~ method, ncol = 2) +
    scale_fill_manual(values = YEAR_COL) +
    scale_x_discrete(labels = YEAR_LAB) +
    scale_y_continuous(breaks = c(0.25, 1, 2, 2.75),
                       labels = c("Novel", "Hypoth.",
                                  "Existing", "Estab."),
                       expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = c(0, 3), clip = "off") +
    labs(x = "Year of first report in GenCC", y = "Mean score") +
    base_theme +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(size = 15, face = "bold"),
          axis.ticks.y = element_blank(),
          aspect.ratio = 2,
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.spacing = unit(0.7, "lines"))
}

panel_year_ap <- function(yb) {
  pd <- position_dodge(width = 0.78)
  x_last <- match(yb$last, YEAR_LEV)
  x_left <- x_last - pd$width / 4
  x_right <- x_last + pd$width / 4
  y_br <- max(yb$bars$hi[yb$bars$bin == yb$last]) + 0.08
  tick <- 0.03
  ggplot(yb$bars, aes(x = bin, y = ap, fill = method)) +
    geom_hline(yintercept = yb$chance, linetype = "dashed",
               linewidth = 0.35, colour = "grey60") +
    geom_col(position = pd, width = 0.72) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = pd,
                  width = 0.15, linewidth = 0.45, colour = "grey25") +
    annotate("segment", x = x_left, xend = x_right, y = y_br, yend = y_br,
             colour = "grey25", linewidth = 0.4) +
    annotate("segment", x = x_left, xend = x_left,
             y = y_br - tick, yend = y_br, colour = "grey25", linewidth = 0.4) +
    annotate("segment", x = x_right, xend = x_right,
             y = y_br - tick, yend = y_br, colour = "grey25", linewidth = 0.4) +
    annotate("text", x = x_last, y = y_br,
             label = sprintf("P = %.3f", yb$p_one),
             vjust = -0.45, size = 4.0, colour = "grey25") +
    scale_fill_manual(values = YEAR_COL) +
    scale_x_discrete(labels = YEAR_LAB) +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       expand = expansion(mult = c(0, 0))) +
    coord_cartesian(ylim = c(0, 1.38), clip = "off") +
    labs(x = "Year of first report in GenCC", y = "Average precision") +
    base_theme +
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.background = element_rect(fill = "white", colour = NA),
          aspect.ratio = 1,
          axis.text.x = element_text(angle = 45, hjust = 1))
}

panel_gencc_found <- function() {
  f <- file.path(RESULTS, "gencc_found_specific_general.tsv")
  if (!file.exists(f)) return(patchwork::plot_spacer())
  dest <- if (IS_HAIKU) file.path(HAIKU, "gencc_pr_roc_bootstrap.tsv") else
    file.path(RESULTS, "gencc_pr_roc_bootstrap.tsv")
  found_auc_panel(gencc_found_data(), dest)
}

panel_asc_found <- function() {
  f <- file.path(RESULTS, "asd_found_specific_general.tsv")
  if (!file.exists(f)) return(patchwork::plot_spacer())
  dest <- if (IS_HAIKU) file.path(HAIKU, "asd_pr_roc_bootstrap.tsv") else
    file.path(RESULTS, "asd_pr_roc_bootstrap.tsv")
  found_auc_panel(asc_found_data(), dest)
}

# --------------------------------------------------------------------------
# D - precision and recall, three models
# --------------------------------------------------------------------------
# Shared body: one curve per model, the average precision in the legend, and a
# dashed line at the share of positives, which is the precision of answering
# yes to everything and the only baseline a PR curve has.
pr_panel <- function(d, xlab = "Recall", pad = 0.04) {
  models <- intersect(names(MCOL), unique(d$model))
  parts <- lapply(models, function(m) {
    s <- filter(d, model == m)
    pr <- pr_curve(s$mean_rank, s$positive)
    called <- s$mean_rank >= 2.5
    list(curve = pr$curve %>% mutate(model = m),
         mark = tibble::tibble(
           model = m,
           recall = sum(called & s$positive == 1) / sum(s$positive),
           precision = sum(called & s$positive == 1) / max(sum(called), 1)),
         # Three decimals: on the easier of the two comparisons all three round
         # to 1.00, and the legend would report a tie that is not one.
         label = sprintf("%s  AP = %.3f", m, pr$ap))
  })

  lines <- bind_rows(lapply(parts, `[[`, "curve"))
  marks <- bind_rows(lapply(parts, `[[`, "mark"))
  labels <- setNames(sapply(parts, `[[`, "label"), models)
  chance <- mean(d$positive[d$model == models[1]])

  # The axis does not start at zero. Everything happens in the top few percent,
  # and a full axis would stack the three curves on one line and hide the only
  # difference between them. The baseline stays in frame so the height is still
  # readable against something.
  ymin <- max(0, min(chance, min(lines$precision)) - pad)

  ggplot(lines, aes(x = recall, y = precision, colour = model)) +
    geom_hline(yintercept = chance, linetype = "dashed", linewidth = 0.35,
               colour = "grey60") +
    annotate("text", x = 0.99, y = chance, vjust = -0.5, hjust = 1, size = 3.1,
             colour = "grey50", label = sprintf("chance %.2f", chance)) +
    geom_step(direction = "vh", linewidth = 0.8) +
    geom_point(data = marks, size = 2.4, show.legend = FALSE) +
    scale_colour_manual(values = MCOL, breaks = models, labels = labels) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
    scale_y_continuous(limits = c(ymin, 1), expand = expansion(mult = 0.01)) +
    labs(x = xlab, y = "Precision") +
    base_theme +
    theme(legend.position = c(0.03, 0.30), legend.justification = c(0, 0),
          legend.background = element_blank(),
          legend.text = element_text(size = 10))
}

panel_D <- function() {
  # The weak end of the ladder against pairs that do not exist. Taking all four
  # curated classes puts every model above 0.998 average precision, because a
  # Definitive association and a random recombination are nothing alike, and the
  # panel then says only that the task can be done. Moderate and Limited are the
  # cases where the literature is thin enough for the distinction to be real
  # work, and they are also where panel c shows the models disagreeing.
  #
  # Unbalanced, the negatives being the control panel c already shows, so the
  # baseline sits at 0.61 rather than 0.50.
  f <- file.path(RESULTS, "bench_gencc_pr_scored.tsv")
  if (!file.exists(f)) return(patchwork::plot_spacer())
  read_tsv(f, show_col_types = FALSE) %>%
    filter(gencc %in% c("Moderate", "Limited", "Absent")) %>%
    mutate(positive = as.integer(gencc != "Absent")) %>%
    pr_panel()
}

# Kept for Extended Data: the same three models on fifty OMIM entries against
# fifty recombined pairs, balanced, and drawn from Mendelian entries only.
panel_omim_pr <- function() {
  # Fifty OMIM gene-disease entries against fifty pairs built by recombining
  # those same genes and phenotypes, kept only where OMIM records nothing of the
  # sort. Balanced on purpose: counting how many recorded associations a model
  # recovers rewards one that answers Established to everything, and only pairs
  # that should be refused separate the models.
  #
  # The score is the mean of the four categories weighted 1 to 4 by the points
  # they hold. Rounding it, which is how a user reads one answer, leaves four
  # values and so four operating points; the marks sit at the threshold between
  # Hypothesized and Existing, where that reading calls a pair known.
  f <- file.path(RESULTS, "bench_omim_pr_scored.tsv")
  if (!file.exists(f)) return(patchwork::plot_spacer())
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(positive = as.integer(omim_class == "In OMIM")) %>%
    pr_panel()
}

# --------------------------------------------------------------------------
# The letters are laid on as tags rather than written as panel titles: the text
# refers to Fig. 2c, so they have to survive the titles being dropped.
fig <- if (IS_HAIKU) {
  top <- (panel_C() | panel_B()) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.box.margin = margin(t = -6, b = 0),
          legend.key.size = unit(0.55, "cm"))
  (top / (panel_gencc_found() | panel_asc_found())) +
    plot_layout(widths = c(1, 1), heights = c(1, 1))
} else {
  yb <- year_bin_prep()
  top <- (panel_C() | panel_B()) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.box.margin = margin(t = -6, b = 0),
          legend.key.size = unit(0.55, "cm"))
  (top / (panel_year_score(yb) | panel_year_ap(yb))) +
    plot_layout(widths = c(1, 1), heights = c(1, 1))
}
fig <- fig +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 20),
        plot.tag.position = c(0, 1))

fig2_stem <- if (IS_HAIKU) "figure2_haiku" else "figure2"
ggsave(file.path(FIGURES, paste0(fig2_stem, ".png")), fig,
       width = 11, height = 11.2, dpi = 600, bg = "white")
ggsave(file.path(FIGURES, paste0(fig2_stem, ".pdf")), fig,
       width = 11, height = 11.2, bg = "white")
cat("Wrote", file.path(FIGURES, paste0(fig2_stem, ".png")), "and .pdf\n")
