#!/usr/bin/env Rscript
# Figure 3 with CHD, parent-preferring top-5 GO terms in panel c.
# Same style as figures/figure3_ad_r50_nocc_mean_top5.png (mean > 2.5,
# no DD, no left tree, rows in owner blocks sorted by OR). Does not
# overwrite that figure.
#
#   python3 scripts/116_fig3_chd_top5.py
#   Rscript scripts/117_figure3_chd_top5.R
#   Rscript scripts/117_figure3_chd_top5.R --loeuf-decile
#     → figures/figure3_ad_r50_nocc_mean_top5_chd_loeufdec.{png,pdf}
#       panel a = Fisher OR of ≥ Existing in the lowest LOEUF decile
#   Rscript scripts/117_figure3_chd_top5.R --loeuf-decile --brain-decile
#     → ..._chd_loeufdec_braindec.{png,pdf} and figures/figure3.{png,pdf}
#       panel b = Fisher OR of ≥ Existing in the top 10% adult brain-specific
#   Rscript scripts/117_figure3_chd_top5.R --loeuf-decile --brain-decile --two-decile
#     → ..._chd_loeufdec_braindec_top2.{png,pdf}
#       a and b use the top two deciles (20%)
#   Rscript scripts/117_figure3_chd_top5.R --loeuf-decile --brain-decile --or
#     → ..._chd_loeufdec_braindec_or.{png,pdf}  (OR not log2 OR; does not
#       overwrite figures/figure3.*)
LOEUF_DEC  <- "--loeuf-decile" %in% commandArgs(TRUE)
BRAIN_DEC  <- "--brain-decile" %in% commandArgs(TRUE)
TWO_DEC    <- "--two-decile" %in% commandArgs(TRUE)
SHOW_OR    <- "--or" %in% commandArgs(TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(cowplot)
  library(grid)
})

root <- normalizePath(file.path(dirname(sub("--file=", "",
         grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
if (is.na(root) || !dir.exists(root)) root <- normalizePath(".")

DATA      <- file.path(root, "data")
RESULTS   <- file.path(root, "results")
AD_RES    <- file.path(RESULTS, "fig3_chd_nocc_mean_top5")
FIGURES   <- file.path(root, "figures")
LOEUF_SRC <- file.path(DATA, "loeuf_scores.tsv")

source(file.path(root, "scripts", "_panel_dendrogram.R"))

PHENOS <- c("Autism Spectrum Disorder", "epilepsy",
            "schizophrenia", "bipolar disorder",
            "Alzheimer's disease",
            "congenital heart disease", "type 2 diabetes")
PH_SHORT <- c("Autism Spectrum Disorder" = "ASD",
              "epilepsy" = "epilepsy",
              "congenital heart disease" = "CHD",
              "schizophrenia" = "schizophrenia",
              "bipolar disorder" = "bipolar disorder",
              "type 2 diabetes" = "type 2 diabetes",
              "Alzheimer's disease" = "Alzheimer's disease")
SCOL <- c(adult = "#2E5F4F", fetal = "#C47A3A")

base_theme <- theme_classic(base_size = 18) +
  theme(
    axis.text    = element_text(colour = "black", size = 16),
    axis.title   = element_text(size = 17),
    plot.title   = element_text(size = 17, hjust = 0.5),
    plot.tag     = element_text(face = "bold", size = 22),
    legend.title = element_blank(),
    legend.text  = element_text(size = 16)
  )

rho_ci <- function(d, x, y) {
  d %>%
    group_by(phenotype) %>%
    summarise(
      n = n(),
      rho = cor({{x}}, {{y}}, method = "spearman"),
      p = cor.test({{x}}, {{y}}, method = "spearman", exact = FALSE)$p.value,
      .groups = "drop") %>%
    mutate(
      se = 1 / sqrt(n - 3),
      lo = tanh(atanh(rho) - 1.96 * se),
      hi = tanh(atanh(rho) + 1.96 * se),
      phenotype = factor(phenotype, levels = PHENOS)
    )
}

p_label <- function(p) {
  vapply(p, function(pi) {
    if (is.na(pi) || pi >= 0.05) return(NA_character_)
    if (pi >= 0.001) return(sprintf("p == %.3f", pi))
    e <- floor(log10(pi))
    m <- pi / 10^e
    sprintf("p == %.1f%%*%%10^{%d}", m, e)
  }, character(1))
}

rho_plot <- function(d, ylab, colour = NULL, xlab = NULL, or_axis = FALSE) {
  d <- d %>%
    mutate(phenotype = factor(phenotype, levels = rev(PHENOS)),
           plab = p_label(p))
  if ("stage" %in% names(d)) {
    d <- d %>%
      group_by(phenotype) %>%
                    mutate(x_lab = if (or_axis) max(hi) * 10^0.08 else max(hi) + 0.22) %>%
      ungroup()
  } else {
    d <- d %>% mutate(x_lab = if (or_axis) hi * 10^0.08 else hi + 0.18)
  }
  x0 <- if (or_axis) 1 else 0
  x_scale <- if (or_axis) {
    scale_x_continuous(
      trans = "log10",
      breaks = c(0.1, 1, 10, 100),
      labels = c("0.1", "1", "10", "100"),
      expand = expansion(mult = c(0.06, 0.80)))
  } else {
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.80)))
  }
  p <- ggplot(d, aes(x = rho, y = phenotype)) +
    geom_vline(xintercept = x0, linewidth = 0.4, colour = "grey55")
  if (!is.null(colour) && "stage" %in% names(d)) {
    d <- d %>%
      mutate(y_pos = as.numeric(phenotype) +
               ifelse(as.character(stage) == "adult", 0.16, -0.16))
    p <- ggplot(d, aes(x = rho, y = y_pos)) +
      geom_vline(xintercept = x0, linewidth = 0.4, colour = "grey55") +
      geom_errorbar(aes(xmin = lo, xmax = hi, colour = stage),
                    width = 0.12, linewidth = 0.55, orientation = "y") +
      geom_point(aes(colour = stage), size = 2.8) +
      geom_text(aes(label = plab, x = x_lab, colour = stage),
                size = 5, hjust = 0, vjust = 0.5,
                parse = TRUE, na.rm = TRUE, show.legend = FALSE) +
      scale_colour_manual(values = colour) +
      scale_y_continuous(breaks = seq_along(levels(d$phenotype)),
                         labels = PH_SHORT[levels(d$phenotype)],
                         limits = c(0.55, length(levels(d$phenotype)) + 0.45))
  } else {
    p <- p +
      geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.13,
                    linewidth = 0.55, orientation = "y", colour = "grey30") +
      geom_point(size = 2.8, colour = "#2E5F4F") +
      geom_text(aes(label = plab, x = x_lab),
                size = 5, hjust = 0, vjust = 0.5,
                parse = TRUE, na.rm = TRUE, colour = "grey25") +
      scale_y_discrete(labels = PH_SHORT)
  }
  p +
    x_scale +
    coord_cartesian(clip = "off") +
    labs(x = ylab, y = xlab) +
    base_theme +
    theme(legend.position = "top",
          legend.margin = margin(0, 0, -6, 0),
          axis.title.x = element_text(margin = margin(t = 6)))
}

scores <- read_tsv(file.path(RESULTS, "fig3_universe550_meancall_scores.tsv"),
                   show_col_types = FALSE) %>%
  filter(phenotype %in% PHENOS, !is.na(lit_score)) %>%
  transmute(gene = toupper(gene), phenotype, score = as.numeric(lit_score))

loeuf <- read_tsv(LOEUF_SRC, show_col_types = FALSE) %>%
  transmute(gene = toupper(gene_symbol), loeuf = as.numeric(LOEUF)) %>%
  filter(!is.na(loeuf)) %>%
  distinct(gene, .keep_all = TRUE)
joined <- scores %>%
  inner_join(loeuf, by = "gene") %>%
  filter(!is.na(score), !is.na(loeuf))
bins <- joined %>%
  distinct(gene, loeuf) %>%
  mutate(decile = ntile(loeuf, 10))
joined <- joined %>% inner_join(bins %>% select(gene, decile), by = "gene")

a_dat <- joined %>% rho_ci(score, loeuf)
write_tsv(a_dat, file.path(AD_RES, "loeuf_vs_verdict.tsv"))

fisher_or <- function(d, top_deciles = 1) {
  bind_rows(lapply(split(d, d$phenotype), function(sub) {
    top <- sub$decile %in% top_deciles
    reach <- sub$score > 2.5
    tab <- matrix(c(sum(reach & top), sum(!reach & top),
                    sum(reach & !top), sum(!reach & !top)),
                  nrow = 2, byrow = TRUE)
    ft <- fisher.test(tab)
    tibble(
      phenotype = sub$phenotype[1],
      n = nrow(sub),
      n_top = sum(top),
      k_top = sum(reach & top),
      k_rest = sum(reach & !top),
      or = unname(ft$estimate),
      lo = ft$conf.int[1],
      hi = ft$conf.int[2],
      p = ft$p.value,
      rho = log2(unname(ft$estimate)),
      lo_l2 = log2(ft$conf.int[1]),
      hi_l2 = log2(ft$conf.int[2])
    )
  })) %>%
    mutate(phenotype = factor(phenotype, levels = PHENOS))
}

a_or <- fisher_or(joined, top_deciles = if (TWO_DEC) 1:2 else 1)
a_or_file <- if (TWO_DEC) {
  "loeuf_decile_enrichment_top2.tsv"
} else {
  "loeuf_decile_enrichment.tsv"
}
write_tsv(a_or, file.path(AD_RES, a_or_file))

if (LOEUF_DEC) {
  a_plot <- if (SHOW_OR) {
    a_or %>% mutate(rho = or)
  } else {
    a_or %>% mutate(rho = rho, lo = lo_l2, hi = hi_l2)
  }
  panel_a <- rho_plot(
    a_plot,
    if (SHOW_OR) {
      if (TWO_DEC) {
        expression(OR~(phantom()>=~Existing~"in lowest two LOEUF deciles"))
      } else {
        expression(OR~(phantom()>=~Existing~"in lowest LOEUF decile"))
      }
    } else if (TWO_DEC) {
      expression(log[2]~OR~(phantom()>=~Existing~"in lowest two LOEUF deciles"))
    } else {
      expression(log[2]~OR~(phantom()>=~Existing~"in lowest LOEUF decile"))
    },
    or_axis = SHOW_OR
  )
} else {
  panel_a <- rho_plot(
    a_dat,
    expression(rho~"(LOEUF vs PESTO)")
  )
}

gtex <- read.delim(gzfile(file.path(DATA, "gtex_median_tpm.gct.gz")),
                   skip = 2, check.names = FALSE, stringsAsFactors = FALSE)
gtex_tissues <- setdiff(names(gtex), c("Name", "Description"))
gtex_brain   <- grep("^Brain - ", gtex_tissues, value = TRUE)
gtex_other   <- setdiff(gtex_tissues, gtex_brain)
gtex_mat <- as.matrix(gtex[, gtex_tissues])
adult <- tibble(
  gene    = toupper(gtex$Description),
  brain   = apply(gtex_mat[, gtex_brain, drop = FALSE], 1, median),
  other   = apply(gtex_mat[, gtex_other, drop = FALSE], 1, median),
  max_tpm = apply(gtex_mat, 1, max)
) %>%
  group_by(gene) %>% slice_max(max_tpm, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(brain_ratio = log2((brain + 0.1) / (other + 0.1))) %>%
  select(gene, brain_ratio)

fetal_raw <- read_csv(file.path(DATA, "fetal_gene_expression_tissue_with_symbols.csv"),
                      show_col_types = FALSE)
fetal_tissues <- setdiff(names(fetal_raw), c("RowID", "gene_symbol"))
fetal_brain   <- c("Cerebrum", "Cerebellum")
fetal_other   <- setdiff(fetal_tissues, fetal_brain)
fetal_mat <- as.matrix(fetal_raw[, fetal_tissues])
storage.mode(fetal_mat) <- "numeric"
fetal <- tibble(
  gene    = toupper(fetal_raw$gene_symbol),
  brain   = apply(fetal_mat[, fetal_brain, drop = FALSE], 1, median),
  other   = apply(fetal_mat[, fetal_other, drop = FALSE], 1, median),
  max_tpm = apply(fetal_mat, 1, max)
) %>%
  group_by(gene) %>% slice_max(max_tpm, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(brain_ratio = log2((brain + 0.1) / (other + 0.1))) %>%
  select(gene, brain_ratio)

b_dat <- bind_rows(
  scores %>% inner_join(adult, by = "gene") %>% rho_ci(score, brain_ratio) %>%
    mutate(stage = "adult"),
  scores %>% inner_join(fetal, by = "gene") %>% rho_ci(score, brain_ratio) %>%
    mutate(stage = "fetal")
) %>%
  mutate(stage = factor(stage, levels = c("adult", "fetal")))
write_tsv(b_dat, file.path(AD_RES, "expression_vs_verdict_by_stage.tsv"))

brain_or_stage <- function(expr, stage) {
  d <- scores %>%
    inner_join(expr, by = "gene") %>%
    filter(!is.na(score), !is.na(brain_ratio))
  bins <- d %>%
    distinct(gene, brain_ratio) %>%
    mutate(decile = ntile(brain_ratio, 10))
  d <- d %>% inner_join(bins %>% select(gene, decile), by = "gene")
  fisher_or(d, top_deciles = if (TWO_DEC) 9:10 else 10) %>%
    mutate(stage = stage)
}
b_or <- bind_rows(
  brain_or_stage(adult, "adult"),
  brain_or_stage(fetal, "fetal")
) %>%
  mutate(stage = factor(stage, levels = c("adult", "fetal")))
b_or_file <- if (TWO_DEC) {
  "brain_decile_enrichment_top2.tsv"
} else {
  "brain_decile_enrichment.tsv"
}
write_tsv(b_or, file.path(AD_RES, b_or_file))

if (BRAIN_DEC) {
  b_plot <- b_or %>%
    filter(stage == "adult") %>%
    select(-stage)
  if (SHOW_OR) {
    b_plot <- b_plot %>% mutate(rho = or)
  } else {
    b_plot <- b_plot %>% mutate(rho = rho, lo = lo_l2, hi = hi_l2)
  }
  panel_b <- rho_plot(
    b_plot,
    if (SHOW_OR) {
      if (TWO_DEC) {
        expression(OR~(phantom()>=~Existing~"in top 20% brain-specific"))
      } else {
        expression(OR~(phantom()>=~Existing~"in top 10% brain-specific"))
      }
    } else if (TWO_DEC) {
      expression(log[2]~OR~(phantom()>=~Existing~"in top 20% brain-specific"))
    } else {
      expression(log[2]~OR~(phantom()>=~Existing~"in top 10% brain-specific"))
    },
    or_axis = SHOW_OR
  )
} else {
  panel_b <- rho_plot(
    b_dat,
    expression(rho~"(brain specific expr. vs PESTO)"),
    colour = SCOL
  )
}

panel_c <- dendrogram_panel(
  AD_RES,
  heat_file = "go_signatures.tsv",
  corr_file = "pheno_go_correlation.tsv",
  row_mode = "blocks",
  within = "or",
  legend = "bottom",
  as_or = SHOW_OR)

narrow <- theme(plot.margin = margin(4, 30, 4, 14, unit = "mm"))
top <- plot_grid(panel_a + narrow, panel_b + narrow, ncol = 2,
                 align = "h", axis = "tb",
                 labels = c("a", "b"), label_size = 24, label_fontface = "bold")
bot <- plot_grid(panel_c, ncol = 1,
                 labels = "c", label_size = 24, label_fontface = "bold")
fig <- plot_grid(top, bot, ncol = 1,
                 rel_heights = c(1, 1.90), align = "none")

stem <- paste0(
  "figure3_ad_r50_nocc_mean_top5_chd",
  if (LOEUF_DEC) "_loeufdec" else "",
  if (BRAIN_DEC) "_braindec" else "",
  if (TWO_DEC) "_top2" else "",
  if (SHOW_OR) "_or" else "")
ggsave(file.path(FIGURES, paste0(stem, ".png")), fig,
       width = 13.2, height = 16, dpi = 1200, bg = "white")
ggsave(file.path(FIGURES, paste0(stem, ".pdf")), fig,
       width = 13.2, height = 16, bg = "white")
message("wrote ", file.path(FIGURES, paste0(stem, ".png")), " and .pdf")
if (LOEUF_DEC && BRAIN_DEC && !TWO_DEC && !SHOW_OR) {
  file.copy(file.path(FIGURES, paste0(stem, ".png")),
            file.path(FIGURES, "figure3.png"), overwrite = TRUE)
  file.copy(file.path(FIGURES, paste0(stem, ".pdf")),
            file.path(FIGURES, "figure3.pdf"), overwrite = TRUE)
  message("wrote ", file.path(FIGURES, "figure3.png"), " and .pdf")
}
if (LOEUF_DEC) print(as.data.frame(a_or)) else print(as.data.frame(a_dat))
if (BRAIN_DEC) print(as.data.frame(b_or))
