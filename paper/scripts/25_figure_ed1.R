#!/usr/bin/env Rscript
# Extended Data Figure 1 - where the Novel calls in Figure 1b come from.
#
#   a  PubMed records retrieved, by literature verdict. Novel is not an empty
#      search: only seven of 252 Novel pairs retrieved nothing.
#   b  literature verdicts by phenotype (n >= 8 after merging Height)
#   c  the same phenotypes, same order, under Max. Open Targets recovers
#      prior GWAS evidence on traits such as height; clinical traits with
#      little rare-variant literature stay Novel.
#
# Reads results/all_runs_auto_v7.tsv, the table behind figure1_v7, so that
# panel c agrees with Figure 1b. Writes figures/figure_ed1.{png,pdf}.
#
# Usage: Rscript scripts/25_figure_ed1.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
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

MIN_PAIRS <- 8
VERDICTS <- c("Novel", "Hypothesized", "Existing", "Established")
VCOL <- c(Novel = "#A9573F", Hypothesized = "#BFA45C",
          Existing = "#5B7F95", Established = "#2E5F4F")

# The two studies name some phenotypes differently; each pair keeps one row
# (scripts/11_build_all_runs_auto.py), and both names are drawn as one.
aliases <- read_tsv(file.path(root, "data", "phenotype_aliases.tsv"),
                    show_col_types = FALSE)

canon_pheno <- function(p) {
  p <- trimws(as.character(p))
  key <- trimws(gsub("[^a-z0-9]+", " ", tolower(p)))
  key <- ifelse(key %in% aliases$alias,
                aliases$canonical[match(key, aliases$alias)], key)
  if (key %in% aliases$canonical) return(aliases$label[match(key, aliases$canonical)])
  if (tolower(p) == "height") return("Height")
  if (tolower(p) == "total cholesterol") return("Total cholesterol")
  if (p == "hip-circumference-mean") return("Hip circumference")
  if (p == "mean corpuscular volume") return("Mean corpuscular volume")
  if (p == "red cell distribution width") return("Red cell distribution width")
  if (p == "mean platelet volume") return("Mean platelet volume")
  p
}

rank_of <- function(x) match(as.character(x), VERDICTS)

d <- read_tsv(file.path(RESULTS, "all_runs_auto_v7.tsv"), show_col_types = FALSE) %>%
  filter(source %in% c("AoU", "BRAVA"), verdict %in% VERDICTS,
         ot_verdict %in% VERDICTS) %>%
  mutate(verdict = factor(verdict, levels = VERDICTS),
         ot_verdict = factor(ot_verdict, levels = VERDICTS),
         max_verdict = factor(VERDICTS[pmax(rank_of(verdict), rank_of(ot_verdict))],
                              levels = VERDICTS),
         phenotype = vapply(phenotype, canon_pheno, character(1)),
         n_pubmed = as.numeric(n_pubmed))

base_theme <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0),
    plot.subtitle = element_text(size = 12, colour = "grey30"),
    axis.text = element_text(colour = "black", size = 13),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    legend.key.size = unit(0.52, "cm")
  )

n_novel <- sum(d$verdict == "Novel")
n_empty <- sum(d$verdict == "Novel" & d$n_pubmed == 0)
med_n <- median(d$n_pubmed[d$verdict == "Novel"], na.rm = TRUE)

pA <- ggplot(d, aes(x = verdict, y = n_pubmed, fill = verdict)) +
  geom_violin(trim = TRUE, colour = NA, alpha = 0.9, scale = "width") +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.14, seed = 1),
             size = 0.28, alpha = 0.22, colour = "grey25") +
  scale_fill_manual(values = VCOL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = "a",
       subtitle = paste0("Only ", n_empty, " of ", n_novel,
                         " Novel pairs retrieved no PubMed records ",
                         "(median ", round(med_n), ")"),
       x = NULL, y = "PubMed records retrieved") +
  base_theme +
  theme(axis.text.x = element_text(angle = 18, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 2)))

keep <- d %>% count(phenotype) %>% filter(n >= MIN_PAIRS)
covered <- sum(keep$n)

order_lv <- d %>%
  filter(phenotype %in% keep$phenotype) %>%
  count(phenotype, verdict) %>%
  group_by(phenotype) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  filter(verdict == "Novel") %>%
  select(phenotype, novel_pct = pct) %>%
  right_join(keep, by = "phenotype") %>%
  mutate(novel_pct = coalesce(novel_pct, 0)) %>%
  arrange(novel_pct) %>%
  pull(phenotype)

make_bars <- function(verdict_col, title, subtitle, show_n, show_y) {
  ph <- d %>%
    filter(phenotype %in% keep$phenotype) %>%
    count(phenotype, verdict = .data[[verdict_col]]) %>%
    complete(phenotype, verdict = VERDICTS, fill = list(n = 0)) %>%
    group_by(phenotype) %>%
    mutate(tot = sum(n), pct = 100 * n / tot) %>%
    ungroup() %>%
    mutate(phenotype = factor(phenotype, levels = order_lv),
           verdict = factor(verdict, levels = VERDICTS))
  labels <- ph %>%
    distinct(phenotype, tot) %>%
    mutate(phenotype = factor(phenotype, levels = order_lv))
  p <- ggplot(ph, aes(x = pct, y = phenotype, fill = verdict)) +
    geom_col(width = 0.78) +
    scale_fill_manual(values = VCOL) +
    scale_x_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 112),
                       expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    labs(title = title, subtitle = subtitle,
         x = "% of associations", y = NULL) +
    base_theme +
    theme(axis.text.y = element_text(size = 11),
          plot.margin = margin(4, 12, 4, 4),
          legend.position = "none")
  if (show_n) {
    p <- p +
      geom_text(data = labels, aes(x = 102, y = phenotype, label = tot),
                inherit.aes = FALSE, hjust = 0, size = 3.4, colour = "grey35") +
      annotate("text", x = 102, y = length(order_lv) + 0.85, label = "n",
               hjust = 0, size = 3.4, colour = "grey35", fontface = "italic")
  }
  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  p
}

pB <- make_bars("verdict", "b", "Literature", TRUE, TRUE) +
  theme(legend.position = "bottom")
pC <- make_bars("max_verdict", "c", "Max", FALSE, FALSE)

fig <- free(pA, type = "label", side = "l") / (pB | pC) +
  plot_layout(heights = c(0.85, 1.65), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(FIGURES, "figure_ed1.png"), fig,
       width = 10.8, height = 10.6, dpi = 600, bg = "white")
ggsave(file.path(FIGURES, "figure_ed1.pdf"), fig,
       width = 10.8, height = 10.6, bg = "white")
cat("Wrote", file.path(FIGURES, "figure_ed1.png"), "and .pdf\n")
cat("Phenotypes in b/c:", length(order_lv), "covering", covered, "of", nrow(d), "\n")
cat("Empty-search Novel:", n_empty, "/", n_novel, "\n")
