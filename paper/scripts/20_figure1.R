#!/usr/bin/env Rscript
# Figure 1 - the pipeline and what it produces on real association results.
#
# Layout:
#   A A
#   B B
#   C D
#
#   A  pipeline schema, with the four-level verdict ladder
#   B  literature, Open Targets, and the stronger of the two, for both
#      cohorts and for the negative-control permutation
#   C  literature verdict against Open Targets verdict
#   D  burden significance against verdict
#
# The table of highest-confidence novel associations left the main text and
# is Extended Data Table 1, not a figure. Do not write figure_ed3 from here:
# that slot is now the GO+PEPPER bars (scripts/33_figure_ed3.R).
#
# Validation against external references is Figure 2.
#
# Usage: Rscript scripts/20_figure1.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggpubr)
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
# The palette of Figure 2: desaturated earth and slate, warm at the unsupported
# end and cool at the settled one, spread far enough in lightness to survive
# greyscale and red-green colour blindness. The two figures show the same four
# categories and had no business showing them in two different sets of colours.
VCOL <- c(Novel = "#A9573F", Hypothesized = "#BFA45C",
          Existing = "#5B7F95", Established = "#2E5F4F")

base_theme <- theme_classic(base_size = 15) +
  theme(
    plot.title    = element_text(face = "bold", size = 18, hjust = 0),
    axis.text     = element_text(colour = "black", size = 15),
    axis.title    = element_text(size = 16),
    legend.title  = element_blank(),
    legend.text   = element_text(size = 13),
    legend.key.size = unit(0.58, "cm")
  )

# The automatic abstracts / abstracts-trait cohort run, with the v7 Open Targets
# match prompt. FIG1_RUNS, FIG1_ABSENT and FIG1_OUT point the same figure at a
# variant table without editing this script; unset, they draw the published
# figure.
fig1_runs <- Sys.getenv("FIG1_RUNS", file.path(RESULTS, "all_runs_auto_v7.tsv"))
d <- read_tsv(fig1_runs, show_col_types = FALSE) %>%
  filter(source %in% c("AoU", "BRAVA"), verdict %in% VERDICTS) %>%
  mutate(verdict = factor(verdict, levels = VERDICTS),
         ot_verdict = factor(ot_verdict, levels = VERDICTS),
         p = suppressWarnings(as.numeric(burden_pvalue)),
         source = factor(source, levels = c("AoU", "BRAVA")))

# Panel b only. The table assigns a shared hit to AoU; the BRaVa bar should
# still show it. AoU's "BMI" is BRaVa's "Body mass index". Panels c and d
# keep the unique rows.
norm_ph <- function(x) {
  p <- trimws(gsub("[^a-z0-9]+", " ", tolower(x)))
  ifelse(p == "bmi", "body mass index", p)
}
aou_src <- read_tsv(file.path(root, "data", "raw", "AoU_results.tsv"),
                    show_col_types = FALSE) %>%
  transmute(gene = toupper(gene_symbol), ph = norm_ph(description)) %>%
  filter(gene != "", !is.na(ph), ph != "") %>% distinct()
brava_raw <- read_tsv(file.path(root, "data", "raw", "Duncan_results.tsv"),
                      show_col_types = FALSE)
brava_src <- bind_rows(
    brava_raw %>% transmute(gene = toupper(external_gene_name),
                            ph = norm_ph(phenotype_full)),
    brava_raw %>% transmute(gene = toupper(external_gene_name),
                            ph = norm_ph(phenotype))) %>%
  filter(gene != "", !is.na(ph), ph != "") %>% distinct()
shared <- inner_join(aou_src, brava_src, by = c("gene", "ph"))
d_key <- d %>% mutate(gene_u = toupper(gene), ph = norm_ph(phenotype))
already_brava <- d_key %>% filter(source == "BRAVA") %>% distinct(gene_u, ph)
to_brava <- d_key %>%
  filter(source == "AoU") %>%
  inner_join(shared, by = c("gene_u" = "gene", "ph" = "ph")) %>%
  anti_join(already_brava, by = c("gene_u", "ph")) %>%
  mutate(source = factor("BRAVA", levels = c("AoU", "BRAVA"))) %>%
  select(-gene_u, -ph)
cat("Panel b: copied", nrow(to_brava), "AoU hits into the BRaVa bars\n")

# --------------------------------------------------------------------------
# A - pipeline schema
# --------------------------------------------------------------------------
source(file.path(root, "scripts", "pipeline_schema.R"))

panel_A <- function() {
  pipeline_schema(VCOL) +
    labs(title = "a") +
    theme(plot.title = element_text(face = "bold", size = 18, hjust = 0))
}

# --------------------------------------------------------------------------
# B - literature, Open Targets, and the stronger of the two
# --------------------------------------------------------------------------
METHODS <- c("Literature", "Open Targets", "Max")
GROUPS <- c("AoU", "BRAVA", "Absent")
GROUP_TITLE <- c(AoU = "AoU + Genebass", BRAVA = "BRaVa",
                 Absent = "Negative controls")
# Three groups of three bars, one empty slot between groups.
GROUP_ORIGIN <- c(AoU = 1L, BRAVA = 5L, Absent = 9L)
rank_of <- function(x) match(as.character(x), VERDICTS)

db <- bind_rows(d, to_brava) %>%
  filter(ot_verdict %in% VERDICTS) %>%
  mutate(max_verdict = VERDICTS[pmax(rank_of(verdict), rank_of(ot_verdict))])
b_cohort <- bind_rows(
    db %>% transmute(source, method = "Literature", verdict),
    db %>% transmute(source, method = "Open Targets", verdict = ot_verdict),
    db %>% transmute(source, method = "Max", verdict = max_verdict))

perm_f <- Sys.getenv("FIG1_ABSENT",
                     file.path(RESULTS, "bench_fig1c_absent50_hgnc_pesto_v7.tsv"))
perm <- read_tsv(perm_f, show_col_types = FALSE) %>%
  filter(verdict %in% VERDICTS, open_targets_verdict %in% VERDICTS) %>%
  mutate(max_verdict = VERDICTS[pmax(rank_of(verdict), rank_of(open_targets_verdict))])
b_absent <- bind_rows(
    perm %>% transmute(source = "Absent", method = "Literature", verdict),
    perm %>% transmute(source = "Absent", method = "Open Targets",
                       verdict = open_targets_verdict),
    perm %>% transmute(source = "Absent", method = "Max", verdict = max_verdict))

b_x <- c(Literature = 0L, `Open Targets` = 1L, Max = 2L)
b_counts <- bind_rows(b_cohort, b_absent) %>%
  mutate(source = factor(source, levels = GROUPS),
         method = factor(method, levels = METHODS),
         verdict = factor(verdict, levels = VERDICTS)) %>%
  count(source, method, verdict) %>%
  complete(source, method, verdict, fill = list(n = 0)) %>%
  group_by(source, method) %>%
  mutate(pct = 100 * n / sum(n), N = sum(n)) %>%
  ungroup() %>%
  mutate(x = b_x[as.character(method)] + GROUP_ORIGIN[as.character(source)])
b_n <- b_counts %>% distinct(source, N)
b_heads <- data.frame(
  x = unname(GROUP_ORIGIN[GROUPS]) + 1,
  label = paste0(GROUP_TITLE[GROUPS], " (n = ",
                 b_n$N[match(GROUPS, as.character(b_n$source))], ")"))

pB <- ggplot(b_counts, aes(x = x, y = pct, fill = verdict)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = ifelse(pct >= 6, n, "")),
            position = position_stack(vjust = 0.5),
            size = 4.0, colour = "white", fontface = "bold") +
  geom_vline(xintercept = c(4, 8), linetype = "dotted", colour = "grey35",
             linewidth = 0.55) +
  geom_text(data = b_heads, aes(x = x, y = 108, label = label),
            fontface = "bold", size = 4.4, inherit.aes = FALSE) +
  scale_fill_manual(values = VCOL) +
  scale_x_continuous(breaks = c(1:3, 5:7, 9:11),
                     labels = rep(METHODS, 3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     limits = c(0, 114)) +
  coord_cartesian(ylim = c(0, 114), clip = "off") +
  labs(title = "b",
       x = NULL, y = "% of associations") +
  base_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 22, hjust = 1),
        plot.margin = margin(t = 8, r = 4, b = 4, l = 6))

# --------------------------------------------------------------------------
# C - literature verdict against Open Targets verdict
# --------------------------------------------------------------------------
mat <- d %>% filter(!is.na(ot_verdict)) %>% count(verdict, ot_verdict, name = "n") %>%
  complete(verdict, ot_verdict, fill = list(n = 0))
nmax <- max(mat$n)

pC <- ggplot(mat, aes(x = ot_verdict, y = verdict, fill = n)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = n, colour = n > nmax * 0.5), size = 4.8, fontface = "bold") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#1e293b"), guide = "none") +
  # A sequential ramp built on the slate of the verdict palette, so the panel
  # reads as part of the same figure. It codes a count rather than a category,
  # which is why it is a ramp and not four colours.
  scale_fill_gradientn(colours = c("#F4F2EE", "#D9E1E5", "#AFC3CC", "#7B9AA9",
                                   "#40606F"),
                       name = "pairs") +
  scale_y_discrete(limits = VERDICTS) +
  labs(title = "c",
       x = "Open Targets verdict", y = "Literature verdict") +
  base_theme +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "right",
        legend.title = element_text(size = 11),
        panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank())

# --------------------------------------------------------------------------
# D - burden significance against Max verdict
# --------------------------------------------------------------------------
dc <- d %>%
  filter(!is.na(p), p > 0, ot_verdict %in% VERDICTS) %>%
  mutate(max_verdict = factor(
           VERDICTS[pmax(rank_of(verdict), rank_of(ot_verdict))],
           levels = VERDICTS),
         logp = -log10(p))

pE <- ggplot(dc, aes(x = max_verdict, y = logp, fill = max_verdict)) +
  geom_violin(trim = TRUE, alpha = 0.92, colour = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.15, seed = 1),
             size = 0.3, alpha = 0.22, colour = "grey25") +
  scale_fill_manual(values = VCOL, guide = "none") +
  scale_y_continuous(trans = "log1p", breaks = c(5, 10, 30, 100, 300),
                     limits = c(5, 320)) +
  labs(title = "d",
       x = NULL, y = expression(-log[10](italic(P))~"burden")) +
  base_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1))

# --------------------------------------------------------------------------
# Extended Data 3 - highest-confidence novel associations
# --------------------------------------------------------------------------
# A table is a poor use of a main-text panel: it is read row by row while the
# --------------------------------------------------------------------------
# coord_fixed on a 100 x 56 canvas: the schema fills the page only if its
# row is tall enough. 8 : 4 : 4 at 16.6 in gives it that room.
fig <- free(panel_A()) / free(pB) / (pC | pE) +
  plot_layout(heights = c(8, 4, 4))

fig1_out <- Sys.getenv("FIG1_OUT", "figure1_v7")
ggsave(file.path(FIGURES, paste0(fig1_out, ".png")), fig, width = 12.5, height = 16.6, dpi = 600, bg = "white")
ggsave(file.path(FIGURES, paste0(fig1_out, ".pdf")), fig, width = 12.5, height = 16.6, bg = "white")
cat("Wrote", file.path(FIGURES, paste0(fig1_out, ".png")), "and .pdf\n")

cat("\nCounts used in panel b:\n")
print(bind_rows(d, to_brava) %>% count(source, verdict) %>%
        pivot_wider(names_from = verdict, values_from = n))
cat("\nDouble-novel associations (unique pairs):",
    sum(d$verdict == "Novel" & d$ot_verdict == "Novel", na.rm = TRUE), "\n")
