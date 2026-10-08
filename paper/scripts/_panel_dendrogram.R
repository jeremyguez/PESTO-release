#!/usr/bin/env Rscript
# The clustered GO heatmap, built once and used twice: on its own by
# scripts/97_figure_dendrogram.R and as panel c of figure 3 in
# scripts/31_figure_loeuf_expr_go.R.
#
# Column tree: 1 - Spearman rho of log2 odds ratio profiles over every GO
# term with K >= 10, average linkage (scripts/96_pheno_dendrogram.py). The
# eighteen rows are the terms figure 3 already selected; the tree above them
# does not use those eighteen, so the grouping is not circular. A black dot
# marks FDR < 0.05.
#
# source() this file, then call dendrogram_panel(results_dir).
# row_mode = "tree" (default) clusters displayed terms by Spearman of
# their log2 OR profiles. row_mode = "blocks" drops the left tree and
# stacks terms by the phenotype with the highest FDR-significant OR,
# then by q; needs a `source` column (go_signatures.tsv).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
})

LOR_CAP <- 5
ALPHA   <- 0.05

dendro_segments <- function(h) {
  n <- length(h$order)
  pos <- numeric(n - 1)
  seg <- list()
  xat <- function(j) if (j < 0) which(h$order == -j) else pos[j]
  yat <- function(j) if (j < 0) 0 else h$height[j]
  for (i in seq_len(n - 1)) {
    l <- h$merge[i, 1]; r <- h$merge[i, 2]
    xl <- xat(l); xr <- xat(r); y <- h$height[i]
    pos[i] <- (xl + xr) / 2
    seg[[length(seg) + 1]] <- data.frame(
      x = c(xl, xl, xr), xend = c(xl, xr, xr),
      y = c(yat(l), y, y), yend = c(y, y, yat(r)))
  }
  do.call(rbind, seg)
}

block_row_levels <- function(d, col_order, within = "q") {
  if (!"source" %in% names(d)) {
    stop("row_mode='blocks' needs a 'source' column; pass go_signatures.tsv")
  }
  if (!"yid" %in% names(d)) {
    d$yid <- d$label
  }
  owner <- d %>%
    filter(as.character(phenotype) == as.character(source)) %>%
    distinct(yid, label, source, q, lor)
  blocks <- intersect(col_order, unique(as.character(owner$source)))
  visual <- character(0)
  gaps <- character(0)
  for (i in seq_along(blocks)) {
    if (i > 1) {
      g <- paste0(".gap", i)
      visual <- c(visual, g)
      gaps <- c(gaps, g)
    }
    chunk <- owner %>% filter(source == blocks[i])
    if (identical(within, "or")) {
      chunk <- chunk %>% arrange(desc(lor), q, label)
    } else {
      chunk <- chunk %>% arrange(q, label)
    }
    visual <- c(visual, chunk %>% pull(yid))
  }
  list(levels = rev(visual), gaps = gaps,
       labels = setNames(owner$label, owner$yid))
}

dendrogram_panel <- function(results_dir,
                             heat_file = "pheno_heatmap450.tsv",
                             corr_file = "pheno_go_correlation450.tsv",
                             row_mode = "tree",
                             within = "q",
                             legend = "right",
                             as_or = FALSE) {
  corr <- as.matrix(read.delim(file.path(results_dir, corr_file),
                               row.names = 1, check.names = FALSE))
  col_h <- hclust(as.dist(1 - corr), method = "average")
  col_order <- colnames(corr)[col_h$order]

  d <- read_tsv(file.path(results_dir, heat_file),
                show_col_types = FALSE) %>%
    mutate(label = sub(" \\(GO:[0-9]+\\)( CC| BP| MF)?$", "", item))
  if (!"yid" %in% names(d)) {
    d$yid <- d$label
  }

  if (row_mode == "tree") {
    wide <- d %>%
      select(label, phenotype, lor) %>%
      pivot_wider(names_from = phenotype, values_from = lor) %>%
      as.data.frame()
    rownames(wide) <- wide$label
    wide <- as.matrix(wide[, colnames(corr), drop = FALSE])
    wide[is.na(wide)] <- 0
    row_h <- hclust(as.dist(1 - cor(t(wide), method = "spearman")),
                    method = "average")
    row_order <- rownames(wide)[row_h$order]
    d$yid <- d$label
    y_scale <- scale_y_discrete(expand = c(0, 0), position = "right")
  } else if (row_mode == "blocks") {
    blocked <- block_row_levels(d, col_order, within = within)
    row_order <- blocked$levels
    y_scale <- scale_y_discrete(
      expand = c(0, 0), position = "right", drop = FALSE,
      labels = function(x) {
        ifelse(grepl("^\\.gap", x), "",
               unname(blocked$labels[x]))
      })
  } else {
    stop("row_mode must be 'tree' or 'blocks'")
  }

  d <- d %>%
    mutate(phenotype = factor(phenotype, levels = col_order),
           yid       = factor(yid, levels = row_order),
           lor_c     = pmax(pmin(lor, LOR_CAP), -LOR_CAP),
           sig       = q < ALPHA)

  fill_guide <- if (identical(legend, "bottom")) {
    guide_colourbar(barwidth = 12, barheight = 0.7, title.position = "left",
                    title.vjust = 1, direction = "horizontal")
  } else {
    guide_colourbar(barheight = 7)
  }
  or_breaks <- log2(c(0.1, 1, 10))
  or_labels <- c("0.1", "1", "10")
  heat <- ggplot(d, aes(x = phenotype, y = yid)) +
    geom_tile(aes(fill = lor_c), colour = "white", linewidth = 0.6) +
    geom_point(data = filter(d, sig), size = 1.5, colour = "grey10") +
    scale_fill_gradient2(
      low = "#2C6E8F", mid = "grey96", high = "#B33A3A", midpoint = 0,
      limits = c(-LOR_CAP, LOR_CAP),
      name = if (as_or) "odds ratio" else expression(log[2]~"odds ratio"),
      breaks = if (as_or) or_breaks else waiver(),
      labels = if (as_or) or_labels else waiver(),
      guide = fill_guide,
      na.value = "white") +
    scale_x_discrete(expand = c(0, 0)) +
    y_scale +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 17) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(colour = "black", face = "bold", size = 16),
          axis.ticks.length.y.right = unit(3, "mm"),
          axis.text.y.right = element_text(colour = "black", size = 15, hjust = 0,
                                           margin = margin(l = 4, unit = "mm")),
          legend.position = legend,
          legend.direction = if (identical(legend, "bottom")) "horizontal" else "vertical",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 14),
          plot.margin = margin(2, 5, 2, 2, unit = "mm"))

  bare <- theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
  top <- ggplot(dendro_segments(col_h)) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = 0.5, colour = "grey25", lineend = "square") +
    scale_x_continuous(limits = c(0.5, ncol(corr) + 0.5), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    bare

  if (row_mode == "blocks") {
    return(top / heat + plot_layout(heights = c(1, 5.5)))
  }

  left <- ggplot(dendro_segments(row_h)) +
    geom_segment(aes(x = y, xend = yend, y = x, yend = xend),
                 linewidth = 0.5, colour = "grey25", lineend = "square") +
    scale_y_continuous(limits = c(0.5, nrow(wide) + 0.5), expand = c(0, 0)) +
    scale_x_reverse(expand = expansion(mult = c(0.04, 0))) +
    bare

  (plot_spacer() + top + plot_layout(widths = c(0.6, 5))) /
    (left + heat + plot_layout(widths = c(0.6, 5))) +
    plot_layout(heights = c(1, 5.5))
}
