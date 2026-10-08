#!/usr/bin/env Rscript
# Extended Data Table 1, third version — the same seven pairs, grouped.
#
# scripts/22_figure_ed1_v2.R spends a whole column on Category, repeating a
# label on every row that only takes three values. Here that label becomes a
# subhead spanning the table, and the width it used goes to the prose. The
# distinction the labels drew inside the Novel/Novel group — no report at all
# against a report too thin to count — moves to the Evidence column, which
# was empty for exactly those two rows.
#
# scripts/22_figure_ed1_v2.R is left in place; both read the same TSV.
#
# Reads results/ed1_v2_pairs.tsv, written by scripts/19_build_ed1_v2_table.py.
# Writes figures/table_ed1_v3.{png,pdf}.
#
# Usage: Rscript scripts/22_figure_ed1_v3.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(grid)
  library(gridExtra)
  library(gtable)
})

root <- normalizePath(file.path(dirname(sub("--file=", "",
         grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
if (is.na(root) || !dir.exists(root)) root <- normalizePath(".")

RESULTS <- file.path(root, "results")
FIGURES <- file.path(root, "figures")
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)

VCOL <- c(Novel = "#b91c1c", Hypothesized = "#b45309", `Hypoth.` = "#b45309",
          Existing = "#1d4ed8", Established = "#15803d")

SIZE      <- 8      # pt, body
SIZE_WHY  <- 7.4    # pt, the prose column
SIZE_HEAD <- 8      # pt, column heads and subheads
PAD_H     <- 2.6    # mm, left and right of every cell
PAD_V     <- 3.6    # mm, total above and below; a third of it is the slack
                    # grid leaves round a line of text, so the visible air
                    # between the last line and the rule is about 0.7 mm
PAD_IN    <- 0.06   # in, white edge around the whole table
SUBHEAD_H <- 6.4    # mm, height of a group subhead row

src <- file.path(RESULTS, "ed1_v2_pairs.tsv")
stopifnot(file.exists(src))
d <- read_tsv(src, show_col_types = FALSE)

wrap <- function(x, width) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

wrap_refs <- function(x, width) {
  vapply(x, function(s) {
    if (is.na(s) || s == "" || s == "-") return(NA_character_)
    parts <- strsplit(s, "; ", fixed = TRUE)[[1]]
    paste(wrap(parts, width), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

short_verdict <- function(x) ifelse(x == "Hypothesized", "Hypoth.", x)

# Same reasoning as v2: plotmath, mantissa quoted so 3.0 keeps its zero, and
# built from the p-value because read_tsv returns p_fmt as a double.
fmt_p <- function(p) {
  e <- floor(log10(p))
  mant <- round(p / 10^e, 1)
  carry <- mant >= 10
  mant[carry] <- mant[carry] / 10
  e[carry] <- e[carry] + 1
  sprintf('"%.1f" %%*%% 10^%d', mant, e)
}

dash <- function(x) ifelse(is.na(x) | x %in% c("", "-", "NA"), "\u2014", x)

# Four patterns. Hypothesized on both is its own group: the pair is claimed
# neither as a published result nor as an Open Targets association, but each
# side has a reason not to call it Novel. The remaining Hypothesized rows
# are Hypothesized in the literature and Novel in Open Targets.
GROUPS <- list(
  list(
    test = function(x) x$lit_verdict == "Novel" & x$ot_verdict == "Novel",
    name = "Novel on both branches"
  ),
  list(
    test = function(x) x$lit_verdict == "Hypothesized" & x$ot_verdict == "Hypothesized",
    name = "Hypothesized on both branches"
  ),
  list(
    test = function(x) x$lit_verdict == "Hypothesized",
    name = "Hypothesized in the literature, Novel in Open Targets"
  ),
  list(
    test = function(x) x$lit_verdict == "Novel" & x$ot_verdict == "Hypothesized",
    name = "Novel in the literature, Hypothesized in Open Targets"
  )
)

grp <- rep(NA_integer_, nrow(d))
for (g in seq_along(GROUPS)) {
  hit <- GROUPS[[g]]$test(d) & is.na(grp)
  grp[hit] <- g
}
stopifnot(!any(is.na(grp)))
d <- d[order(grp, seq_len(nrow(d))), ]
grp <- sort(grp)

# Category's parenthetical is the only place the Novel/Novel rows say why they
# are still Novel, and Evidence is empty on exactly those rows.
paren <- ifelse(grepl("\\([^)]*\\)$", d$category),
                sub("^.*\\(([^)]*)\\)$", "\\1", d$category), NA_character_)
# Not named `evidence`: inside transmute the column of that name would mask it.
evidence_shown <- ifelse(is.na(d$evidence) | d$evidence == "NA",
                         paren, d$evidence)

body <- d %>%
  transmute(
    Gene = gene,
    Phenotype = wrap(phenotype, 16),
    `P (burden)` = fmt_p(burden_pvalue),
    Literature = short_verdict(lit_verdict),
    `Open Targets` = short_verdict(ot_verdict),
    Evidence = wrap(dash(evidence_shown), 16),
    `Why it is worth following up` = wrap(interest, 52),
    `Key reports` = dash(wrap_refs(refs, 24))
  ) %>%
  as.data.frame(stringsAsFactors = FALSE)

# A subhead is a blank row in the matrix, so it adds no width to any column;
# its text is drawn afterwards as one grob spanning every column.
blank <- body[1, ]
blank[] <- ""
rows <- list()
is_sub <- logical(0)
sub_of <- integer(0)
for (g in seq_along(GROUPS)) {
  rows[[length(rows) + 1]] <- blank
  is_sub <- c(is_sub, TRUE)
  sub_of <- c(sub_of, g)
  idx <- which(grp == g)
  rows[[length(rows) + 1]] <- body[idx, ]
  is_sub <- c(is_sub, rep(FALSE, length(idx)))
  sub_of <- c(sub_of, rep(NA_integer_, length(idx)))
}
tbl <- do.call(rbind, rows)

n <- nrow(tbl)
m <- ncol(tbl)
col_gene <- which(names(tbl) == "Gene")
col_p    <- which(names(tbl) == "P (burden)")
col_why  <- which(names(tbl) == "Why it is worth following up")
col_lit  <- which(names(tbl) == "Literature")
col_ot   <- which(names(tbl) == "Open Targets")

# Every per-cell attribute is a matrix passed to the theme, because the grob
# computes its cell sizes from these.
fs <- matrix(SIZE, n, m)
fs[, col_why] <- SIZE_WHY
ink <- matrix("grey15", n, m)
ink[, col_gene] <- "black"
ink[, col_lit] <- unname(VCOL[as.character(tbl$Literature)])
ink[, col_ot]  <- unname(VCOL[as.character(tbl$`Open Targets`)])
face <- matrix("plain", n, m)
face[, col_gene] <- "bold.italic"
face[, col_lit] <- "bold"
face[, col_ot]  <- "bold"
parse_cell <- matrix(FALSE, n, m)
parse_cell[!is_sub, col_p] <- TRUE   # blank subhead rows would fail to parse
parse_head <- rep(FALSE, m)
parse_head[col_p] <- TRUE
names(tbl)[col_p] <- "bold(italic(P))~bold('(burden)')"

tt <- ttheme_default(
  core = list(
    fg_params = list(parse = parse_cell, fontsize = fs, col = ink,
                     fontface = face, hjust = 0, x = 0.02, lineheight = 1.15),
    bg_params = list(fill = NA, col = NA)
  ),
  colhead = list(
    fg_params = list(parse = parse_head, fontsize = SIZE_HEAD, col = "black",
                     fontface = "bold", hjust = 0, x = 0.02),
    bg_params = list(fill = NA, col = NA)
  ),
  padding = unit(c(PAD_H, PAD_V), "mm")
)
gt <- tableGrob(tbl, rows = NULL, theme = tt)

# Cells run from one line to eight, so grid's default vertical centring
# scatters the short entries down the middle of each row. Everything hangs from
# the top instead, on the row's first baseline.
for (i in which(vapply(gt$grobs, inherits, logical(1), "text"))) {
  gt$grobs[[i]]$y <- unit(1, "npc") - unit(PAD_V / 2, "mm")
  gt$grobs[[i]]$vjust <- 1
}

# gt row 1 is the column head; body row i is gt row i + 1.
sub_rows <- which(is_sub) + 1L
gt$heights[sub_rows] <- unit(SUBHEAD_H, "mm")

for (k in seq_along(sub_rows)) {
  g <- sub_of[which(is_sub)][k]
  gt <- gtable_add_grob(gt,
    textGrob(GROUPS[[g]]$name, x = unit(0, "npc"), y = unit(0.9, "mm"),
             hjust = 0, vjust = 0,
             gp = gpar(fontsize = SIZE_HEAD, fontface = "bold", col = "black")),
    t = sub_rows[k], b = sub_rows[k], l = 1, r = ncol(gt),
    z = 90, clip = "off", name = paste0("subhead", k))
}

# Horizontal rules only, in three weights: the frame of the table, the start of
# a group, and the seam between two rows inside one group.
hline <- function(gt, row, y, lwd, col) {
  gtable_add_grob(gt,
    segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                 y0 = unit(y, "npc"), y1 = unit(y, "npc"),
                 gp = gpar(lwd = lwd, col = col, lineend = "butt")),
    t = row, b = row, l = 1, r = ncol(gt), z = 100, clip = "off",
    name = paste0("rule", row, "_", y))
}
gt <- hline(gt, 1, 1, 1.1, "black")   # above the column head
gt <- hline(gt, 1, 0, 0.7, "black")   # under the column head
for (k in seq_along(sub_rows)[-1]) gt <- hline(gt, sub_rows[k], 1, 0.5, "grey45")
for (i in which(!is_sub)) {
  if (i < n && !is_sub[i + 1L]) gt <- hline(gt, i + 1L, 0, 0.3, "grey82")
}
gt <- hline(gt, n + 1L, 0, 1.1, "black")   # foot

# A cell's width is the width of its text on the device that draws it, so the
# canvas is measured on that same device: measuring on pdf() and drawing on
# png() sizes the page to the wrong font metrics.
measure <- function(open_dev) {
  open_dev(nullfile())
  on.exit(dev.off())
  c(convertWidth(sum(gt$widths), "in", valueOnly = TRUE),
    convertHeight(sum(gt$heights), "in", valueOnly = TRUE))
}

draw_table <- function(wh) {
  grid.newpage()
  pushViewport(viewport(x = unit(PAD_IN, "in"), y = unit(PAD_IN, "in"),
                        width = unit(wh[1], "in"), height = unit(wh[2], "in"),
                        just = c("left", "bottom")))
  grid.draw(gt)
  popViewport()
}

render <- function(open_dev, path) {
  wh <- measure(open_dev)
  open_dev(path, width = wh[1] + 2 * PAD_IN, height = wh[2] + 2 * PAD_IN)
  draw_table(wh)
  dev.off()
  cat(sprintf("Wrote %s (%.2f x %.2f in)\n", path,
              wh[1] + 2 * PAD_IN, wh[2] + 2 * PAD_IN))
}

as_png <- function(filename, width = 1, height = 1) {
  png(filename, width = width, height = height, units = "in", res = 600,
      bg = "white", type = "cairo")
}
as_pdf <- function(filename, width = 1, height = 1) {
  cairo_pdf(filename, width = width, height = height, bg = "white")
}

render(as_png, file.path(FIGURES, "table_ed1_v3.png"))
render(as_pdf, file.path(FIGURES, "table_ed1_v3.pdf"))
