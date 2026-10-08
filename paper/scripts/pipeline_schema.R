#!/usr/bin/env Rscript
# The pipeline schema, drawn rather than generated.
#
# The figure carried a raster image produced by a language model. It looked
# right, but nothing in it could be checked against the code, its four evidence
# levels were near the palette of the other panels without matching it, and it
# went into the PDF as pixels. This draws the same diagram from coordinates: the
# levels take their colours from the same constant the rest of the figure uses,
# and the output is vector.
#
# Everything below is placed on one grid. Steps sit on four columns of equal
# width and equal gap, both branches use those columns, and the two rows are
# equidistant from the middle line where they rejoin, so the fork and the join
# are mirror images. Nothing is nudged by eye: change a constant and the whole
# diagram moves with it.
#
# Sourced by scripts/20_figure1.R for panel a. Run directly to write it alone:
#   Rscript scripts/pipeline_schema.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggforce)
  library(grid)
})

# --------------------------------------------------------------------------
# Grid. x runs 0 to 100, y 0 to 56, and coord_fixed holds that ratio. The
# proportion matters: at a flatter shape the boxes shrink against the width of
# the page and the whole diagram reads thin.
# --------------------------------------------------------------------------
CANVAS_H <- 56

BOX_W <- 11.4; BOX_H <- 10       # every workflow step, identical
GAP <- 2.4                       # between steps, identical
COL1 <- 24.75
COLS <- COL1 + (0:3) * (BOX_W + GAP)
CENTRE <- mean(range(c(COLS - BOX_W / 2, COLS + BOX_W / 2)))

ROW_LIT <- 41                    # branch 1 centre line
ROW_DB <- 13                     # branch 2 centre line
MID <- (ROW_LIT + ROW_DB) / 2    # where they fork from and rejoin

IN_X <- 7; IN_W <- 11            # input, wider because its label is
FORK_X <- 15.75                  # vertical stem of the fork
JOIN_X <- 73.4                   # vertical stem of the join

SCALE_X <- c(75.4, 98.4)         # evidence scale
STRIP_W <- 5.8                   # its icon column
SCALE_Y0 <- 5; BAND_H <- 11

HEAD_DY <- 3.8                   # every heading, same distance above its box
HEAD_LIT <- ROW_LIT + BOX_H / 2 + HEAD_DY
HEAD_DB <- ROW_DB + BOX_H / 2 + HEAD_DY
ANNOT_DY <- 1.55                 # every model name, same distance below

# One stroke weight, one arrowhead, one radius, one navy, everywhere.
STROKE <- 0.55
ARROWHEAD <- arrow(length = unit(0.17, "cm"), type = "closed")
RADIUS <- unit(2.5, "pt")
INK <- "#1F3648"
GREEN <- "#3F6B4A"               # branch 2, taken from its database box
MUTED <- "grey40"

FILL <- c(llm = "#DCE6F0", source = "#E4DEEC", db = "#DCE7DB",
          encoder = "#F0E4CC", input = "#F1F1EF")
EDGE <- c(llm = "#61819A", source = "#7C6E96", db = "#5F8562",
          encoder = "#9A7B4A", input = "#9C9C98")

darken <- function(hex, f = 0.74) {
  rgb(t(col2rgb(hex) * f), maxColorValue = 255)
}

# Points of a rectangle, in the order geom_shape wants them.
rect_pts <- function(id, xmin, xmax, ymin, ymax, ...) {
  data.frame(id = id, x = c(xmin, xmax, xmax, xmin),
             y = c(ymin, ymin, ymax, ymax), ..., stringsAsFactors = FALSE)
}

pipeline_schema <- function(vcol) {
  steps <- data.frame(
    x = c(IN_X, COLS, COLS),
    y = c(MID, rep(ROW_LIT, 4), rep(ROW_DB, 4)),
    w = c(IN_W, rep(BOX_W, 8)),
    label = c("(Gene,\nPhenotype)",
              "Synonym\nexpansion", "PubMed\nsearch",
              "Relevance\nfilter", "Novelty\nassessment",
              "Open Targets\nlookup", "Shortlist\nrelated traits",
              "Select\nrelevant traits", "Trait\nmatching"),
    kind = c("input",
             "llm", "source", "llm", "llm",
             "db", "encoder", "llm", "llm"),
    num = c(NA, 1, 2, 3, 4, 5, 6, 7, 8),
    model = c(NA,
              "Claude Haiku", NA, "Claude Haiku", "Claude Opus",
              NA, "BioLORD", "Claude Haiku", "Claude Opus"),
    stringsAsFactors = FALSE)
  steps$xmin <- steps$x - steps$w / 2
  steps$xmax <- steps$x + steps$w / 2
  steps$ymin <- steps$y - BOX_H / 2
  steps$ymax <- steps$y + BOX_H / 2
  boxes <- do.call(rbind, lapply(seq_len(nrow(steps)), function(i) {
    s <- steps[i, ]
    rect_pts(i, s$xmin, s$xmax, s$ymin, s$ymax, kind = s$kind)
  }))

  # Connectors. Every one of them runs from a box edge to a box edge at the
  # centre line of its row, or along a stem; there are no other heights.
  lit <- steps[steps$y == ROW_LIT, ]
  db <- steps[steps$y == ROW_DB, ]
  between <- rbind(
    data.frame(x = head(lit$xmax, -1), xend = tail(lit$xmin, -1), y = ROW_LIT),
    data.frame(x = head(db$xmax, -1), xend = tail(db$xmin, -1), y = ROW_DB))
  between$yend <- between$y

  fork_stem <- data.frame(x = c(IN_X + IN_W / 2, FORK_X), xend = c(FORK_X, FORK_X),
                          y = c(MID, ROW_DB), yend = c(MID, ROW_LIT))
  fork_in <- data.frame(x = FORK_X, xend = c(min(lit$xmin), min(db$xmin)),
                        y = c(ROW_LIT, ROW_DB), yend = c(ROW_LIT, ROW_DB))
  join_stem <- data.frame(
    x = c(max(lit$xmax), max(db$xmax), JOIN_X, JOIN_X),
    xend = c(JOIN_X, JOIN_X, JOIN_X, JOIN_X),
    y = c(ROW_LIT, ROW_DB, ROW_LIT, ROW_DB),
    yend = c(ROW_LIT, ROW_DB, MID, MID))
  join_in <- data.frame(x = JOIN_X, xend = SCALE_X[1], y = MID, yend = MID)

  # The scale. Equal blocks by construction, all four labelled in white.
  bands <- data.frame(
    level = factor(names(vcol), levels = names(vcol)),
    fill = unname(vcol),
    ymin = SCALE_Y0 + BAND_H * (seq_along(vcol) - 1),
    descr = c("No paper links the\ngene to the phenotype",
              "Mechanistic or\nindirect argument only",
              "Human genetic\nevidence reported",
              "Replicated, accepted\ngene-disease pair"),
    stringsAsFactors = FALSE)
  bands$ymax <- bands$ymin + BAND_H
  bands$ink <- "white"
  band_boxes <- do.call(rbind, lapply(seq_len(nrow(bands)), function(i) {
    b <- bands[i, ]
    rbind(rect_pts(i, SCALE_X[1], SCALE_X[2], b$ymin, b$ymax, fill = b$fill),
          rect_pts(i + 100, SCALE_X[1], SCALE_X[1] + STRIP_W, b$ymin, b$ymax,
                   fill = darken(b$fill)))
  }))
  # In place of the icons of the drawn version, four bars filled up to the rank
  # of the level: the same ladder the scale is, said twice.
  bars <- do.call(rbind, lapply(seq_len(nrow(bands)), function(i) {
    b <- bands[i, ]
    data.frame(id = paste(i, 1:4),
               x0 = SCALE_X[1] + (STRIP_W - 3 * 0.9 - 0.62) / 2 + (0:3) * 0.9,
               h = c(1.3, 2.0, 2.7, 3.4), lit = 1:4 <= i,
               base = b$ymin + BAND_H / 2 - 1.7, ink = b$ink,
               stringsAsFactors = FALSE)
  }))
  bar_boxes <- do.call(rbind, lapply(seq_len(nrow(bars)), function(i) {
    r <- bars[i, ]
    rect_pts(r$id, r$x0, r$x0 + 0.62, r$base, r$base + r$h,
             ink = r$ink, lit = r$lit)
  }))
  TEXT_X <- SCALE_X[1] + STRIP_W + 1.6

  keys <- data.frame(x = 1.8 + (0:3) * 17.6,
                     kind = c("llm", "source", "db", "encoder"),
                     lab = c("LLM agent", "Literature source",
                             "External database", "Biomedical\ntext encoder"),
                     stringsAsFactors = FALSE)
  key_boxes <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    rect_pts(i, keys$x[i], keys$x[i] + 3.4, 1.0, 3.4, kind = keys$kind[i])
  }))

  seg <- function(d, head = FALSE) {
    geom_segment(data = d, aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = STROKE, colour = INK,
                 arrow = if (head) ARROWHEAD else NULL)
  }

  ggplot() +
    annotate("text", x = CENTRE, y = HEAD_LIT, label = "BRANCH 1 - LITERATURE",
             size = 5.0, fontface = "bold", colour = INK) +
    annotate("text", x = CENTRE, y = HEAD_DB, label = "BRANCH 2 - OPEN TARGETS",
             size = 5.0, fontface = "bold", colour = INK) +
    annotate("text", x = IN_X, y = MID + BOX_H / 2 + 2.6, label = "INPUT",
             size = 4.2, fontface = "bold", colour = INK) +
    seg(fork_stem) + seg(join_stem) +
    seg(between, head = TRUE) + seg(fork_in, head = TRUE) +
    seg(join_in, head = TRUE) +
    geom_shape(data = boxes, aes(x = x, y = y, group = id, fill = kind,
                                 colour = kind),
               radius = RADIUS, linewidth = 0.8) +
    geom_text(data = steps, aes(x = x, y = y, label = label), size = 3.9,
              lineheight = 0.95, colour = "#16232E") +
    geom_text(data = subset(steps, !is.na(model)),
              aes(x = x, y = ymin - ANNOT_DY, label = model), size = 3.9,
              fontface = "italic", colour = MUTED) +
    # the scale
    annotate("text", x = mean(SCALE_X), y = HEAD_LIT + 1.2,
             label = "EVIDENCE STRENGTH", size = 5.0, fontface = "bold",
             colour = INK) +
    geom_shape(data = band_boxes, aes(x = x, y = y, group = id, fill = I(fill)),
               radius = RADIUS, colour = NA) +
    geom_segment(data = bands[-1, ], aes(x = SCALE_X[1], xend = SCALE_X[2],
                                         y = ymin, yend = ymin),
                 colour = "white", linewidth = 0.9) +
    geom_shape(data = bar_boxes, aes(x = x, y = y, group = id, fill = I(ink),
                                     alpha = I(ifelse(lit, 1, 0.32))),
               radius = unit(0.6, "pt"), colour = NA) +
    geom_text(data = bands, aes(x = TEXT_X, y = ymax - 4.2, label = level,
                                colour = I(ink)),
              hjust = 0, size = 5.4, fontface = "bold") +
    geom_text(data = bands, aes(x = TEXT_X, y = ymin + 3.4, label = descr,
                                colour = I(ink)),
              hjust = 0, size = 4.2, lineheight = 1.05) +
    # the key, set apart the way the drawn version had it
    geom_shape(data = rect_pts(1, 1.2, max(keys$x) + 15.2, 0.15, 4.35),
               aes(x = x, y = y, group = id), fill = NA, colour = "grey78",
               radius = RADIUS, linewidth = 0.4) +
    geom_shape(data = key_boxes, aes(x = x, y = y, group = id, fill = kind,
                                     colour = kind),
               radius = unit(1.5, "pt"), linewidth = 0.5) +
    geom_text(data = keys, aes(x = x + 4.2, y = 2.2, label = lab), hjust = 0,
              size = 4.0, colour = "grey15", lineheight = 0.92) +
    scale_fill_manual(values = FILL, guide = "none") +
    scale_colour_manual(values = EDGE, guide = "none") +
    coord_fixed(ratio = 1, xlim = c(0, 100), ylim = c(0, CANVAS_H),
                expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(4, 10, 2, 2))
}

if (sys.nframe() == 0) {
  VCOL <- c(Novel = "#A9573F", Hypothesized = "#BFA45C",
            Existing = "#5B7F95", Established = "#2E5F4F")
  root <- normalizePath(file.path(dirname(sub("--file=", "",
          grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
  p <- pipeline_schema(VCOL)
  h <- 12.5 * CANVAS_H / 100
  ggsave(file.path(root, "figures", "pipeline_schema.png"), p,
         width = 12.5, height = h, dpi = 320, bg = "white")
  ggsave(file.path(root, "figures", "pipeline_schema.pdf"), p,
         width = 12.5, height = h, bg = "white")
  cat(sprintf("Wrote figures/pipeline_schema.png and .pdf (%.2f x %.2f in)\n",
              12.5, h))
}
