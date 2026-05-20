#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(png))

read_image_grob <- function(path, label) {
  if (!is.null(path) && file.exists(path)) {
    img <- readPNG(path)
    grob <- rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))
  } else {
    grob <- textGrob(label, gp = gpar(col = "red", cex = 1.5, fontface = "bold"))
  }

  # Add border for visual separation
  grobTree(
      rectGrob(gp = gpar(
        fill = NA,
        col = "NavyBlue",
        lwd = 2,
        lty = "dashed"
      )),
      grob
    )
}

parser <- ArgumentParser(description = "Assemble QC plots into a single panel")
parser$add_argument("--sample_id", required = TRUE, help = "Sample identifier")
parser$add_argument("--knee_plot", required = TRUE, help = "Path to knee plot image")
parser$add_argument("--histogram", required = TRUE, help = "Path to doublet score histogram image")
parser$add_argument("--violin_qc1", required = TRUE, help = "Path to QC1 violin comparison image")
parser$add_argument("--mito_qc2", required = TRUE, help = "Path to mito filtering QC2 image")
parser$add_argument("--soupx_combined", required = TRUE, help = "Path to SoupX combined plot")
parser$add_argument("--umap", required = TRUE, help = "Path to UMAP image")
parser$add_argument("--dotplot", required = TRUE, help = "Path to marker genes dotplot")
parser$add_argument("--output", required = TRUE, help = "Output PNG path")
args <- parser$parse_args()

plots <- list(
  read_image_grob(args$knee_plot, "Missing knee plot"),
  read_image_grob(args$histogram, "Missing histogram"),
  read_image_grob(args$violin_qc1, "Missing violin QC1"),
  read_image_grob(args$mito_qc2, "Missing mito QC2"),
  read_image_grob(args$soupx_combined, "Missing SoupX combined"),
  read_image_grob(args$umap, "Missing UMAP"),
  read_image_grob(args$dotplot, "Missing marker dotplot")
)

################

row1 <- arrangeGrob(
  plots[[1]], plots[[2]],
  ncol = 2,
  widths = c(2, 3)  # effective to this row only
)

row2 <- arrangeGrob(
  plots[[3]], plots[[4]],
  ncol = 2,
  widths = c(2, 1)  # effective to this row only
)

final_grobs <- list(
  row1,               # row1 col1
  row2,               # row1 col2
  plots[[5]], plots[[6]], plots[[7]]
)


dotplot_img <- readPNG(args$dotplot)


################ Layout: lock the dotplot row to its real aspect ratio ################

# 1. Compute the dotplot's true physical height (inches) so it isn't squashed.
dotplot_img <- readPNG(args$dotplot)
ratio <- dim(dotplot_img)[1] / dim(dotplot_img)[2]  # height / width

# Canvas width in inches. Below we render with width=2200 px at res=450,
# so the canvas is 2200/450 inches wide.
canvas_width_in <- 2200 / 450

# If the dotplot is stretched to full canvas width, this is the matching height
# that preserves its 1:1 aspect ratio.
dotplot_height_in <- canvas_width_in * ratio

# 2. Row heights for grid.arrange. unit.c() lets us mix relative ("null") and
# absolute ("inch") units so the dotplot row stays exact while the others
# divvy up the remaining vertical space proportionally.
layout_heights <- unit.c(
  unit(1, "null"),    # Row 1: weight 1
  unit(2, "null"),    # Row 2: weight 2
  unit(1.3, "null"),  # Row 3: weight 1.3
  unit(1.8, "null"),  # Row 4: weight 1.8 (UMAP, usually larger)
  unit(dotplot_height_in, "inch")  # Row 5: pinned to the true dotplot height
)

# 3. Render.
png(args$output, width = 2200, height = 3200, res = 450)

grid.arrange(
  grobs = final_grobs,
  ncol = 1,
  heights = layout_heights
)

dev.off()

cat("Combined QC panel saved to", args$output, "\n")
