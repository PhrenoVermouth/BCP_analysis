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


################ Core modification: precise height control ################

# 1. Calculate the true physical height (in inches) required for a dot plot
dotplot_img <- readPNG(args$dotplot)
# Get the aspect ratio
ratio <- dim(dotplot_img)[1] / dim(dotplot_img)[2]

# Get the canvas width (in inches) # Your settings are width=2200, res=450, so the width is 2200/450 inches
canvas_width_in <- 2200 / 450

# Calculate how tall (in inches) a dotplot should be if it fills the width # This ensures a 1:1 aspect ratio without distortion
dotplot_height_in <- canvas_width_in * ratio

# 2. Define the height list
# Using unit.c to combine different units
layout_heights <- unit.c(
  unit(1, "null"),   # Row 1: Weight 1 (divide up the remaining space)
  unit(2, "null"), # Row 2: Weight 2
  unit(1.3, "null"), # Row 3: Weight 1.3
  unit(1.8, "null"), # Row 4: Weight 1.8 (UMAP is usually relatively large; allocate more points to it.)
  unit(dotplot_height_in, "inch") # Row 5: Force-lock to the true height in inches
)

# 3. Drawing
png(args$output, width = 2200, height = 3200, res = 450)

grid.arrange(
  grobs = final_grobs,
  ncol = 1,
  heights = layout_heights # A height list using mixed units)

dev.off()

cat("Combined QC panel saved to", args$output, "\n")
