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
    rectGrob(gp = gpar(fill = NA, col = "grey50", lwd = 2)),
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

layout <- rbind(
  c(1, 2),  # knee + histogram
  c(3, 4),  # doublet removal vs QC2 (vertical) side-by-side
  c(5, 5),  # SoupX
  c(6, 6),  # UMAP + bars
  c(7, 7)   # marker dotplot
)

png(args$output, width = 2200, height = 3200, res = 220)
  grid.arrange(grobs = plots, layout_matrix = layout, widths = c(1, 1))
  grid.text(paste0("QC Overview - ", args$sample_id),
            x = unit(0.5, "npc"), y = unit(0.99, "npc"),
            gp = gpar(fontface = "bold", cex = 1.2))
dev.off()

cat("Combined QC panel saved to", args$output, "\n")
