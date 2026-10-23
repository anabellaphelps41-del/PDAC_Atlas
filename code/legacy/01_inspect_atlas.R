#!/usr/bin/env Rscript
# 01_inspect_atlas.R
# Load the full PDAC atlas once and report everything needed to plan subsetting.
# Nothing is modified. Outputs a text report + a compressed metadata CSV.
#
# Usage: Rscript 01_inspect_atlas.R <path/to/atlas.rds> <outdir>

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args    <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript 01_inspect_atlas.R <atlas.rds> <outdir>")
obj_path <- args[1]
outdir   <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

sv <- packageVersion("Seurat")
cat("=== ENVIRONMENT ===\n")
cat("R:              ", R.version.string, "\n")
cat("Seurat:         ", as.character(sv), "\n")
cat("SeuratObject:   ", as.character(packageVersion("SeuratObject")), "\n")
cat("Matrix:         ", as.character(packageVersion("Matrix")), "\n\n")

cat("=== LOADING ===\n")
cat("file:  ", obj_path, "\n")
cat("size:  ", round(file.size(obj_path) / 1e9, 2), "GB on disk\n")
t0 <- Sys.time()
obj <- readRDS(obj_path)
cat("readRDS wall time: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min\n")
cat("in-memory size:    ", round(as.numeric(object.size(obj)) / 1e9, 2), " GB\n")
cat("R heap in use:     ", round(sum(gc()[, 2]) / 1024, 2), " GB\n\n")

cat("=== OBJECT ===\n")
cat("class: ", paste(class(obj), collapse = ", "), "\n")
print(obj)
cat("\ncells: ", ncol(obj), "\n")
cat("features (default assay): ", nrow(obj), "\n\n")

cat("=== ASSAYS ===\n")
print(Assays(obj))
cat("default assay: ", DefaultAssay(obj), "\n")
for (a in Assays(obj)) {
  cat("\n-- assay:", a, "--\n")
  ao <- obj[[a]]
  cat("  class: ", paste(class(ao), collapse = ", "), "\n")
  cat("  features: ", nrow(ao), "\n")
  # Seurat v5 uses layers; v4 uses slots
  if (inherits(ao, "Assay5")) {
    cat("  layers: ", paste(SeuratObject::Layers(ao), collapse = ", "), "\n")
  } else {
    for (sl in c("counts", "data", "scale.data")) {
      d <- tryCatch(dim(methods::slot(ao, sl)), error = function(e) NULL)
      cat("  slot ", sl, ": ",
          if (is.null(d) || any(d == 0)) "empty" else paste(d, collapse = " x "), "\n", sep = "")
    }
  }
}

cat("\n=== REDUCTIONS ===\n")
print(Reductions(obj))
for (r in Reductions(obj)) {
  cat("  ", r, ": ", paste(dim(Embeddings(obj, r)), collapse = " x "), "\n", sep = "")
}

cat("\n=== GRAPHS ===\n")
print(Graphs(obj))

cat("\n=== METADATA COLUMNS ===\n")
md <- obj@meta.data
cat("dim: ", paste(dim(md), collapse = " x "), "\n\n")
info <- data.frame(
  column   = colnames(md),
  class    = vapply(md, function(x) class(x)[1], character(1)),
  n_unique = vapply(md, function(x) length(unique(x)), integer(1)),
  n_na     = vapply(md, function(x) sum(is.na(x)), integer(1)),
  stringsAsFactors = FALSE
)
print(info, row.names = FALSE)

cat("\n=== CATEGORICAL LEVELS (<= 150 unique) ===\n")
cat("Look here for your hierarchical annotation columns.\n")
for (cn in colnames(md)) {
  v <- md[[cn]]
  if (!(is.character(v) || is.factor(v) || is.logical(v))) next
  v <- as.character(v)
  u <- unique(v)
  if (length(u) > 150) {
    cat("\n== ", cn, " : ", length(u), " levels (too many, skipped) ==\n", sep = "")
    next
  }
  cat("\n== ", cn, " (", length(u), " levels) ==\n", sep = "")
  print(sort(table(v), decreasing = TRUE))
}

cat("\n=== NUMERIC SUMMARIES ===\n")
for (cn in colnames(md)) {
  v <- md[[cn]]
  if (!is.numeric(v)) next
  cat("\n== ", cn, " ==\n", sep = "")
  print(summary(v))
}

# Small, downloadable artifact: metadata only. No expression data.
meta_out <- file.path(outdir, "atlas_metadata.csv.gz")
cat("\n=== WRITING METADATA ===\n")
md$cell_barcode <- rownames(md)
write.csv(md, gzfile(meta_out), row.names = FALSE)
cat("wrote: ", meta_out, " (", round(file.size(meta_out) / 1e6, 1), " MB)\n", sep = "")

# Also dump the level tables in machine-readable form
lev_out <- file.path(outdir, "annotation_levels.txt")
sink(lev_out)
for (cn in colnames(md)) {
  v <- md[[cn]]
  if (!(is.character(v) || is.factor(v))) next
  u <- unique(as.character(v))
  if (length(u) > 150) next
  cat("### ", cn, "\n", sep = "")
  cat(paste0("  ", sort(u)), sep = "\n")
  cat("\n")
}
sink()
cat("wrote: ", lev_out, "\n", sep = "")

cat("\n=== DONE ===\n")
sessionInfo()
