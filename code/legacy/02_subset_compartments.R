#!/usr/bin/env Rscript
# 02_subset_compartments.R
# Load the full PDAC atlas ONCE, carve out compartments, strip everything
# that will be recomputed downstream, and write small per-compartment .rds files.
#
# Run this once. Every subclustering run afterward loads a subset, not the atlas.
#
# Usage: Rscript 02_subset_compartments.R <atlas.rds> <outdir> <compartments.tsv>
#
# compartments.tsv is TAB-separated with a header and 3 columns:
#   compartment    annot_col              labels
#   tnk            cell_type_level2       T cell;NK cell;CD8 T;CD4 T
#   caf            cell_type_level2       Fibroblast;CAF;myCAF;iCAF
#
# Fill it in from the output of 01_inspect_atlas.R. Labels must match EXACTLY,
# including case and spaces. The script fails loudly if any label is missing.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 02_subset_compartments.R <atlas.rds> <outdir> <compartments.tsv>")
}
obj_path  <- args[1]
outdir    <- args[2]
cfg_path  <- args[3]

# Keep the raw counts layer only. Normalization, HVGs, scaling, PCA, neighbors
# and UMAP must all be recomputed within each compartment anyway -- carrying the
# atlas-level versions forward is both wrong and the reason the file is 48 GB.
KEEP_ASSAY  <- "RNA"
KEEP_LAYERS <- c("counts")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cfg <- read.delim(cfg_path, stringsAsFactors = FALSE, comment.char = "#")
required <- c("compartment", "annot_col", "labels")
if (!all(required %in% colnames(cfg))) {
  stop("compartments.tsv must have columns: ", paste(required, collapse = ", "))
}
cat("=== CONFIG ===\n"); print(cfg, row.names = FALSE); cat("\n")

cat("=== LOADING ATLAS ===\n")
t0 <- Sys.time()
obj <- readRDS(obj_path)
cat("readRDS: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min\n", sep = "")
cat("cells: ", ncol(obj), " | in-memory: ",
    round(as.numeric(object.size(obj)) / 1e9, 2), " GB\n\n", sep = "")

md <- obj@meta.data
is_v5 <- inherits(obj[[KEEP_ASSAY]], "Assay5")
cat("assay class: ", if (is_v5) "Assay5 (Seurat v5)" else "Assay (Seurat v4)", "\n\n", sep = "")

# --- validate the whole config BEFORE doing any expensive work -----------------
for (i in seq_len(nrow(cfg))) {
  col <- cfg$annot_col[i]
  if (!col %in% colnames(md)) {
    stop("annot_col '", col, "' not found in metadata. Available: ",
         paste(colnames(md), collapse = ", "))
  }
  labs    <- trimws(strsplit(cfg$labels[i], ";", fixed = TRUE)[[1]])
  present <- unique(as.character(md[[col]]))
  missing <- setdiff(labs, present)
  if (length(missing)) {
    stop("compartment '", cfg$compartment[i], "': these labels are not in '", col, "': ",
         paste(missing, collapse = " | "),
         "\n  Did you mean any of: ", paste(head(sort(present), 40), collapse = " | "))
  }
}
cat("config validated: all annotation columns and labels exist.\n\n")

manifest <- list()

for (i in seq_len(nrow(cfg))) {
  cmp  <- cfg$compartment[i]
  col  <- cfg$annot_col[i]
  labs <- trimws(strsplit(cfg$labels[i], ";", fixed = TRUE)[[1]])

  cat("=== COMPARTMENT: ", cmp, " ===\n", sep = "")
  cells <- rownames(md)[as.character(md[[col]]) %in% labs]
  cat("cells selected: ", length(cells), " (",
      round(100 * length(cells) / ncol(obj), 2), "% of atlas)\n", sep = "")
  if (length(cells) == 0) stop("compartment '", cmp, "' selected 0 cells.")
  if (length(cells) < 500) {
    warning("compartment '", cmp, "' has only ", length(cells),
            " cells -- subclustering will be unstable.")
  }

  sub <- subset(obj, cells = cells)

  # Strip to raw counts of one assay; drop reductions, graphs, extra assays.
  sub <- if (is_v5) {
    DietSeurat(sub, assays = KEEP_ASSAY, layers = KEEP_LAYERS,
               dimreducs = NULL, graphs = NULL)
  } else {
    # Seurat v4 objects dislike an empty data slot, so keep it. The real space
    # hog is scale.data (dense, genes x cells) -- that is what we are dropping.
    DietSeurat(sub, assays = KEEP_ASSAY, counts = TRUE, data = TRUE,
               scale.data = FALSE, dimreducs = NULL, graphs = NULL)
  }
  DefaultAssay(sub) <- KEEP_ASSAY

  # v5 objects split per-dataset keep one counts layer per batch. Join them so
  # downstream integration decisions are explicit rather than inherited.
  if (is_v5) {
    lyr <- SeuratObject::Layers(sub[[KEEP_ASSAY]])
    if (length(lyr) > 1) {
      cat("joining ", length(lyr), " counts layers\n", sep = "")
      sub[[KEEP_ASSAY]] <- JoinLayers(sub[[KEEP_ASSAY]])
    }
  }

  # Drop genes detected in zero cells of this compartment.
  # `layer=` is Seurat v5; v4 uses `slot=`. Passing the wrong one is a hard error.
  cts <- if (is_v5) {
    GetAssayData(sub, assay = KEEP_ASSAY, layer = "counts")
  } else {
    GetAssayData(sub, assay = KEEP_ASSAY, slot = "counts")
  }
  # rowSums(cts) > 0 rather than rowSums(cts > 0): counts are non-negative, so
  # this is equivalent and avoids materializing a second sparse matrix.
  keep_g <- Matrix::rowSums(cts) > 0
  cat("genes: ", nrow(sub), " -> ", sum(keep_g), " (dropped all-zero)\n", sep = "")
  rm(cts); gc(verbose = FALSE)
  sub <- subset(sub, features = rownames(sub)[keep_g])
  gc(verbose = FALSE)

  out_rds <- file.path(outdir, paste0(cmp, "_counts.rds"))
  saveRDS(sub, out_rds, compress = TRUE)
  cat("wrote: ", out_rds, " (", round(file.size(out_rds) / 1e9, 2), " GB)\n", sep = "")

  out_md <- file.path(outdir, paste0(cmp, "_metadata.csv.gz"))
  m <- sub@meta.data; m$cell_barcode <- rownames(m)
  write.csv(m, gzfile(out_md), row.names = FALSE)
  cat("wrote: ", out_md, "\n", sep = "")

  manifest[[cmp]] <- data.frame(
    compartment = cmp, annot_col = col,
    labels = paste(labs, collapse = ";"),
    n_cells = ncol(sub), n_genes = nrow(sub),
    rds_gb = round(file.size(out_rds) / 1e9, 3),
    stringsAsFactors = FALSE
  )

  rm(sub); gc(verbose = FALSE)
  cat("\n")
}

man <- do.call(rbind, manifest)
write.csv(man, file.path(outdir, "compartment_manifest.csv"), row.names = FALSE)
cat("=== MANIFEST ===\n"); print(man, row.names = FALSE)
cat("\n=== DONE ===\n")
sessionInfo()
