#!/usr/bin/env Rscript
# rebuild_02_subset.R
# Step 2 of the PDAC atlas rebuild: split the harmonized object into the three
# compartments that get their own pipeline (T/NK, stromal/CAF, B/plasma).
#
# Input : /ix/tbruno/anp795/PDAC_Atlas/rebuild/atlas_harmonized_counts_v2.rds
# Output: subset_TNK.rds, subset_stromal.rds, subset_Bcell.rds  (+ PROVENANCE.txt)
#
# Notes:
#   - Subsetting uses celltype_level2, NOT level1. level1/level2 disagree on
#     doublet status for 31,521 cells; level2 is the deliberate choice.
#   - "Mixed/possible doublet" (30,140) is its own level2 value and is therefore
#     excluded from all three subsets by construction, not by a filter.
#   - The data layer is NOT recomputed. LogNormalize(scale.factor=1e4) is a
#     per-cell operation, so values are identical whether computed on the full
#     atlas or within a subset. Re-running it would be a no-op.
#   - percent.mt IS recomputed here, purely as a diagnostic, to establish which
#     datasets structurally lack MT genes.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

IN  <- "/ix/tbruno/anp795/PDAC_Atlas/rebuild/atlas_harmonized_counts_v2.rds"
OUT <- "/ix/tbruno/anp795/PDAC_Atlas/rebuild"

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

# --------------------------------------------------------------------------
# Load and validate
# --------------------------------------------------------------------------
msg("loading ", IN)
obj <- readRDS(IN)
msg("loaded: ", format(nrow(obj), big.mark = ","), " genes x ",
    format(ncol(obj), big.mark = ","), " cells")

md <- obj@meta.data

# Metadata was pruned 212 -> 21 columns in step 1. Fail loudly if the columns
# this script depends on did not survive, rather than silently doing the wrong
# thing (e.g. keeping Neutrophil-like cells because the column is gone).
need <- c("celltype_level2", "celltype_level3", "dataset_id", "nCount_RNA")
missing <- setdiff(need, colnames(md))
if (length(missing)) {
  stop("required metadata column(s) absent after step-1 pruning: ",
       paste(missing, collapse = ", "),
       "\nStop and recover them before subsetting.")
}

msg("celltype_level2 tally:")
print(sort(table(md$celltype_level2), decreasing = TRUE))

# --------------------------------------------------------------------------
# Diagnostic: is percent.mt salvageable at all?
# --------------------------------------------------------------------------
mt_genes <- grep("^MT-", rownames(obj), value = TRUE)
msg("MT- genes present in harmonized matrix: ", length(mt_genes))

if (length(mt_genes) > 0) {
  mt_counts <- Matrix::colSums(LayerData(obj, layer = "counts")[mt_genes, , drop = FALSE])
  pmt <- 100 * mt_counts / md$nCount_RNA
  pmt[!is.finite(pmt)] <- NA_real_
  obj$percent.mt.recomputed <- pmt

  mt_report <- do.call(rbind, lapply(split(pmt, md$dataset_id), function(x) {
    data.frame(n           = length(x),
               pct_zero    = round(100 * mean(x == 0, na.rm = TRUE), 1),
               median_pmt  = round(median(x, na.rm = TRUE), 2))
  }))
  msg("recomputed percent.mt, per dataset:")
  print(mt_report)
  msg("READ THIS: any dataset at ~100% zeros has no MT genes in its source ",
      "matrix. percent.mt cannot be rescued for those datasets by recomputation ",
      "- it is permanently unusable there, and that gets written down, not fixed.")
} else {
  obj$percent.mt.recomputed <- NA_real_
  msg("no MT- genes at all. percent.mt is dead atlas-wide. Record and move on.")
}

# --------------------------------------------------------------------------
# Global exclusion: Neutrophil-like (single-dataset artifact, 85.6% GSE155698)
# %in% rather than != so that NA labels are retained, not silently dropped.
# --------------------------------------------------------------------------
neut <- md$celltype_level3 %in% "Neutrophil-like"
msg("Neutrophil-like cells flagged for exclusion: ", format(sum(neut), big.mark = ","))
if (sum(neut) > 0) {
  msg("their celltype_level2 distribution (expect these to sit under Myeloid):")
  print(table(md$celltype_level2[neut]))
}
keep_global <- !neut

# --------------------------------------------------------------------------
# Subset definitions and expected counts (from the step-1 compartment tally,
# before Neutrophil-like removal). A shortfall means Neutrophil-like cells were
# living inside that compartment - which is itself a finding worth chasing.
# --------------------------------------------------------------------------
defs <- list(
  TNK     = list(level2 = c("T/NK ambiguous", "T/NK"), expected = 193408L),
  stromal = list(level2 = "Stromal",                   expected =  41625L),
  Bcell   = list(level2 = "B/plasma",                  expected =  17975L)
)

written <- list()

for (nm in names(defs)) {
  d      <- defs[[nm]]
  in_cmp <- md$celltype_level2 %in% d$level2
  cells  <- rownames(md)[in_cmp & keep_global]

  msg("---- ", nm, " ----")
  msg("compartment total: ", format(sum(in_cmp), big.mark = ","),
      " | expected: ",       format(d$expected,  big.mark = ","),
      " | after Neutrophil-like exclusion: ", format(length(cells), big.mark = ","))

  if (sum(in_cmp) != d$expected) {
    msg("*** compartment total does not match the expected count. Do not ",
        "proceed past this until you know why. ***")
  }
  if (length(cells) == 0) stop("empty subset for ", nm, " - check level2 spelling")

  s <- subset(obj, cells = cells)
  DefaultAssay(s) <- "RNA"

  msg("cells per dataset_id:")
  print(table(s$dataset_id))
  msg("NOTE: any dataset contributing <100 cells here is a poor Harmony batch ",
      "in step 3. Decide whether to merge or drop it before running Harmony.")

  msg("celltype_level3 composition:")
  print(sort(table(s$celltype_level3), decreasing = TRUE))

  f <- file.path(OUT, paste0("subset_", nm, ".rds"))
  msg("writing ", f, " (compress=FALSE)")
  saveRDS(s, f, compress = FALSE)
  msg("wrote ", round(file.size(f) / 2^30, 2), " GB")

  written[[nm]] <- data.frame(
    subset    = nm,
    criterion = paste0('celltype_level2 %in% c("', paste(d$level2, collapse = '", "'), '")'),
    cells     = length(cells),
    genes     = nrow(s),
    file      = basename(f),
    stringsAsFactors = FALSE
  )

  rm(s); gc(verbose = FALSE)
}

# --------------------------------------------------------------------------
# Provenance file - answers "which annotation did this subset come from?"
# --------------------------------------------------------------------------
prov <- file.path(OUT, "PROVENANCE.txt")
tbl  <- do.call(rbind, written)

writeLines(c(
  "PDAC Atlas rebuild - compartment subsets",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("by:", "rebuild_02_subset.R"),
  paste("source:", IN),
  "",
  "ANNOTATION COLUMN USED FOR SUBSETTING: celltype_level2",
  "",
  "Why level2 and not level1: level1 and level2 disagree on doublet status for",
  "31,521 cells (level1 calls 15,576 doublets, level2 calls 47,097; the extra",
  "came out of Immune 23,381 and Stromal/vascular 8,140). level2 is the",
  "deliberate choice.",
  "",
  "The object carries four competing annotation schemes:",
  "  celltype_level1/2/3, celltype_level3_auto   (hierarchical pass - THIS ONE)",
  "  celltype_broad_reviewed, celltype_fine_reviewed  (earlier reviewed pass)",
  "  celltype_broad, celltype_fine                    (earlier auto pass)",
  "  Level.1-3.Annotation, broad_celltypes, celltypes, new_celltypes (inherited)",
  "Only compartment-level calls (level1/level2) survived review. celltype_level3,",
  "celltype_fine, the 49 clusters and the UMAP did NOT and are not used here.",
  "",
  "GLOBAL EXCLUSION: celltype_level3 == 'Neutrophil-like'",
  "  Single-dataset artifact: 85.6% from GSE155698 alone (22.5% of that whole",
  "  dataset), absent from five of nine datasets. Cannot be used cross-dataset;",
  "  ambient-driven expression would manufacture false ligand-receptor hits.",
  "",
  "NOT INCLUDED IN ANY SUBSET: 'Mixed/possible doublet' (30,140 cells).",
  "  These are their own celltype_level2 value, so they are excluded by",
  "  construction. Their doublet status is UNRESOLVED - see open items.",
  "",
  "Normalization: the data layer is carried over from step 1 unchanged.",
  "  LogNormalize(scale.factor=1e4) is per-cell, so subsetting does not",
  "  invalidate it. Do not re-normalize.",
  "",
  "Subsets written:",
  paste(capture.output(print(tbl, row.names = FALSE)), collapse = "\n"),
  "",
  "Everything saved with compress=FALSE (larger on disk, loads in minutes)."
), prov)

msg("wrote ", prov)
msg("done")
sessionInfo()
