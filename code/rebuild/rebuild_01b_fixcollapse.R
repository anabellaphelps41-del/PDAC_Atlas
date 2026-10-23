#!/usr/bin/env Rscript
# rebuild_01b_fixcollapse.R
#
# Step 1 left 315 features unmerged that should have been collapsed:
#   66  make.unique() ".N" duplicates whose base symbol also exists
#        (PSMA2.1/PSMA2, NDUFB8.1/NDUFB8, TMBIM4.1/TMBIM4, UQCR11.1/UQCR11 ...)
#   250 HGNC prev_symbol pairs where both old and new survived as separate rows
#        (ATP5B/ATP5F1B, SELM/SELENOM, SEPT7/SEPTIN7, DARC/ACKR1 ...)
#
# Evidence: perfectly complementary per-dataset detection rates. ATP5B is
# detected in 53.2/32.2/48.1% of cells in three datasets and exactly 0% in the
# other three; ATP5F1B is exactly the inverse. Same gene, split across two rows,
# behaving as a perfect dataset marker with no biological content.
#
# This operates on the COUNTS layer of atlas_harmonized_counts.rds. Counts sum;
# that is why a full rebuild from the 52 GB original is unnecessary. The old
# log-normalized data layer is discarded and recomputed, because normalized
# values cannot be summed.
#
# In : rebuild/atlas_harmonized_counts.rds
#      rebuild/harmonization_fix_collapse.tsv
# Out: rebuild/atlas_harmonized_counts_v2.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

DIR  <- "/ix/tbruno/anp795/PDAC_Atlas/rebuild"
IN   <- file.path(DIR, "atlas_harmonized_counts.rds")
FIX  <- file.path(DIR, "harmonization_fix_collapse_v3.tsv")
OUT  <- file.path(DIR, "atlas_harmonized_counts_v2.rds")

msg <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep=""); flush.console() }

# --------------------------------------------------------------------------
msg("reading collapse table ", FIX)
fx <- read.delim(FIX, stringsAsFactors = FALSE)
stopifnot(all(c("feature", "collapse_into") %in% colnames(fx)))
fx <- fx[!duplicated(fx$feature), ]
msg("features to merge: ", nrow(fx))
print(table(fx$reason))

# A feature must never be both a source and a target - that would make the
# result depend on evaluation order.
chain <- intersect(fx$feature, fx$collapse_into)
if (length(chain)) {
  stop("chained collapse detected (feature is also a target): ",
       paste(head(chain, 20), collapse = ", "),
       "\nResolve to a single target per group before rerunning.")
}

msg("loading ", IN)
obj <- readRDS(IN)
msg("loaded: ", format(nrow(obj), big.mark=","), " genes x ",
    format(ncol(obj), big.mark=","), " cells")

cts <- LayerData(obj, layer = "counts")
md  <- obj@meta.data
rm(obj); gc(verbose = FALSE)
msg("counts extracted, object released")

# --------------------------------------------------------------------------
# Sanity: every source and target must exist
# --------------------------------------------------------------------------
g <- rownames(cts)
miss_src <- setdiff(fx$feature, g)
miss_tgt <- setdiff(fx$collapse_into, g)
if (length(miss_src)) msg("WARNING: ", length(miss_src), " source features not in matrix, skipping: ",
                          paste(head(miss_src, 10), collapse=", "))
if (length(miss_tgt)) msg("NOTE: ", length(miss_tgt), " targets absent - those become renames")
fx <- fx[fx$feature %in% g, ]
msg("applicable merges: ", nrow(fx))

total_before <- sum(cts@x)
cells_before <- Matrix::colSums(cts)
msg("total counts before: ", format(total_before, big.mark=",", scientific=FALSE))

# --------------------------------------------------------------------------
# Aggregate: build feature -> target map, then sum rows within each group
# --------------------------------------------------------------------------
tgt <- g
i   <- match(fx$feature, g)
tgt[i] <- fx$collapse_into

f   <- factor(tgt, levels = unique(tgt))
msg("features: ", format(length(g), big.mark=","), " -> ",
    format(nlevels(f), big.mark=","), "  (", length(g) - nlevels(f), " merged away)")

ind <- Matrix::sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
                            dims = c(length(f), nlevels(f)),
                            dimnames = list(NULL, levels(f)))

msg("aggregating (t(indicator) %*% counts)")
new_cts <- Matrix::t(ind) %*% cts
new_cts <- as(new_cts, "CsparseMatrix")
rm(ind, cts); gc(verbose = FALSE)

# --------------------------------------------------------------------------
# Verification - the whole point of this script is that nothing is lost
# --------------------------------------------------------------------------
total_after <- sum(new_cts@x)
cells_after <- Matrix::colSums(new_cts)

msg("total counts after : ", format(total_after, big.mark=",", scientific=FALSE))
if (!isTRUE(all.equal(total_before, total_after))) {
  stop("TOTAL COUNTS CHANGED - aggregation is wrong. before=", total_before,
       " after=", total_after)
}
if (!isTRUE(all.equal(unname(cells_before), unname(cells_after)))) {
  stop("PER-CELL TOTALS CHANGED - aggregation is wrong.")
}
msg("verified: total and per-cell UMI counts identical before and after")

if (any(duplicated(rownames(new_cts)))) {
  stop("duplicate rownames after aggregation - ",
       sum(duplicated(rownames(new_cts))), " found")
}
msg("verified: no duplicate feature names")

# Spot-check a few known pairs: the target should now be nonzero in datasets
# where it previously read exactly 0.
spot <- intersect(c("ATP5F1B","SELENOM","SEPTIN7","ACKR1","PSMA2","MICOS10","SNU13"),
                  rownames(new_cts))
if (length(spot) && "dataset_id" %in% colnames(md)) {
  msg("spot check - % of cells detecting, per dataset (expect NO zero columns):")
  sub <- new_cts[spot, , drop = FALSE]
  dsf <- factor(md$dataset_id)
  ind2 <- Matrix::sparseMatrix(i = seq_along(dsf), j = as.integer(dsf), x = 1,
                               dims = c(length(dsf), nlevels(dsf)),
                               dimnames = list(NULL, levels(dsf)))
  det <- as.matrix((sub > 0) %*% ind2)
  det <- sweep(det, 2, Matrix::colSums(ind2), "/") * 100
  print(round(det, 1))
  if (any(det == 0)) {
    msg("*** a spot-check gene is still exactly 0 in some dataset. Either that ",
        "pair was not in the collapse table, or a third symbol for the same ",
        "gene is still separate. Investigate before proceeding. ***")
  }
}

# --------------------------------------------------------------------------
# Rebuild and renormalize
# --------------------------------------------------------------------------
msg("creating Seurat object")
new <- CreateSeuratObject(counts = new_cts, meta.data = md,
                          min.cells = 0, min.features = 0)
rm(new_cts); gc(verbose = FALSE)

msg("NormalizeData(LogNormalize, scale.factor = 1e4)")
new <- NormalizeData(new, normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)

# nCount_RNA is recomputed by CreateSeuratObject; it must match the original
if ("nCount_RNA" %in% colnames(md)) {
  d <- max(abs(new$nCount_RNA - md$nCount_RNA))
  msg("max |nCount_RNA change| vs step-1 metadata: ", d, " (expect 0)")
}

msg("final: ", format(nrow(new), big.mark=","), " genes x ",
    format(ncol(new), big.mark=","), " cells")
msg("writing ", OUT, " (compress=FALSE)")
saveRDS(new, OUT, compress = FALSE)
msg("wrote ", round(file.size(OUT)/2^30, 2), " GB")

writeLines(c(
  "atlas_harmonized_counts_v2.rds",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "script: rebuild_01b_fixcollapse.R",
  paste("source:", IN),
  "",
  "Fixes two harmonization failures left by rebuild_01_harmonize.R:",
  "  1. make.unique() '.N' suffixes never stripped, so e.g. PSMA2.1 stayed a",
  "     separate feature from PSMA2.",
  "  2. HGNC prev_symbol pairs where both old and new symbols survived as",
  "     separate rows (ATP5B/ATP5F1B, SELM/SELENOM, SEPT7/SEPTIN7, DARC/ACKR1).",
  "",
  paste("Features merged:", nrow(fx)),
  paste("Feature count:", "49,570 ->", nrow(new)),
  "",
  "Collapse rule: prev_symbol pairs ONLY, and only where the old symbol is not",
  "itself currently approved for a different gene. All 606 alias_symbol-based",
  "candidates were REJECTED - HGNC alias fields list unrelated genes (AMN/ABCD1,",
  "AR/AKR1B1, MIF/AMH) and collapsing on them would merge distinct genes.",
  "",
  "Verified: total UMI and per-cell UMI identical before and after aggregation;",
  "no duplicate feature names.",
  "",
  "NOT fixed here: ~400-600 residual structural absences per compartment, which",
  "are genuine reference-version differences (clone-ID features, lncRNAs).",
  "Those belong in the HVG exclusion mask, not in a collapse table."
), file.path(DIR, "PROVENANCE_v2.txt"))

msg("wrote ", file.path(DIR, "PROVENANCE_v2.txt"))
msg("done")
sessionInfo()
