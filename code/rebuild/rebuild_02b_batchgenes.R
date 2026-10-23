#!/usr/bin/env Rscript
# rebuild_02b_batchgenes.R
# Data-driven screen for genes whose detection rate is driven by dataset of
# origin rather than biology. Runs on the three compartment subsets and writes
# an exclusion list for rebuild_03 to union with the hand-written HVG mask.
#
# Motivation: TRAC is detected in 0.0% of CRA001160 and GSE155698 T/NK cells
# while CD3E is at 82.8% and 75.7% in the same cells. That is structural
# absence from the source reference, not low expression. Symbol harmonization
# (step 1) fixed genes that were MISNAMED across datasets; it cannot fix genes
# that were MISSING. Those genes are perfect dataset markers with zero
# biological content and are eligible for HVG selection.
#
# Two artifact classes are screened:
#   STRUCTURAL - detected in one dataset, exactly zero in another (TRAC-like)
#   GRADED     - detected everywhere but at wildly different rates, the
#                signature of dataset-specific ambient RNA (e.g. COL1A1 at
#                21-28% in T cells from three datasets, 1.9-4.7% in three others)
#
# Output: batchgene_detection_<subset>.tsv  (full genes x datasets table)
#         exclude_batchgenes.rds            (named list, one vector per subset)
#         batchgene_report.txt

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

source("/ix/tbruno/anp795/PDAC_Atlas/code/rebuild/io.R")   # read_obj() / save_obj()

DIR <- Sys.getenv("PDAC_REBUILD_DIR", "/ix/tbruno/anp795/PDAC_Atlas/rebuild_v3")

# A dataset needs enough cells before "zero detection" means anything. At
# n cells and true detection rate p, P(observing zero) = (1-p)^n. At n=150 and
# p=0.05 that is 4.6e-4 -- safe. At n=28 (GSE211644 in stromal) it is 0.24, so
# a quarter of genuinely-expressed genes would read as false structural zeros.
# Datasets below this threshold are excluded from the screen entirely.
MIN_CELLS   <- 150L

PRESENT_CUT <- 0.05   # "detected" in the high dataset
GRADED_CUT  <- 10     # max/min detection ratio for the graded class
GRADED_MIN  <- 0.10   # graded class must reach this in at least one dataset

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

subsets <- c(TNK = "subset_TNK.qs2",
             stromal = "subset_stromal.qs2",
             Bcell = "subset_Bcell.qs2")

exclude <- list()
report  <- character()
add     <- function(...) report <<- c(report, paste0(...))

for (nm in names(subsets)) {
  f <- file.path(DIR, subsets[[nm]])
  msg("---- ", nm, " ----")
  msg("loading ", f)
  s <- read_obj(f)

  d  <- LayerData(s, layer = "data")
  ds <- factor(s$dataset_id)

  tally <- table(ds)
  keep_ds <- names(tally)[tally >= MIN_CELLS]
  drop_ds <- setdiff(names(tally), keep_ds)

  msg("cells per dataset: ", paste(sprintf("%s=%d", names(tally), as.integer(tally)),
                                   collapse = "  "))
  if (length(drop_ds)) {
    msg("excluded from screen (< ", MIN_CELLS, " cells): ",
        paste(drop_ds, collapse = ", "))
  }
  if (length(keep_ds) < 2) {
    msg("fewer than two usable datasets - screen not meaningful, skipping")
    exclude[[nm]] <- character(0)
    rm(s, d); gc(verbose = FALSE)
    next
  }

  # Per-dataset detection rate without materializing per-dataset submatrices:
  # (data > 0) %*% indicator gives genes x datasets counts of detecting cells.
  cells <- colnames(s)[ds %in% keep_ds]
  dk    <- d[, cells, drop = FALSE]
  dsk   <- droplevels(ds[ds %in% keep_ds])

  ind <- Matrix::sparseMatrix(i = seq_along(dsk),
                              j = as.integer(dsk),
                              x = 1,
                              dims = c(length(dsk), nlevels(dsk)),
                              dimnames = list(NULL, levels(dsk)))

  nz  <- as.matrix((dk > 0) %*% ind)
  det <- sweep(nz, 2, Matrix::colSums(ind), "/")

  mx <- apply(det, 1, max)
  mn <- apply(det, 1, min)

  structural <- rownames(det)[mx >= PRESENT_CUT & mn == 0]
  graded     <- rownames(det)[mn > 0 & mx >= GRADED_MIN & (mx / mn) >= GRADED_CUT]

  msg("genes screened: ", format(nrow(det), big.mark = ","))
  msg("STRUCTURAL (>=", PRESENT_CUT * 100, "% in one dataset, exactly 0 in another): ",
      length(structural))
  msg("GRADED (>=", GRADED_CUT, "x detection spread, nonzero everywhere): ",
      length(graded))

  if (length(structural)) {
    top <- head(sort(mx[structural], decreasing = TRUE), 40)
    msg("top structural absences by max detection rate:")
    print(round(det[names(top), , drop = FALSE] * 100, 1))
  }
  if (length(graded)) {
    topg <- head(sort(mx[graded], decreasing = TRUE), 25)
    msg("top graded genes by max detection rate:")
    print(round(det[names(topg), , drop = FALSE] * 100, 1))
  }

  # Full table for review
  out <- data.frame(gene = rownames(det),
                    round(det * 100, 3),
                    max_pct = round(mx * 100, 3),
                    min_pct = round(mn * 100, 3),
                    class = ifelse(rownames(det) %in% structural, "STRUCTURAL",
                             ifelse(rownames(det) %in% graded, "GRADED", "")),
                    check.names = FALSE, stringsAsFactors = FALSE)
  out <- out[order(-out$max_pct), ]
  tf  <- file.path(DIR, paste0("batchgene_detection_", nm, ".tsv"))
  write.table(out, tf, sep = "\t", quote = FALSE, row.names = FALSE)
  msg("wrote ", tf)

  exclude[[nm]] <- union(structural, graded)

  add("== ", nm, " ==")
  add("  datasets screened: ", paste(keep_ds, collapse = ", "))
  if (length(drop_ds)) add("  datasets too small to screen: ", paste(drop_ds, collapse = ", "))
  add("  structural absences: ", length(structural))
  add("  graded/ambient:      ", length(graded))
  add("  total excluded:      ", length(exclude[[nm]]))
  add("")

  rm(s, d, dk, ind, nz, det); gc(verbose = FALSE)
}

saveRDS(exclude, file.path(DIR, "exclude_batchgenes.rds"))
msg("wrote ", file.path(DIR, "exclude_batchgenes.rds"))

writeLines(c(
  "Batch-driven gene screen",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("script:", "rebuild_02b_batchgenes.R"),
  "",
  paste0("MIN_CELLS=", MIN_CELLS, "  PRESENT_CUT=", PRESENT_CUT,
         "  GRADED_CUT=", GRADED_CUT, "  GRADED_MIN=", GRADED_MIN),
  "",
  "STRUCTURAL = detected in >=5% of cells in one dataset and in exactly 0 cells",
  "  in another dataset that has >=150 cells. Means the gene is absent from that",
  "  dataset's source reference. Cannot be fixed by symbol harmonization or by",
  "  Harmony. Must not be used as a marker and must not enter HVG selection.",
  "",
  "GRADED = nonzero in every dataset but with >=10x spread in detection rate.",
  "  Usually dataset-specific ambient RNA. Judgement call - review the tsv",
  "  before excluding, since a real biological difference between cohorts can",
  "  look like this. Do not exclude a gene you intend to use as a marker.",
  "",
  "REVIEW THE TSV BEFORE TRUSTING THE GRADED LIST.",
  "",
  report
), file.path(DIR, "batchgene_report.txt"))

msg("wrote ", file.path(DIR, "batchgene_report.txt"))
msg("done")
