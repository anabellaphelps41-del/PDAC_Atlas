#!/usr/bin/env Rscript
# =============================================================================
# rebuild_00_provenance_audit.R          (v2 - barcode-prefix aware)
#
# Metadata-only audit of the compartment subsets. Drops the counts with
# rm()+gc() straight after each readRDS, so peak memory is one subset load.
#
# WHY v2: rebuild_02_subset.R requires only celltype_level2, celltype_level3,
# dataset_id and nCount_RNA, and step 1 pruned metadata from 212 to 21 columns.
# There may therefore be NO sample or patient column at all. If so, the cell
# barcode prefix is the only surviving sample identifier, so this script parses
# it. A merged Seurat object normally carries "<sample>_<16bp barcode>-1".
#
# Three questions, from the 2026-08-29 source-methods audit:
#
#   Q1  PBMC samples inside GSE155698?  *** BLOCKS THE T/NK RERUN ***
#       Steele et al. 2020 (doi:10.1038/s43018-020-00121-4) deposited 16 PDA
#       tissue + 3 adjacent normal + 16 patient PBMC + 4 healthy-volunteer PBMC
#       under ONE accession. GSE155698 gives 46,523 T/NK cells; with GSE211644
#       dropped that is 35% of the subset. Blood T cells are naive/CCR7+/SELL+
#       dominant where tumour-infiltrating T cells are exhausted and tissue
#       resident. If blood is in here it is a second cause of the Step 3
#       under-correction, independent of GSE211644, and the rerun should
#       exclude it too.
#
#   Q2  PUMCH (= CRA001160) re-contained inside GSE211644?
#       Schalck et al. 2022 (doi:10.1158/2159-8290.CD-21-1248): "UMI matrices
#       for the MDA2 and PUMCH cohorts were also generated from the provided
#       FASTQ files". Does not block the rerun - dropping GSE211644 removes any
#       duplicates - but it blocks patient counts and the manuscript table.
#
#   Q3  Which GSE154778 supplementary product was ingested?
#       45.1% of its cells have percent.mt == 0, median 1.32%. Lin et al.
#       applied no mitochondrial filter, so their QC did not cause it.
#       Zeros spread EVENLY across samples -> the merged dgeMtx.csv.gz, which
#       is gene-filtered. Zeros CONCENTRATED in some samples -> shallow
#       libraries, and percent.mt is real but unusable in those samples.
#
# Usage: Rscript rebuild_00_provenance_audit.R [subset ...]   default: all three
# =============================================================================

suppressPackageStartupMessages({ library(Seurat); library(SeuratObject) })

DIR     <- "/ix/tbruno/anp795/PDAC_Atlas/rebuild"
args    <- commandArgs(trailingOnly = TRUE)
SUBSETS <- if (length(args)) args else c("TNK", "stromal", "Bcell")

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
rep_lines <- character(); add <- function(...) rep_lines <<- c(rep_lines, paste0(...))

# Broad on purpose. A false positive costs one line of reading; a false
# negative costs a rerun.
PBMC_RX  <- "PBMC|BLOOD|_PB$|^PB[_-]|PERIPH|BUFFY|HEALTHY|HV[0-9]|DONOR|NORMAL_PB"
PUMCH_RX <- "PUMCH|CRA0|PENG|^T[0-9]{1,2}$|^N[0-9]{1,2}$"

# Strip a trailing 10x barcode to leave the sample prefix.
# Handles "SAMPLE_AAACCTGAGAAACCAT-1", "SAMPLE-AAACCTGAGAAACCAT", "AAAC...-1".
barcode_prefix <- function(x) {
  p <- sub("[ACGTN]{14,20}(-[0-9]+)?$", "", x)
  p <- sub("[._-]+$", "", p)
  p[!nzchar(p)] <- "(no prefix)"
  p
}

looks_like_sample <- function(x) {
  if (is.numeric(x)) return(FALSE)
  u <- length(unique(x)); u >= 2 && u <= 500
}

first <- TRUE
for (SUB in SUBSETS) {
  f <- file.path(DIR, sprintf("subset_%s.rds", SUB))
  if (!file.exists(f)) { msg("MISSING, skipping: ", f); next }

  msg("==== ", SUB, " ====")
  msg("loading ", f, "  (", sprintf("%.2f GB", file.size(f) / 1e9), " on disk)")
  s  <- readRDS(f)
  md <- s@meta.data
  cn <- colnames(s)
  rm(s); gc(verbose = FALSE)
  msg("metadata: ", nrow(md), " cells x ", ncol(md), " columns")

  cat("\n--- metadata columns ---\n")
  print(data.frame(column = colnames(md),
                   class  = vapply(md, function(x) class(x)[1], ""),
                   n_uniq = vapply(md, function(x) length(unique(x)), 1L),
                   row.names = NULL))
  if (first) {
    writeLines(c("Metadata columns in the compartment subsets",
                 paste("subset inspected:", SUB), "",
                 paste(sprintf("%2d. %-30s %-10s n_uniq=%d", seq_along(md), colnames(md),
                               vapply(md, function(x) class(x)[1], ""),
                               vapply(md, function(x) length(unique(x)), 1L)),
                       collapse = "\n"), "",
                 "First 10 cell names (barcode prefixes are the fallback sample ID):",
                 paste(head(cn, 10), collapse = "\n")),
               file.path(DIR, "audit_metadata_columns.txt"))
    msg("wrote audit_metadata_columns.txt")
    first <- FALSE
  }

  stopifnot("dataset_id" %in% colnames(md))
  ds <- as.character(md$dataset_id)
  cat("\n--- cells per dataset_id ---\n"); print(table(ds))
  add("== ", SUB, " ==  cells: ", nrow(md))

  # ---- build the candidate sample columns, INCLUDING the barcode prefix ----
  cand <- setdiff(names(md)[vapply(md, looks_like_sample, TRUE)], "dataset_id")
  sample_cols <- list()
  for (cl in cand) sample_cols[[cl]] <- as.character(md[[cl]])

  bp <- barcode_prefix(cn)
  cat("\n--- barcode prefix parse ---\n")
  cat("  example cell names : ", paste(head(cn, 3), collapse = " | "), "\n")
  cat("  parsed prefixes    : ", paste(head(bp, 3), collapse = " | "), "\n")
  cat("  unique prefixes    : ", length(unique(bp)), "\n")
  if (length(unique(bp)) >= 2 && length(unique(bp)) <= 500) {
    sample_cols[["barcode_prefix"]] <- bp
    msg("barcode prefix looks usable as a sample identifier (",
        length(unique(bp)), " values)")
  } else {
    msg("barcode prefix NOT usable (", length(unique(bp)),
        " unique values). If no metadata column carries sample identity either,",
        " Q1 and Q2 cannot be answered from the subsets and you must go back to",
        " combined_9_datasets_hierarchical_annotation_seurat.rds (212 columns).")
    add("  WARNING ", SUB, ": no usable sample identifier found")
  }

  if (!length(sample_cols)) { add(""); rm(md); gc(verbose = FALSE); next }
  msg("sample-level candidates: ", paste(names(sample_cols), collapse = ", "))

  for (cl in names(sample_cols)) {
    v  <- sample_cols[[cl]]
    ct <- table(ds, v)
    if (max(apply(ct > 0, 1, sum)) < 2) next   # one value per dataset = dataset_id again

    write.table(as.data.frame.matrix(ct),
                file.path(DIR, sprintf("audit_%s_dataset_by_%s.tsv", SUB, make.names(cl))),
                sep = "\t", quote = FALSE, col.names = NA)

    cat(sprintf("\n--- '%s': values within each dataset (%d unique) ---\n", cl, ncol(ct)))
    for (d in rownames(ct)) {
      vals <- colnames(ct)[ct[d, ] > 0]
      cat(sprintf("  %-12s n=%-4d : %s\n", d, length(vals),
                  paste(head(sort(vals), 30), collapse = ", ")))
      if (length(vals) > 30) cat("                  ... and ", length(vals) - 30, " more\n")
    }

    # ---------------- Q1 -----------------------------------------------------
    hit <- grep(PBMC_RX, colnames(ct), ignore.case = TRUE, value = TRUE)
    if (length(hit)) {
      sc <- ct[, hit, drop = FALSE]; sc <- sc[rowSums(sc) > 0, , drop = FALSE]
      cat("\n  *** Q1  BLOOD-LIKE SAMPLE NAMES via '", cl, "' ***\n", sep = "")
      print(sc)
      for (d in rownames(sc)) {
        n <- sum(sc[d, ]); tot <- sum(ct[d, ])
        cat(sprintf("      %s: %s of %s cells (%.1f%%) in blood-like samples\n",
                    d, format(n, big.mark = ","), format(tot, big.mark = ","), 100 * n / tot))
        add("  Q1 ", SUB, "/", cl, ": ", d, " ", n, "/", tot,
            sprintf(" (%.1f%%) blood-like -> ", 100 * n / tot), paste(hit, collapse = ", "))
      }
      cat("      DECISION RULE: if this is >5% of GSE155698 in T/NK, exclude those\n")
      cat("      samples from the Step 3 rerun alongside GSE211644.\n")
    } else {
      cat(sprintf("\n  Q1 via '%s': no blood-like sample names matched.\n", cl))
    }

    # ---------------- Q2 -----------------------------------------------------
    if (all(c("GSE211644", "CRA001160") %in% rownames(ct))) {
      v211 <- colnames(ct)[ct["GSE211644", ] > 0]
      cra  <- colnames(ct)[ct["CRA001160", ] > 0]
      shared <- intersect(v211, cra)
      pum    <- grep(PUMCH_RX, v211, ignore.case = TRUE, value = TRUE)
      cat("\n  *** Q2  GSE211644 sample names via '", cl, "' ***\n", sep = "")
      cat("      values in GSE211644          : ", length(v211), "\n")
      cat("      shared verbatim with CRA001160: ",
          if (length(shared)) paste(shared, collapse = ", ") else "(none)", "\n")
      cat("      matching PUMCH-like pattern   : ",
          if (length(pum)) paste(pum, collapse = ", ") else "(none)", "\n")
      add("  Q2 ", SUB, "/", cl, ": ", length(shared), " shared with CRA001160, ",
          length(pum), " PUMCH-like")
    }
  }

  # ---------------- Q3 -------------------------------------------------------
  if ("percent.mt" %in% colnames(md)) {
    agg <- do.call(rbind, lapply(split(md$percent.mt, ds), function(x)
      data.frame(n = length(x), pct_zero = round(100 * mean(x == 0, na.rm = TRUE), 2),
                 median = round(median(x, na.rm = TRUE), 3))))
    agg$dataset_id <- rownames(agg)
    cat("\n--- Q3  percent.mt by dataset ---\n")
    print(agg[, c("dataset_id", "n", "pct_zero", "median")])

    for (cl in names(sample_cols)) {
      v <- sample_cols[[cl]]; k <- ds == "GSE154778"
      if (sum(k) == 0 || length(unique(v[k])) < 2) next
      per <- do.call(rbind, lapply(split(md$percent.mt[k], v[k]), function(x)
        data.frame(n = length(x), pct_zero = round(100 * mean(x == 0, na.rm = TRUE), 2),
                   median = round(median(x, na.rm = TRUE), 3))))
      per$sample <- rownames(per)
      cat(sprintf("\n--- Q3  GSE154778 percent.mt by '%s' ---\n", cl))
      print(per[, c("sample", "n", "pct_zero", "median")])
      cat("    READ: zeros EVEN across samples -> merged GSE154778_dgeMtx.csv.gz\n")
      cat("          was ingested and is gene-filtered; percent.mt unrecoverable.\n")
      cat("          Zeros CONCENTRATED -> shallow libraries; percent.mt is real\n")
      cat("          but unusable in those samples only.\n")
      write.table(per[, c("sample", "n", "pct_zero", "median")],
                  file.path(DIR, "audit_percent_mt_by_sample.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      add("  Q3 GSE154778 zero-fraction ",
          sprintf("%.1f%%-%.1f%% across %d samples (spread %.1f pts)",
                  min(per$pct_zero), max(per$pct_zero), nrow(per),
                  diff(range(per$pct_zero))))
      break
    }
  }

  add(""); rm(md); gc(verbose = FALSE)
}

writeLines(c("Provenance audit - metadata only",
             paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
             "script: rebuild_00_provenance_audit.R (v2, barcode-prefix aware)", "",
             "Q1 PBMC in GSE155698  -- BLOCKS the T/NK rerun",
             "Q2 PUMCH/CRA001160 duplication in GSE211644 -- blocks patient-level claims",
             "Q3 which GSE154778 supplementary product was ingested", "",
             rep_lines),
           file.path(DIR, "audit_report.txt"))
msg("wrote ", file.path(DIR, "audit_report.txt")); msg("done")
