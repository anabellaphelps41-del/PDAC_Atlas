#!/usr/bin/env Rscript
# =============================================================================
# rebuild_01d_decontx.R   --   ambient-RNA removal (DecontX) + per-cell QC floor
#
# Input : <UP>/atlas_harmonized_counts_v3.qs2   (rebuild_01_harmonize.R)
#         <UP>/sample_manifest.tsv               (rebuild_00m_manifest.R)
#         <UP>/doublets_scDblFinder_v3.tsv       (rebuild_01c_doublets.R)
#         <DIR>/cell_routing.tsv                 (rebuild_02a_routing.R, for z)
# Output: <DIR>/decontx_counts.qs2     genes x cells, DecontX counts rounded to integers
#         <DIR>/decontx_percell.tsv    contamination, raw QC metrics, QC floor verdict
#         <DIR>/decontx_by_sample.tsv
#         <DIR>/PROVENANCE_01d_decontx.txt
# (UP = rebuild_v3, DIR = rebuild_v4; see io.R)
#
# WHY ATLAS LEVEL, PER SAMPLE
#   DecontX models each cell as its own cell type's expression profile plus
#   contamination drawn from the OTHER cell types in the same library. Acinar,
#   ductal or fibroblast reads in a B cell can only be recognised as
#   contamination if acinar, ductal and fibroblast cells are in the same run.
#   So this runs before the compartment split, one sample at a time, with every
#   compartment present. Running it on the subsets would not work.
#
# CELLS IN SCOPE
#   samples with include_in_atlas == TRUE, minus scDblFinder doublets, minus
#   celltype_level2 'Mixed/possible doublet' and 'Unassigned' (none of these
#   enter any subset).
#
# CELL-TYPE LABELS (z)
#   DecontX treats each label's average profile as "native". A label that
#   lumps very different cells breaks this: with level2 'B/plasma' as one
#   label, the profile is plasma-like (immunoglobulin-heavy) and ordinary B
#   cells came out ~98% "contamination" in the GSE205013:P03 test, while
#   plasma cells were 3%. So labels are lineage-level but homogeneous:
#     1. cells of the v3 TNK/stromal/Bcell subsets: lineage_reviewed from the
#        reviewed res-0.6 calls (Tcell, NK, Bcell, plasma, fibroblast, mural,
#        schwann, ...), via cell_routing.tsv
#     2. other cells: celltype_level3 (macrophage, monocyte, neutrophil, mast,
#        endothelial, ductal, acinar, endocrine)
#     3. anything left: celltype_level2
#   Labels are lineages, never clusters: a cluster defined by its ambient RNA
#   (e.g. B cells carrying fibroblast reads) must share a label with the clean
#   cells of its lineage, or its contamination would count as native.
#
# NO BACKGROUND PROFILE
#   No raw (unfiltered) droplet matrices are available for these datasets, so
#   DecontX estimates contamination from the cells alone. GSE205013 and
#   GSE212966 were already SoupX-corrected by their authors; they are run too,
#   so one method applies to every dataset.
#
# ROUNDING
#   DecontX counts are fractional. They are rounded to integers (and the zeros
#   that creates are dropped) so count-based methods can use them.
#
# QC FLOOR
#   Per dataset x compartment (level2, T/NK merged): a cell fails if it is more than
#   3 MADs below the median on log10 nFeature or log10 nCount (raw counts).
#   A per-dataset floor removes almost nothing, because the authors already
#   filtered at 100-555 genes. Compartments differ in complexity (neutrophils
#   and T cells vs fibroblasts), so one floor per dataset would cut the
#   low-complexity types. percent.mt is not used (GSE154778's MT column is
#   broken). Groups below 50 cells use the dataset-wide floor. The floor is a
#   lower bound only; doublets are handled by 01c.
#
# TEST RUNS
#   PDAC_DECONTX_SAMPLES="ds:sample;ds:sample" runs only those samples and
#   writes to <DIR>/test_01d/. QC floors are still computed on all cells.
#
# RESUME
#   Each sample is written to <DIR>/tmp_01d/<sample>.qs2 as it finishes, and a
#   rerun skips samples already there. The folder is removed after assembly.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(celda)
  library(BiocParallel)
})

source("/ix/tbruno/anp795/PDAC_Atlas/code/rebuild/io.R")   # read_obj(), save_obj(), dirs

UP   <- UP_DIR
IN   <- file.path(UP, "atlas_harmonized_counts_v3.qs2")
MAN  <- file.path(UP, "sample_manifest.tsv")
DBL  <- file.path(UP, "doublets_scDblFinder_v3.tsv")

# ';' or ',' separated. Use ';' with sbatch --export, which splits on commas.
ONLY <- strsplit(Sys.getenv("PDAC_DECONTX_SAMPLES"), "[;,[:space:]]+")[[1]]
ONLY <- ONLY[nzchar(ONLY)]
TEST <- length(ONLY) > 0
DIR  <- if (TEST) file.path(V4_DIR, "test_01d") else V4_DIR
TMP  <- file.path(DIR, "tmp_01d")
dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
OUTM <- file.path(DIR, "decontx_counts.qs2")
OUTC <- file.path(DIR, "decontx_percell.tsv")
OUTS <- file.path(DIR, "decontx_by_sample.tsv")
PROV <- file.path(DIR, "PROVENANCE_01d_decontx.txt")

MIN_CELLS  <- 100     # samples below this keep raw counts
MIN_TYPES  <- 2       # DecontX needs at least two cell types (with >= 2 cells) in a sample
QC_NMAD    <- 3
QC_MIN_GRP <- 50
SEED       <- 12345   # DecontX default
NTHREADS   <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))

options(width = 200)
msg <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = ""); flush.console() }
t0 <- Sys.time()

# --------------------------------------------------------------------------
# Load, attach manifest and doublet calls, define scope
# --------------------------------------------------------------------------
msg("loading ", IN)
obj <- read_obj(IN)
md  <- obj@meta.data
cts <- LayerData(obj, assay = "RNA", layer = "counts")
rm(obj); gc(verbose = FALSE)
msg("atlas: ", format(nrow(cts), big.mark = ","), " genes x ", format(ncol(cts), big.mark = ","), " cells")

md$sample_key <- paste(md$dataset_id, md$sample_id, sep = ":")
man <- fread(MAN)
mi  <- match(md$sample_key, man$sample_key)
if (anyNA(mi)) stop("cells whose sample is not in the manifest: ", sum(is.na(mi)))
inc <- as.logical(man$include_in_atlas[mi])

dbl <- fread(DBL)
di  <- match(rownames(md), dbl$cell)
if (anyNA(di)) stop("cells missing from the doublet table: ", sum(is.na(di)))
dcl <- dbl$scDblFinder_class[di]

l2 <- as.character(md$celltype_level2)
scope <- inc & dcl != "doublet" & !is.na(l2) & !(l2 %in% c("Mixed/possible doublet", "Unassigned"))
msg("in scope: ", format(sum(scope), big.mark = ","), " cells (included samples, singlets, ",
    "not Mixed/possible doublet or Unassigned)")

ROUT <- file.path(V4_DIR, "cell_routing.tsv")
if (!file.exists(ROUT)) stop("run rebuild_02a_routing.R first: ", ROUT)
rt <- fread(ROUT, select = c("cell", "lineage"))
L3_MAP <- c("Macrophage" = "macrophage", "Monocyte/inflammatory myeloid" = "monocyte",
            "Neutrophil-like" = "neutrophil", "Mast cell" = "mast", "Endothelial" = "endothelial",
            "Ductal/epithelial" = "ductal", "Acinar" = "acinar", "Endocrine" = "endocrine")
Z_MAP  <- c("acinar_epithelial" = "acinar")      # same lineage, one label
z_src  <- rep("level2", nrow(md))
z      <- l2
l3     <- unname(L3_MAP[as.character(md$celltype_level3)])
z[!is.na(l3)] <- l3[!is.na(l3)];  z_src[!is.na(l3)] <- "level3"
ri     <- match(rownames(md), rt$cell)
z[!is.na(ri)] <- rt$lineage[ri[!is.na(ri)]];  z_src[!is.na(ri)] <- "reviewed lineage"
z[z %in% names(Z_MAP)] <- Z_MAP[z[z %in% names(Z_MAP)]]
msg("z labels in scope (source):")
print(table(z = z[scope], source = z_src[scope]))

cs   <- cts[, scope]
rm(cts); gc(verbose = FALSE)
pc   <- data.table(cell = colnames(cs), dataset_id = md$dataset_id[scope],
                   sample_key = md$sample_key[scope], z_label = z[scope], z_source = z_src[scope],
                   compartment = l2[scope])
rm(md); gc(verbose = FALSE)

# --------------------------------------------------------------------------
# Raw QC metrics and the per dataset x compartment floor (all scope cells)
# --------------------------------------------------------------------------
pc[, nCount_raw   := Matrix::colSums(cs)]
pc[, nFeature_raw := Matrix::colSums(cs > 0)]
pc[, `:=`(lf = log10(pmax(nFeature_raw, 1)), lc = log10(pmax(nCount_raw, 1)))]

floor_of <- function(v) {          # median - k * MAD; no floor if MAD is 0
  m <- stats::mad(v)
  if (!is.finite(m) || m == 0) -Inf else stats::median(v) - QC_NMAD * m
}
ds_fl <- pc[, .(f_ds = floor_of(lf), c_ds = floor_of(lc)), by = dataset_id]
pc[, qc_group := fifelse(compartment %in% c("T/NK ambiguous", "T/NK"), "T/NK", compartment)]
gr_fl <- pc[, .(n_grp = .N, f_grp = floor_of(lf), c_grp = floor_of(lc)), by = .(dataset_id, qc_group)]
gr_fl <- merge(gr_fl, ds_fl, by = "dataset_id")
gr_fl[, `:=`(qc_floor_log10_nFeature = ifelse(n_grp >= QC_MIN_GRP, f_grp, f_ds),
             qc_floor_log10_nCount   = ifelse(n_grp >= QC_MIN_GRP, c_grp, c_ds),
             floor_source = ifelse(n_grp >= QC_MIN_GRP, "dataset x compartment", "dataset (group < 50)"))]
pc <- merge(pc, gr_fl[, .(dataset_id, qc_group, qc_floor_log10_nFeature, qc_floor_log10_nCount)],
            by = c("dataset_id", "qc_group"), sort = FALSE)
pc[, qc_pass := lf >= qc_floor_log10_nFeature & lc >= qc_floor_log10_nCount]
pc[, c("lf", "lc") := NULL]
setkey(pc, NULL)
pc <- pc[match(colnames(cs), cell)]          # back to matrix column order
stopifnot(identical(pc$cell, colnames(cs)))

qc_tab <- pc[, .(cells = .N, fail = sum(!qc_pass), pct_fail = round(100 * mean(!qc_pass), 2)),
             by = .(dataset_id, qc_group)][order(dataset_id, qc_group)]
msg("QC floor: ", format(sum(!pc$qc_pass), big.mark = ","), " of ",
    format(nrow(pc), big.mark = ","), " cells fail (",
    sprintf("%.2f%%", 100 * mean(!pc$qc_pass)), ")")
print(dcast(qc_tab, qc_group ~ dataset_id, value.var = "pct_fail"))

# --------------------------------------------------------------------------
# DecontX per sample
# --------------------------------------------------------------------------
samples <- unique(pc$sample_key)
if (TEST) {
  miss <- setdiff(ONLY, samples)
  if (length(miss)) stop("test samples not in scope: ", paste(miss, collapse = ", "))
  samples <- ONLY
  msg("*** TEST RUN: ", length(samples), " sample(s); output in ", DIR, " ***")
}
cols_by <- split(seq_len(nrow(pc)), pc$sample_key)[samples]
safe    <- function(k) gsub("[^A-Za-z0-9._-]", "_", k)
tmp_of  <- function(k) file.path(TMP, paste0(safe(k), ".qs2"))

one <- function(k) {
  f <- tmp_of(k)
  if (file.exists(f)) return(data.table(sample_key = k, status = "cached", seconds = NA_real_))
  ii <- cols_by[[k]]
  zz <- pc$z_label[ii]
  tz <- table(zz)
  t1 <- Sys.time()
  status <- if (length(ii) < MIN_CELLS) "not_run_small_sample"
            else if (sum(tz >= 2) < MIN_TYPES) "not_run_one_cell_type" else "ok"
  out <- list(status = status, cells = pc$cell[ii],
              contamination = rep(NA_real_, length(ii)), counts = NULL, delta = NA_real_)
  if (status == "ok") {
    res <- tryCatch(decontX(cs[, ii], z = zz, seed = SEED, verbose = FALSE),
                    error = function(e) e)
    if (inherits(res, "error")) {
      out$status <- paste0("error: ", conditionMessage(res))
    } else {
      m <- res$decontXcounts
      if (!inherits(m, "dgCMatrix")) m <- as(m, "CsparseMatrix")
      m@x <- round(m@x)
      out$counts        <- drop0(m)
      out$contamination <- as.numeric(res$contamination)
      out$delta         <- paste(round(res$estimates[[1]]$delta, 3), collapse = "/")
    }
  }
  qs2::qs_save(out, f, nthreads = 1)
  data.table(sample_key = k, status = out$status,
             seconds = as.numeric(difftime(Sys.time(), t1, units = "secs")))
}

n_cached <- sum(file.exists(vapply(samples, tmp_of, "")))
msg("DecontX on ", length(samples), " samples (", n_cached, " already done) with ",
    NTHREADS, " workers")
runlog <- rbindlist(bplapply(samples, one,
                             BPPARAM = MulticoreParam(workers = min(NTHREADS, length(samples)))))
print(runlog[order(-seconds)][1:min(.N, 10)])

# --------------------------------------------------------------------------
# Assemble: DecontX counts where it ran, raw counts otherwise
# --------------------------------------------------------------------------
msg("assembling")
pc[, `:=`(decontX_status = NA_character_, decontX_contamination = NA_real_)]
parts <- vector("list", length(samples)); names(parts) <- samples
deltas <- character(length(samples)); names(deltas) <- samples
for (k in samples) {
  r  <- qs2::qs_read(tmp_of(k), nthreads = 1)
  ii <- cols_by[[k]]
  stopifnot(identical(r$cells, pc$cell[ii]))
  pc[ii, `:=`(decontX_status = r$status, decontX_contamination = r$contamination)]
  parts[[k]] <- if (is.null(r$counts)) cs[, ii] else r$counts
  deltas[k]  <- as.character(r$delta)
}
bad <- unique(pc$decontX_status[!pc$decontX_status %in% c("ok", NA)])
if (length(bad)) msg("*** samples NOT decontaminated (raw counts kept): ",
                     paste(unique(pc[decontX_status %in% bad, sample_key]), collapse = ", "),
                     " | ", paste(bad, collapse = "; "))

dx <- do.call(cbind, unname(parts))
rm(parts); gc(verbose = FALSE)
keep_cells <- pc$cell[unlist(cols_by, use.names = FALSE)]
dx <- dx[, keep_cells]
stopifnot(identical(rownames(dx), rownames(cs)))
pc <- pc[match(keep_cells, cell)]
pc[, nCount_decontX := Matrix::colSums(dx)]
pc[, pct_counts_removed := round(100 * (1 - nCount_decontX / nCount_raw), 2)]

msg("writing ", OUTM)
save_obj(dx, OUTM)
msg("wrote ", round(file.size(OUTM) / 2^30, 2), " GB")
fwrite(pc, OUTC, sep = "\t")
msg("wrote ", OUTC)

by_s <- pc[, .(n_cells = .N, status = paste(unique(decontX_status), collapse = ","),
               median_contamination = round(median(decontX_contamination, na.rm = TRUE), 4),
               p90_contamination = round(quantile(decontX_contamination, 0.9, na.rm = TRUE), 4),
               pct_counts_removed = round(100 * (1 - sum(nCount_decontX) / sum(nCount_raw)), 2),
               pct_qc_fail = round(100 * mean(!qc_pass), 2)),
           by = .(dataset_id, sample_key)][order(dataset_id, sample_key)]
by_s[, delta := deltas[sample_key]]
fwrite(by_s, OUTS, sep = "\t")
msg("wrote ", OUTS)

# --------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------
by_ds <- pc[, .(cells = .N, median_contamination = round(median(decontX_contamination, na.rm = TRUE), 4),
                pct_counts_removed = round(100 * (1 - sum(nCount_decontX) / sum(nCount_raw)), 2),
                qc_fail = sum(!qc_pass)), by = dataset_id][order(dataset_id)]
by_z  <- dcast(pc[, .(m = round(median(decontX_contamination, na.rm = TRUE), 3)),
                  by = .(z_label, dataset_id)], z_label ~ dataset_id, value.var = "m")
fl_tab <- gr_fl[, .(dataset_id, qc_group, n_grp, floor_source,
                    min_nFeature = round(10^qc_floor_log10_nFeature),
                    min_nCount   = round(10^qc_floor_log10_nCount))][order(dataset_id, qc_group)]
print(by_ds); print(by_z)

show <- function(x) capture.output(print(x, row.names = FALSE))
writeLines(c(
  "PROVENANCE - ambient RNA removal (DecontX) and per-cell QC floor",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "script   : rebuild_01d_decontx.R",
  paste("input    :", IN),
  paste("output   :", OUTM),
  if (TEST) paste("TEST RUN : samples", paste(samples, collapse = ", ")) else NULL,
  sprintf("method   : celda::decontX %s per sample_key, z = reviewed lineage / level3 / level2 (see header), no background, seed %d",
          as.character(packageVersion("celda")), SEED),
  sprintf("cells    : %s in scope (included samples, singlets, not Mixed/possible doublet or Unassigned)",
          format(nrow(pc), big.mark = ",")),
  sprintf("runtime  : %.1f min with %d workers", as.numeric(difftime(Sys.time(), t0, units = "mins")), NTHREADS),
  "",
  "BY DATASET", show(by_ds), "",
  "MEDIAN CONTAMINATION BY CELL TYPE (z) x DATASET", show(by_z), "",
  sprintf("QC FLOOR: log10 nFeature and log10 nCount >= median - %d MAD per dataset x compartment", QC_NMAD),
  sprintf("  %s of %s cells fail (%.2f%%). Floors (back-transformed):",
          format(sum(!pc$qc_pass), big.mark = ","), format(nrow(pc), big.mark = ","),
          100 * mean(!pc$qc_pass)),
  show(fl_tab), "",
  "PERCENT FAILING, compartment x dataset", show(dcast(qc_tab, qc_group ~ dataset_id, value.var = "pct_fail")), "",
  "DecontX LABELS (z) BY SOURCE", show(as.data.frame.matrix(table(pc$z_label, pc$z_source))), "",
  "NOTES",
  "  - counts are rounded DecontX estimates; raw counts stay in the RNA assay downstream.",
  "  - samples where DecontX did not run keep raw counts (decontX_status says why).",
  "  - per-cell table: decontx_percell.tsv; per-sample: decontx_by_sample.tsv"
), PROV)
msg("wrote ", PROV)

if (!TEST) unlink(TMP, recursive = TRUE)
msg("done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
sessionInfo()
