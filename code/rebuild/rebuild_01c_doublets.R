#!/usr/bin/env Rscript
# =============================================================================
# rebuild_01c_doublets.R   --   scDblFinder per sample on the whole atlas
#
# Input : <DIR>/atlas_harmonized_counts_v3.rds
#         <DIR>/sample_manifest.tsv
# Output: <DIR>/doublets_scDblFinder_v3.tsv     one row per cell
#         <DIR>/doublets_by_sample_v3.tsv
#         <DIR>/PROVENANCE_01c_doublets.txt
#
# Runs on ALL retained cells (before compartment subsetting and before the
# manifest exclusions), per sample: doublets are a property of a library, and
# heterotypic doublets can only be simulated from the other cell types that
# were actually in that library. Samples are keyed dataset:sample because
# sample names repeat across datasets (P01-P10 in GSE154778 and GSE205013).
#
# Caveats recorded in provenance:
#   - every source dataset was already QC-filtered by its authors, some with
#     their own doublet removal; this catches what remains.
#   - GSE211644 was ingested as T cells only, so heterotypic doublets cannot be
#     simulated there; its calls are weak. It is excluded at step 2 anyway.
#   - samples below MIN_CELLS are not scored (class "not_scored").
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
  library(data.table)
})

ROOT <- "/ix/tbruno/anp795/PDAC_Atlas"
DIR  <- Sys.getenv("PDAC_REBUILD_DIR", file.path(ROOT, "rebuild_v3"))
IN   <- file.path(DIR, "atlas_harmonized_counts_v3.rds")
MAN  <- file.path(DIR, "sample_manifest.tsv")
OUT  <- file.path(DIR, "doublets_scDblFinder_v3.tsv")
OUTS <- file.path(DIR, "doublets_by_sample_v3.tsv")
PROV <- file.path(DIR, "PROVENANCE_01c_doublets.txt")

MIN_CELLS         <- 100
DBR_PER1K         <- 0.008    # scDblFinder default
CELLS_PER_CAPTURE <- 10000    # cap for the expected-rate calculation
SEED              <- 42
NTHREADS          <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))

options(width = 200)
msg <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = ""); flush.console() }

msg("loading ", IN)
obj <- readRDS(IN)
md  <- obj@meta.data
cts <- LayerData(obj, assay = "RNA", layer = "counts")
rm(obj); gc(verbose = FALSE)
msg("cells: ", format(ncol(cts), big.mark = ","), " | threads: ", NTHREADS)

skey  <- paste(md$dataset_id, md$sample_id, sep = ":")
man   <- fread(MAN)
miss  <- setdiff(unique(skey), man$sample_key)
if (length(miss)) stop("samples missing from manifest: ", paste(head(miss, 20), collapse = ", "))

n_per <- table(skey)
small <- names(n_per)[n_per < MIN_CELLS]
msg("samples: ", length(n_per), " | below ", MIN_CELLS, " cells (not scored): ", length(small),
    if (length(small)) paste0(" -> ", paste(small, collapse = ", ")) else "")

score <- names(n_per)[n_per >= MIN_CELLS]
sel   <- skey %in% score

# Expected doublet rate per sample. scDblFinder's default (0.8% per 1,000
# cells) assumes one 10x capture per sample; several deposited samples are
# pooled captures (GSE205013:P03 has 26,962 cells), which would ask for >20%
# doublets. Cap the cell count at one capture's worth.
dbr_of <- function(n) DBR_PER1K * min(n, CELLS_PER_CAPTURE) / 1000

msg("scDblFinder on ", format(sum(sel), big.mark = ","), " cells in ", length(score), " samples")
cols_by <- split(which(sel), skey[sel])
one <- function(k) {
  ii  <- cols_by[[k]]
  sce <- SingleCellExperiment(list(counts = cts[, ii]))
  set.seed(SEED)
  sce <- scDblFinder(sce, dbr = dbr_of(length(ii)))
  data.table(i = ii, score = sce$scDblFinder.score, class = as.character(sce$scDblFinder.class))
}
out <- rbindlist(bplapply(names(cols_by), one,
                          BPPARAM = MulticoreParam(workers = NTHREADS, RNGseed = SEED)))
rm(cts); gc(verbose = FALSE)

res <- data.table(cell = rownames(md), dataset_id = md$dataset_id, sample_key = skey,
                  celltype_level2 = md$celltype_level2,
                  scDblFinder_score = NA_real_, scDblFinder_class = "not_scored")
res[out$i, `:=`(scDblFinder_score = out$score, scDblFinder_class = out$class)]
fwrite(res, OUT, sep = "\t")
msg("wrote ", OUT)

by_s <- res[, .(n_cells = .N, n_doublet = sum(scDblFinder_class == "doublet"),
                pct_doublet = round(100 * mean(scDblFinder_class == "doublet"), 2)),
            by = .(dataset_id, sample_key)][order(dataset_id, sample_key)]
fwrite(by_s, OUTS, sep = "\t")

by_ds <- res[, .(n_cells = .N, n_doublet = sum(scDblFinder_class == "doublet"),
                 pct_doublet = round(100 * mean(scDblFinder_class == "doublet"), 2),
                 not_scored = sum(scDblFinder_class == "not_scored")), by = dataset_id][order(dataset_id)]
vs_l2 <- dcast(res, celltype_level2 ~ scDblFinder_class, fun.aggregate = length, value.var = "cell")
vs_l2[, pct_doublet := round(100 * doublet / (doublet + singlet), 1)]

print(by_ds); print(vs_l2)
show <- function(x) capture.output(print(x, row.names = FALSE))
writeLines(c(
  "PROVENANCE - doublet detection",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "script   : rebuild_01c_doublets.R",
  paste("input    :", IN),
  paste("output   :", OUT),
  sprintf("method   : scDblFinder %s, run per dataset:sample, dbr = %.3f per 1k cells (cells capped at %d), seed %d",
          as.character(packageVersion("scDblFinder")), DBR_PER1K, CELLS_PER_CAPTURE, SEED),
  sprintf("not scored: %d samples below %d cells (%d cells)", length(small), MIN_CELLS,
          sum(res$scDblFinder_class == "not_scored")),
  "",
  "BY DATASET", show(by_ds), "",
  "AGAINST INHERITED celltype_level2 (does 'Mixed/possible doublet' agree?)", show(vs_l2), "",
  "CAVEATS",
  "  - source datasets were QC-filtered by their authors before deposit.",
  "  - GSE211644 holds T cells only; heterotypic doublets cannot be simulated there.",
  "  - per-sample tables: doublets_by_sample_v3.tsv"
), PROV)
msg("wrote ", PROV)
msg("done")
sessionInfo()
