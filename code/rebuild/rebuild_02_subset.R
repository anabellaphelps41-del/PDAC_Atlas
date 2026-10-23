#!/usr/bin/env Rscript
# rebuild_02_subset.R   (v3)
# Step 2 of the PDAC atlas rebuild: split the harmonized object into the three
# compartments that get their own pipeline (T/NK, stromal/CAF, B/plasma).
#
# Input : <DIR>/atlas_harmonized_counts_v3.rds   (rebuild_01_harmonize.R)
#         <DIR>/sample_manifest.tsv               (rebuild_00m_manifest.R)
#         <DIR>/doublets_scDblFinder_v3.tsv       (rebuild_01c_doublets.R)
# Output: subset_TNK.rds, subset_stromal.rds, subset_Bcell.rds  (+ PROVENANCE_02_subset.txt)
#
# Notes:
#   - Subsetting uses celltype_level2, NOT level1. level1/level2 disagree on
#     doublet status for 31,521 cells; level2 is the deliberate choice.
#   - "Mixed/possible doublet" (30,140) is its own level2 value and is therefore
#     excluded from all three subsets by construction, not by a filter.
#   - Cells are removed, in this order of attribution, when:
#       manifest      include_in_atlas == FALSE (blood, metastasis, non-PDAC,
#                     GSE211644 - the reason is carried per cell)
#       doublet       scDblFinder_class == "doublet"
#       neutrophil    celltype_level3 == "Neutrophil-like"
#   - Step 1 writes counts only. Each subset is normalized here
#     (LogNormalize, scale.factor 1e4). That is a per-cell operation, so the
#     values are identical to normalizing the whole atlas.
#   - percent.mt IS recomputed here, purely as a diagnostic, to establish which
#     datasets structurally lack MT genes.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

ROOT <- "/ix/tbruno/anp795/PDAC_Atlas"
DIR  <- Sys.getenv("PDAC_REBUILD_DIR", file.path(ROOT, "rebuild_v3"))
IN   <- file.path(DIR, "atlas_harmonized_counts_v3.rds")
MAN  <- file.path(DIR, "sample_manifest.tsv")
DBL  <- file.path(DIR, "doublets_scDblFinder_v3.tsv")
OUT  <- DIR

options(width = 200)
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
need <- c("celltype_level2", "celltype_level3", "dataset_id", "sample_id", "nCount_RNA")
missing <- setdiff(need, colnames(md))
if (length(missing)) {
  stop("required metadata column(s) absent: ", paste(missing, collapse = ", "))
}

# --------------------------------------------------------------------------
# Attach manifest + doublet calls
# --------------------------------------------------------------------------
man <- fread(MAN)
md$sample_key <- paste(md$dataset_id, md$sample_id, sep = ":")
mi <- match(md$sample_key, man$sample_key)
if (anyNA(mi)) stop("cells whose sample is not in the manifest: ", sum(is.na(mi)))
for (cc in c("patient_key", "patient_id_source", "tissue_class", "primary_or_metastatic",
             "patient_diagnosis", "treatment_status", "include_in_atlas", "exclude_reason")) {
  md[[cc]] <- man[[cc]][mi]
}

dbl <- fread(DBL)
di  <- match(rownames(md), dbl$cell)
if (anyNA(di)) stop("cells missing from the doublet table: ", sum(is.na(di)))
md$scDblFinder_score <- dbl$scDblFinder_score[di]
md$scDblFinder_class <- dbl$scDblFinder_class[di]
obj@meta.data <- md

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
} else {
  obj$percent.mt.recomputed <- NA_real_
  msg("no MT- genes at all. percent.mt is dead atlas-wide. Record and move on.")
}

# --------------------------------------------------------------------------
# Exclusions. Each removed cell is attributed to the first reason that applies.
# %in% rather than != so that NA labels are retained, not silently dropped.
# --------------------------------------------------------------------------
why <- rep("", nrow(md))
why[!md$include_in_atlas]                                   <- paste0("manifest: ", md$exclude_reason[!md$include_in_atlas])
why[!nzchar(why) & md$scDblFinder_class %in% "doublet"]     <- "scDblFinder doublet"
why[!nzchar(why) & md$celltype_level3 %in% "Neutrophil-like"] <- "Neutrophil-like (GSE155698 artifact)"
keep_global <- !nzchar(why)
msg("cells kept atlas-wide: ", format(sum(keep_global), big.mark = ","), " of ",
    format(nrow(md), big.mark = ","))

# --------------------------------------------------------------------------
# Subset definitions. Expected counts are the step-1 compartment tallies
# (unchanged since v2: same cells, same labels). A mismatch means the input
# object is not what this script thinks it is.
# --------------------------------------------------------------------------
defs <- list(
  TNK     = list(level2 = c("T/NK ambiguous", "T/NK"), expected = 193408L),
  stromal = list(level2 = "Stromal",                   expected =  41625L),
  Bcell   = list(level2 = "B/plasma",                  expected =  17975L)
)

written  <- list()
excluded <- list()

for (nm in names(defs)) {
  d      <- defs[[nm]]
  in_cmp <- md$celltype_level2 %in% d$level2
  cells  <- rownames(md)[in_cmp & keep_global]

  msg("---- ", nm, " ----")
  msg("compartment total: ", format(sum(in_cmp), big.mark = ","),
      " | expected: ",       format(d$expected,  big.mark = ","),
      " | after exclusions: ", format(length(cells), big.mark = ","))
  if (sum(in_cmp) != d$expected) {
    if (!nzchar(Sys.getenv("PDAC_ALLOW_COUNT_MISMATCH")))
      stop("compartment total does not match the expected count for ", nm)
    msg("*** count mismatch tolerated (PDAC_ALLOW_COUNT_MISMATCH set; test runs only) ***")
  }
  if (length(cells) == 0) stop("empty subset for ", nm, " - check level2 spelling")

  ex <- as.data.frame(table(reason = why[in_cmp & !keep_global]), stringsAsFactors = FALSE)
  ex <- ex[order(-ex$Freq), ]
  msg("excluded from ", nm, ":")
  print(ex, row.names = FALSE)
  excluded[[nm]] <- cbind(subset = nm, ex)

  s <- subset(obj, cells = cells)
  DefaultAssay(s) <- "RNA"
  s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

  msg("cells per dataset_id:")
  print(table(s$dataset_id))
  msg("cells per tissue_class:")
  print(table(s$tissue_class))

  f <- file.path(OUT, paste0("subset_", nm, ".rds"))
  msg("writing ", f, " (compress=FALSE)")
  saveRDS(s, f, compress = FALSE)
  msg("wrote ", round(file.size(f) / 2^30, 2), " GB")

  written[[nm]] <- data.frame(
    subset    = nm,
    criterion = paste0('celltype_level2 %in% c("', paste(d$level2, collapse = '", "'), '")'),
    compartment = sum(in_cmp),
    kept      = length(cells),
    genes     = nrow(s),
    file      = basename(f),
    stringsAsFactors = FALSE
  )

  rm(s); gc(verbose = FALSE)
}

# --------------------------------------------------------------------------
# Provenance file - answers "which cells, and why not the others?"
# --------------------------------------------------------------------------
prov <- file.path(OUT, "PROVENANCE_02_subset.txt")
tbl  <- do.call(rbind, written)
ext  <- do.call(rbind, excluded)
show <- function(x) capture.output(print(x, row.names = FALSE))

writeLines(c(
  "PDAC Atlas rebuild - compartment subsets (v3)",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("by:", "rebuild_02_subset.R"),
  paste("source:", IN),
  paste("manifest:", MAN),
  paste("doublets:", DBL),
  "",
  "ANNOTATION COLUMN USED FOR SUBSETTING: celltype_level2",
  "  (level1/level2 disagree on doublet status for 31,521 cells; level2 is the",
  "  deliberate choice.) 'Mixed/possible doublet' (30,140) is excluded by construction.",
  "",
  "EXCLUSIONS (first applicable reason per cell)",
  "  manifest     include_in_atlas == FALSE: blood, metastasis, non-PDAC tumour, GSE211644",
  "  doublet      scDblFinder_class == 'doublet' (per-sample, whole atlas)",
  "  neutrophil   celltype_level3 == 'Neutrophil-like' (85.6% GSE155698; single-dataset artifact)",
  "",
  show(ext),
  "",
  "Normalization: LogNormalize(scale.factor=1e4) per subset; counts unchanged from step 1.",
  "",
  "Subsets written:",
  show(tbl),
  "",
  "New metadata columns: sample_key, patient_key, patient_id_source, tissue_class,",
  "  primary_or_metastatic, patient_diagnosis, treatment_status, scDblFinder_score/class,",
  "  percent.mt.recomputed.",
  "",
  "Everything saved with compress=FALSE (loads in minutes; the filesystem compresses ~2.5x)."
), prov)

msg("wrote ", prov)
msg("done")
sessionInfo()
