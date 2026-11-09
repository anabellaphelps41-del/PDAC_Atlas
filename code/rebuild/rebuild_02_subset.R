#!/usr/bin/env Rscript
# rebuild_02_subset.R   (v4)
# Step 2 of the PDAC atlas rebuild: split the harmonized atlas into the five
# compartments that get their own pipeline (T/NK, stromal, B/plasma, myeloid,
# endothelial), each carrying raw counts (RNA) and DecontX counts (decontX).
#
# Input : <UP>/atlas_harmonized_counts_v3.qs2   (rebuild_01_harmonize.R)
#         <UP>/sample_manifest.tsv               (rebuild_00m_manifest.R)
#         <UP>/doublets_scDblFinder_v3.tsv       (rebuild_01c_doublets.R)
#         <DIR>/decontx_counts.qs2, decontx_percell.tsv   (rebuild_01d_decontx.R)
#         <DIR>/cell_routing.tsv                 (rebuild_02a_routing.R)
# Output: <DIR>/subset_{TNK,stromal,Bcell,myeloid,endothelial}.qs2
#         <DIR>/PROVENANCE_02_subset.txt
# (UP = rebuild_v3, DIR = rebuild_v4; see io.R)
#
# Notes:
#   - Compartments start from celltype_level2, NOT level1 (level1/level2
#     disagree on doublet status for 31,521 cells; level2 is the deliberate
#     choice). 'Mixed/possible doublet' and 'Unassigned' are in no compartment.
#   - Then the reviewed res-0.6 calls (cell_routing.tsv) move misrouted cells
#     between subsets and drop the clusters marked drop. Routing is by barcode.
#   - Cells are removed, first applicable reason, when:
#       manifest   include_in_atlas == FALSE (blood, metastasis, non-PDAC, GSE211644)
#       doublet    scDblFinder_class == "doublet"
#       qc_floor   below the per dataset x compartment floor (rebuild_01d)
#       reviewed   cluster call 'drop' (code/calls)
#   - CHANGED FROM v3: celltype_level3 'Neutrophil-like' cells are no longer
#     removed. v3 dropped them as a GSE155698 artifact (85.6% of them), but
#     that figure counted blood samples. In the included samples GSE155698 has
#     57%, GSE212966 27% and GSE205013 15%, and the myeloid subset now exists
#     to hold them. Its lineage gate decides about them like any other cells.
#   - Both assays are LogNormalized (scale.factor 1e4). decontX is the DEFAULT
#     assay: steps 02b-04b select genes, cluster, find markers and call SingleR
#     on decontaminated data. RNA keeps the raw counts; nCount_RNA and
#     nFeature_RNA (raw) remain the QC metrics.
#   - percent.mt is recomputed on raw counts as a diagnostic only.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

source("/ix/tbruno/anp795/PDAC_Atlas/code/rebuild/io.R")   # read_obj(), save_obj(), dirs

UP   <- UP_DIR
DIR  <- V4_DIR
IN   <- file.path(UP, "atlas_harmonized_counts_v3.qs2")
MAN  <- file.path(UP, "sample_manifest.tsv")
DBL  <- file.path(UP, "doublets_scDblFinder_v3.tsv")
DXM  <- file.path(DIR, "decontx_counts.qs2")
DXC  <- file.path(DIR, "decontx_percell.tsv")
ROUT <- file.path(DIR, "cell_routing.tsv")
OUT  <- DIR

options(width = 200)
msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

# --------------------------------------------------------------------------
# Subset definitions. Expected counts are the atlas-wide level2 tallies after
# step 1 (538,766 cells). A mismatch means the input object is not what this
# script thinks it is.
# --------------------------------------------------------------------------
defs <- list(
  TNK         = list(level2 = c("T/NK ambiguous", "T/NK"), expected = 193408L),
  stromal     = list(level2 = "Stromal",                   expected =  41625L),
  Bcell       = list(level2 = "B/plasma",                  expected =  17975L),
  myeloid     = list(level2 = c("Myeloid", "Mast"),        expected =  95284L),
  endothelial = list(level2 = "Vascular",                  expected =  15201L)
)
stopifnot(identical(names(defs), SUBSETS))
V3_SUBSETS <- c("TNK", "stromal", "Bcell")   # the ones the reviewed calls cover

# --------------------------------------------------------------------------
# Load and validate
# --------------------------------------------------------------------------
msg("loading ", IN)
obj <- read_obj(IN)
msg("loaded: ", format(nrow(obj), big.mark = ","), " genes x ",
    format(ncol(obj), big.mark = ","), " cells")

md <- obj@meta.data
need <- c("celltype_level2", "celltype_level3", "dataset_id", "sample_id", "nCount_RNA")
missing <- setdiff(need, colnames(md))
if (length(missing)) stop("required metadata column(s) absent: ", paste(missing, collapse = ", "))

man <- fread(MAN)
md$sample_key <- paste(md$dataset_id, md$sample_id, sep = ":")
mi <- match(md$sample_key, man$sample_key)
if (anyNA(mi)) stop("cells whose sample is not in the manifest: ", sum(is.na(mi)))
for (cc in c("patient_key", "patient_id_source", "tissue_class", "primary_or_metastatic",
             "patient_diagnosis", "treatment_status", "include_in_atlas", "exclude_reason")) {
  md[[cc]] <- man[[cc]][mi]
}
md$include_in_atlas <- as.logical(md$include_in_atlas)

dbl <- fread(DBL)
di  <- match(rownames(md), dbl$cell)
if (anyNA(di)) stop("cells missing from the doublet table: ", sum(is.na(di)))
md$scDblFinder_score <- dbl$scDblFinder_score[di]
md$scDblFinder_class <- dbl$scDblFinder_class[di]

dxc <- fread(DXC)
xi  <- match(rownames(md), dxc$cell)                     # NA = out of 01d scope
for (cc in c("decontX_contamination", "decontX_status", "qc_pass",
             "qc_floor_log10_nFeature", "qc_floor_log10_nCount", "nCount_decontX")) {
  md[[cc]] <- dxc[[cc]][xi]
}

rt <- fread(ROUT, colClasses = list(character = c("to_subset", "followup", "celltype_label")))
ri <- match(rownames(md), rt$cell)
md$v3_subset         <- rt$from_subset[ri]
md$v3_cluster        <- rt$v3_cluster[ri]
md$v3_action         <- rt$action[ri]
md$v3_celltype_label <- rt$celltype_label[ri]
md$v3_followup       <- rt$followup[ri]
md$v3_to_subset      <- rt$to_subset[ri]

# Compartment from level2, then the routing
cmp <- rep(NA_character_, nrow(md))
for (nm in names(defs)) cmp[md$celltype_level2 %in% defs[[nm]]$level2] <- nm
md$compartment_level2 <- cmp

# --------------------------------------------------------------------------
# Diagnostic: percent.mt on raw counts
# --------------------------------------------------------------------------
mt_genes <- grep("^MT-", rownames(obj), value = TRUE)
msg("MT- genes present in harmonized matrix: ", length(mt_genes))
if (length(mt_genes) > 0) {
  mt_counts <- Matrix::colSums(LayerData(obj, layer = "counts")[mt_genes, , drop = FALSE])
  pmt <- 100 * mt_counts / md$nCount_RNA
  pmt[!is.finite(pmt)] <- NA_real_
  md$percent.mt.recomputed <- pmt
} else {
  md$percent.mt.recomputed <- NA_real_
}
obj@meta.data <- md

# --------------------------------------------------------------------------
# Exclusions. Each removed cell is attributed to the first reason that applies.
# --------------------------------------------------------------------------
why <- rep("", nrow(md))
why[!md$include_in_atlas] <- paste0("manifest: ", md$exclude_reason[!md$include_in_atlas])
why[!nzchar(why) & md$scDblFinder_class %in% "doublet"] <- "scDblFinder doublet"
in_any <- !is.na(cmp) | (!is.na(md$v3_to_subset) & nzchar(md$v3_to_subset))
lost <- !nzchar(why) & in_any & is.na(md$qc_pass)
if (any(lost)) stop(sum(lost), " compartment cells are missing from decontx_percell.tsv")
why[!nzchar(why) & md$qc_pass %in% FALSE] <- "QC floor (rebuild_01d)"
rdrop <- !nzchar(why) & md$v3_action %in% "drop"
why[rdrop] <- paste0("reviewed call: drop (v3 ", md$v3_subset[rdrop], " cluster ", md$v3_cluster[rdrop], ")")
keep_global <- !nzchar(why)

# Final subset membership: level2 compartment, overridden by routing
dest <- cmp
routed <- !is.na(md$v3_action) & md$v3_action == "route"
dest[routed] <- md$v3_to_subset[routed]
md$subset_v4 <- ifelse(keep_global, dest, NA_character_)
obj$subset_v4 <- md$subset_v4

# Cells of a v3 subset with no routing row (dropped at v3 step 3, never reviewed)
unrev <- keep_global & cmp %in% V3_SUBSETS & is.na(md$v3_action)
unrev_tab <- as.data.frame(table(subset = cmp[unrev], dataset_id = md$dataset_id[unrev]))
unrev_tab <- unrev_tab[unrev_tab$Freq > 0, ]
msg("cells in v3 subsets without a reviewed call (kept, unreviewed): ", sum(unrev))
if (nrow(unrev_tab)) print(unrev_tab, row.names = FALSE)
if (sum(unrev) > 500) stop("more unreviewed cells than expected (the v3 step-3 batch drop was 49); ",
                           "the routing table may not match this atlas")
md$v3_action[unrev] <- "unreviewed"
obj$v3_action <- md$v3_action

msg("cells kept atlas-wide (in some subset): ",
    format(sum(!is.na(md$subset_v4)), big.mark = ","))

# --------------------------------------------------------------------------
# DecontX matrix
# --------------------------------------------------------------------------
msg("loading ", DXM)
dx <- read_obj(DXM)
if (!identical(rownames(dx), rownames(obj))) stop("decontX genes do not match the atlas")
msg("decontX matrix: ", format(nrow(dx), big.mark = ","), " x ", format(ncol(dx), big.mark = ","))

# --------------------------------------------------------------------------
# Build subsets
# --------------------------------------------------------------------------
written <- list(); excluded <- list(); flows <- list()

for (nm in names(defs)) {
  d      <- defs[[nm]]
  in_cmp <- cmp %in% nm
  cells  <- rownames(md)[md$subset_v4 %in% nm]

  msg("---- ", nm, " ----")
  msg("level2 compartment: ", format(sum(in_cmp), big.mark = ","),
      " | expected: ", format(d$expected, big.mark = ","),
      " | final subset: ", format(length(cells), big.mark = ","))
  if (sum(in_cmp) != d$expected) {
    if (!nzchar(Sys.getenv("PDAC_ALLOW_COUNT_MISMATCH")))
      stop("compartment total does not match the expected count for ", nm)
    msg("*** count mismatch tolerated (PDAC_ALLOW_COUNT_MISMATCH set; test runs only) ***")
  }
  if (length(cells) == 0) stop("empty subset for ", nm)

  ex <- as.data.frame(table(reason = why[in_cmp & !keep_global]), stringsAsFactors = FALSE)
  ex <- ex[order(-ex$Freq), ]
  excluded[[nm]] <- if (nrow(ex)) cbind(subset = nm, ex) else NULL
  print(ex, row.names = FALSE)

  # routing flows into and out of this subset
  out_n <- sum(in_cmp & keep_global & routed & md$v3_to_subset != nm)
  in_n  <- sum(!in_cmp & md$subset_v4 %in% nm)
  flows[[nm]] <- data.frame(subset = nm, level2_after_exclusions = sum(in_cmp & keep_global),
                            routed_out = out_n, routed_in = in_n, final = length(cells))

  missing_dx <- setdiff(cells, colnames(dx))
  if (length(missing_dx)) stop(length(missing_dx), " ", nm, " cells have no DecontX counts")

  s <- subset(obj, cells = cells)
  DefaultAssay(s) <- "RNA"
  s <- NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  s[["decontX"]] <- CreateAssay5Object(counts = dx[, colnames(s)])
  s <- NormalizeData(s, assay = "decontX", normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)
  s$nCount_decontX   <- Matrix::colSums(LayerData(s, assay = "decontX", layer = "counts"))
  s$nFeature_decontX <- Matrix::colSums(LayerData(s, assay = "decontX", layer = "counts") > 0)
  DefaultAssay(s) <- "decontX"

  msg("cells per dataset_id:"); print(table(s$dataset_id))
  msg("cells per v3 decision:"); print(table(s$v3_action, useNA = "ifany"))
  msg("median DecontX contamination: ", round(median(s$decontX_contamination, na.rm = TRUE), 4))

  f <- file.path(OUT, paste0("subset_", nm, ".qs2"))
  msg("writing ", f, " (qs2)")
  save_obj(s, f)
  msg("wrote ", round(file.size(f) / 2^30, 2), " GB")

  written[[nm]] <- data.frame(
    subset = nm,
    level2 = paste(d$level2, collapse = " + "),
    compartment = sum(in_cmp),
    kept = length(cells),
    median_contamination = round(median(s$decontX_contamination, na.rm = TRUE), 4),
    file = basename(f), stringsAsFactors = FALSE)
  rm(s); gc(verbose = FALSE)
}

# --------------------------------------------------------------------------
# Provenance
# --------------------------------------------------------------------------
prov <- file.path(OUT, "PROVENANCE_02_subset.txt")
show <- function(x) capture.output(print(x, row.names = FALSE))

writeLines(c(
  "PDAC Atlas rebuild - compartment subsets (v4)",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "by: rebuild_02_subset.R",
  paste("source:", IN), paste("manifest:", MAN), paste("doublets:", DBL),
  paste("decontX:", DXM), paste("routing:", ROUT),
  "",
  "COMPARTMENTS: celltype_level2, then the reviewed res-0.6 calls (code/calls)",
  "  move misrouted cells by barcode and drop the clusters marked drop.",
  "  'Mixed/possible doublet' and 'Unassigned' are in no compartment.",
  "",
  "EXCLUSIONS (first applicable reason per cell)",
  "  manifest   include_in_atlas == FALSE: blood, metastasis, non-PDAC tumour, GSE211644",
  "  doublet    scDblFinder_class == 'doublet' (per sample, whole atlas)",
  "  qc_floor   below median - 3 MAD of log10 nFeature/nCount, per dataset x compartment (01d)",
  "  reviewed   cluster call 'drop' in code/calls/*_reviewed.tsv",
  "  NOT excluded any more: celltype_level3 'Neutrophil-like' (see script header).",
  "",
  show(do.call(rbind, excluded)),
  "",
  "ROUTING FLOWS", show(do.call(rbind, flows)),
  "",
  sprintf("UNREVIEWED: %d cells of v3 subsets had no reviewed call (kept):", sum(unrev)),
  if (nrow(unrev_tab)) show(unrev_tab) else "  none",
  "",
  "ASSAYS: RNA = raw counts; decontX = rounded DecontX counts (DEFAULT assay).",
  "  Both LogNormalize(scale.factor=1e4). QC metrics (nCount_RNA, nFeature_RNA)",
  "  are raw; nCount_decontX / nFeature_decontX and decontX_contamination are added.",
  "",
  "Subsets written:", show(do.call(rbind, written)),
  "",
  "Metadata added: sample_key, patient_key, patient_id_source, tissue_class,",
  "  primary_or_metastatic, patient_diagnosis, treatment_status, scDblFinder_score/class,",
  "  decontX_contamination/status, qc_pass + floors, compartment_level2, subset_v4,",
  "  v3_subset, v3_cluster, v3_action, v3_to_subset, v3_celltype_label, v3_followup,",
  "  percent.mt.recomputed."
), prov)

msg("wrote ", prov)
msg("done")
sessionInfo()
