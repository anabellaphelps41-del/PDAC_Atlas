#!/usr/bin/env Rscript
# =============================================================================
# rebuild_04_markers.R   --   Step 4 of the PDAC atlas rebuild
#
# Marker discovery + contamination triage for one Harmony-integrated subset.
#
# Usage:  Rscript rebuild_04_markers.R <subset> [resolution]
#         <subset>    in {stromal, TNK, Bcell}   (default: stromal)
#         [resolution] one of the swept values   (default: 0.6)
#
# Input : subset_<subset>_harmony.rds  (from rebuild_03_pipeline.R)
# Output: step4_<subset>_res<r>_allmarkers.tsv     full FindAllMarkers table
#         step4_<subset>_res<r>_top20.tsv          top 20 per cluster
#         step4_<subset>_res<r>_decision.tsv  <<<  ONE ROW PER CLUSTER, the
#                                                  table you actually annotate from
#         step4_<subset>_res<r>_qc.pdf
#         PROVENANCE_step4_<subset>_res<r>.txt
#
# WHAT THIS DOES NOT DO: it does not name cell types. It assembles the evidence
# needed to name them -- markers, lineage purity, dataset composition, QC -- so
# that every call is made against numbers rather than against a UMAP shape.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(ggplot2)
})

ts <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                paste0(..., collapse = "")))

if (!requireNamespace("presto", quietly = TRUE)) {
  ts("NOTE: 'presto' is not installed. FindAllMarkers will use the slow Wilcoxon path.")
  ts("      On a LOGIN node, this is worth 10-50x on the big subsets:")
  ts("        R -e 'install.packages(\"presto\", repos=c(\"https://immunogenomics.r-universe.dev\",\"https://cloud.r-project.org\"))'")
}

# =============================== parameters ==================================
args   <- commandArgs(trailingOnly = TRUE)
SUBSET <- if (length(args) >= 1) args[1] else "stromal"
RES    <- if (length(args) >= 2) as.numeric(args[2]) else 0.6
stopifnot(SUBSET %in% c("stromal", "TNK", "Bcell"))

DIR      <- Sys.getenv("PDAC_REBUILD_DIR", "/ix/tbruno/anp795/PDAC_Atlas/rebuild_v3")
IN       <- file.path(DIR, sprintf("subset_%s_harmony.rds", SUBSET))
tag      <- sprintf("%s_res%s", SUBSET, RES)
F_ALL    <- file.path(DIR, sprintf("step4_%s_allmarkers.tsv", tag))
F_TOP    <- file.path(DIR, sprintf("step4_%s_top20.tsv", tag))
F_DEC    <- file.path(DIR, sprintf("step4_%s_decision.tsv", tag))
F_PDF    <- file.path(DIR, sprintf("step4_%s_qc.pdf", tag))
F_PROV   <- file.path(DIR, sprintf("PROVENANCE_step4_%s.txt", tag))

MAX_CELLS_PER_IDENT <- 5000   # FindAllMarkers subsample; see provenance note
MIN_PCT             <- 0.25
LOGFC               <- 0.25
PURITY_FLAG         <- 0.80   # cluster >this fraction one dataset -> flag
FOREIGN_PCT_FLAG    <- 40     # foreign lineage detected in >this % of cells -> flag
NATIVE_PCT_FLAG     <- 40     # native lineage below this % -> flag
SEED                <- 42
set.seed(SEED)

# ===================== lineage panels (pct-expressed, not markers) ===========
# These are used to ask "is this cluster the lineage the subset claims?" -- a
# purity question, separate from marker discovery.
lineage <- list(
  fibroblast  = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  mural       = c("RGS5", "ACTA2", "NOTCH3", "MCAM", "MYH11", "CNN1"),
  myeloid     = c("PTPRC", "CD68", "LYZ", "AIF1", "C1QA", "C1QB", "CSF1R"),
  Tcell       = c("CD3D", "CD3E", "CD3G", "CD2", "IL7R"),
  NK          = c("NKG7", "GNLY", "KLRD1", "FCGR3A"),
  Bcell       = c("MS4A1", "CD79A", "CD79B", "BANK1"),
  plasma      = c("MZB1", "XBP1", "TNFRSF17", "DERL3"),
  epithelial  = c("EPCAM", "KRT19", "KRT8", "KRT18", "MUC1"),
  endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  RBC         = c("HBB", "HBA1", "HBA2"),
  acinar      = c("PRSS1", "CTRB1", "CPA1", "REG1A"),   # ambient tracer, not a cell type
  MHCII       = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD74")
)
NATIVE <- list(stromal = "fibroblast", TNK = "Tcell", Bcell = "Bcell")[[SUBSET]]

# ===================== module score programmes ==============================
programmes <- list(
  stromal = list(
    myCAF = c("ACTA2","TAGLN","MYL9","POSTN","COL10A1","COL11A1","THBS2","CTHRC1","INHBA"),
    iCAF  = c("IL6","CXCL12","CXCL14","PDGFRA","CFD","C3","C7","HAS1","LMNA","DPT","APOD"),
    apCAF = c("CD74","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","SLPI","SAA1"),
    stellate = c("LRAT","NGFR","CYGB","DES","RGS5","PDGFRB")
  ),
  TNK = list(
    naive     = c("CCR7","SELL","TCF7","LEF1"),
    cytotoxic = c("GZMA","GZMB","GZMK","GZMH","PRF1","NKG7","GNLY"),
    Treg      = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF4","TNFRSF18"),
    exhausted = c("PDCD1","LAG3","HAVCR2","TIGIT","TOX","ENTPD1","CXCL13"),
    NKcell    = c("NCAM1","FCGR3A","KLRD1","TYROBP","FCER1G","NCR1"),
    prolif    = c("MKI67","TOP2A","STMN1","TYMS")
  ),
  Bcell = list(
    naive   = c("TCL1A","FCER2","IL4R"),
    memory  = c("CD27","AIM2","TNFRSF13B","CD80"),
    GC      = c("BCL6","AICDA","RGS13","MEF2B","LMO2"),
    plasma  = c("MZB1","XBP1","SDC1","TNFRSF17","PRDM1","DERL3","JCHAIN"),
    prolif  = c("MKI67","TOP2A","STMN1")
  )
)[[SUBSET]]

# ============================ 1. load ========================================
ts("subset = ", SUBSET, " | resolution = ", RES)
ts("loading ", IN)
sub <- readRDS(IN)
ts("loaded: ", format(nrow(sub), big.mark = ","), " genes x ",
   format(ncol(sub), big.mark = ","), " cells")

res_col <- paste0("RNA_snn_res.", RES)
if (!res_col %in% colnames(sub@meta.data)) {
  stop(sprintf("'%s' not in metadata. Available: %s", res_col,
               paste(grep("^RNA_snn_res", colnames(sub@meta.data), value = TRUE),
                     collapse = ", ")))
}
Idents(sub) <- res_col
sub$cluster <- Idents(sub)
nclust <- nlevels(Idents(sub))
ts("clusters at res ", RES, ": ", nclust)

# ============================ 2. FindAllMarkers ==============================
ts("FindAllMarkers (wilcox, only.pos, min.pct=", MIN_PCT, ", logfc=", LOGFC,
   ", max.cells.per.ident=", MAX_CELLS_PER_IDENT, ")")
mk <- FindAllMarkers(sub, only.pos = TRUE, min.pct = MIN_PCT,
                     logfc.threshold = LOGFC, test.use = "wilcox",
                     max.cells.per.ident = MAX_CELLS_PER_IDENT,
                     random.seed = SEED, verbose = FALSE)
mk <- mk[order(mk$cluster, -mk$avg_log2FC), ]
write.table(mk, F_ALL, sep = "\t", quote = FALSE, row.names = FALSE)
ts("wrote ", F_ALL, " (", nrow(mk), " rows)")

top20 <- do.call(rbind, lapply(split(mk, mk$cluster), function(d) head(d, 20)))
write.table(top20, F_TOP, sep = "\t", quote = FALSE, row.names = FALSE)
ts("wrote ", F_TOP)

ts("top 12 markers per cluster:")
for (cl in levels(Idents(sub))) {
  g <- head(mk$gene[mk$cluster == cl], 12)
  cat(sprintf("  cluster %-3s n=%-7s %s\n", cl,
              format(sum(Idents(sub) == cl), big.mark = ","),
              paste(g, collapse = ", ")))
}

# ============================ 3. lineage purity ==============================
ts("computing per-cluster lineage detection rates")
dat  <- GetAssayData(sub, assay = "RNA", layer = "data")
idx  <- split(seq_len(ncol(sub)), Idents(sub))

pct_panel <- function(genes) {
  genes <- intersect(genes, rownames(sub))
  if (length(genes) == 0) return(setNames(rep(NA_real_, length(idx)), names(idx)))
  d <- dat[genes, , drop = FALSE]
  vapply(idx, function(ii) mean(Matrix::rowMeans(d[, ii, drop = FALSE] > 0)) * 100,
         numeric(1))
}
lin <- as.data.frame(lapply(lineage, pct_panel))
rownames(lin) <- names(idx)
ts("per-cluster lineage detection (% of cells, panel mean):")
print(round(lin, 1))

# ============================ 4. module scores ===============================
ts("AddModuleScore for ", length(programmes), " programmes")
for (nm in names(programmes)) {
  gset <- intersect(programmes[[nm]], rownames(sub))
  if (length(gset) < 3) { ts("  skipping ", nm, " (only ", length(gset), " genes present)"); next }
  sub <- AddModuleScore(sub, features = list(gset), name = paste0("score_", nm),
                        seed = SEED, ctrl = min(50, length(gset) * 10))
  colnames(sub@meta.data)[colnames(sub@meta.data) == paste0("score_", nm, "1")] <-
    paste0("score_", nm)
}
score_cols <- grep("^score_", colnames(sub@meta.data), value = TRUE)
score_by_cl <- if (length(score_cols)) {
  as.data.frame(lapply(score_cols, function(cc)
    vapply(idx, function(ii) mean(sub@meta.data[[cc]][ii]), numeric(1))))
} else data.frame(row.names = names(idx))
if (ncol(score_by_cl)) colnames(score_by_cl) <- sub("^score_", "", score_cols)
if (ncol(score_by_cl)) {
  ts("mean programme score per cluster:")
  print(round(score_by_cl, 3))
}

# ============================ 5. decision table ==============================
ts("assembling decision table")
comp      <- prop.table(table(Idents(sub), sub$dataset_id), 1)
top_ds    <- colnames(comp)[apply(comp, 1, which.max)]
top_ds_f  <- apply(comp, 1, max)
med <- function(v) vapply(idx, function(ii) median(v[ii], na.rm = TRUE), numeric(1))

dec <- data.frame(
  cluster          = names(idx),
  n_cells          = vapply(idx, length, integer(1)),
  pct_of_subset    = round(100 * vapply(idx, length, integer(1)) / ncol(sub), 2),
  dominant_dataset = top_ds,
  dominant_frac    = round(as.numeric(top_ds_f), 3),
  median_nFeature  = round(med(sub$nFeature_RNA)),
  median_nCount    = round(med(sub$nCount_RNA)),
  median_pct_mt    = if ("percent.mt" %in% colnames(sub@meta.data))
                       round(med(sub$percent.mt), 2) else NA,
  round(lin, 1),
  round(score_by_cl, 3),
  top_markers      = vapply(names(idx), function(cl)
                       paste(head(mk$gene[mk$cluster == cl], 8), collapse = ";"),
                       character(1)),
  check.names = FALSE, row.names = NULL
)

# --- automated flags. These are prompts to look, not verdicts. --------------
foreign <- setdiff(names(lineage), c(NATIVE, "MHCII", "acinar"))
dec$flag_batch      <- ifelse(dec$dominant_frac > PURITY_FLAG,
                              sprintf("PURITY>%.0f%%_%s", 100 * PURITY_FLAG, dec$dominant_dataset), "")
dec$flag_lowcomplex <- ifelse(dec$median_nFeature < 1000, "LOW_COMPLEXITY", "")
dec$flag_highcount  <- ifelse(dec$median_nCount > 2.5 * median(sub$nCount_RNA), "HIGH_UMI_DOUBLET?", "")
dec$flag_mt         <- ifelse(!is.na(dec$median_pct_mt) & dec$median_pct_mt > 10, "HIGH_MT", "")
dec$flag_mt_zero    <- ifelse(!is.na(dec$median_pct_mt) & dec$median_pct_mt < 0.5, "MT_NEAR_ZERO_QC_ARTEFACT?", "")
dec$flag_foreign    <- vapply(seq_len(nrow(dec)), function(i) {
  hits <- foreign[!is.na(as.numeric(lin[i, foreign])) &
                    as.numeric(lin[i, foreign]) > FOREIGN_PCT_FLAG]
  nat  <- as.numeric(lin[i, NATIVE])
  if (length(hits) && !is.na(nat) && nat < NATIVE_PCT_FLAG)
    paste0("FOREIGN:", paste(hits, collapse = "+")) else ""
}, character(1))
dec$flags <- trimws(apply(dec[, grep("^flag_", names(dec))], 1,
                          function(r) paste(r[nzchar(r)], collapse = " | ")))
dec <- dec[, !grepl("^flag_[a-z]", names(dec))]

write.table(dec, F_DEC, sep = "\t", quote = FALSE, row.names = FALSE)
ts("wrote ", F_DEC)

ts("FLAGGED CLUSTERS (inspect before naming):")
fl <- dec[nzchar(dec$flags), c("cluster","n_cells","pct_of_subset","dominant_dataset",
                               "dominant_frac","median_nFeature","median_pct_mt","flags")]
if (nrow(fl)) print(fl, row.names = FALSE) else cat("   none\n")
cat("\n   READ: a flag is a prompt to look, not a verdict. A FOREIGN flag on a\n")
cat("   MHC-II-high cluster with low collagen is the classic apCAF false\n")
cat("   positive -- apCAF must be MHC-II positive AND collagen positive in the\n")
cat("   SAME cells. A PURITY flag on a cluster whose datasets share a reference\n")
cat("   build is batch; on a single dataset it may be cohort biology.\n\n")

# ============================ 6. plots =======================================
ts("writing ", F_PDF)
pdf(F_PDF, width = 12, height = 8.5)

print(DimPlot(sub, reduction = "umap", label = TRUE, repel = TRUE,
              raster = TRUE, pt.size = 0.3) + NoLegend() +
        ggtitle(sprintf("%s - res %s (%d clusters)", SUBSET, RES, nclust)))

lin_long <- data.frame(
  cluster = factor(rep(rownames(lin), times = ncol(lin)), levels = rownames(lin)),
  panel   = factor(rep(colnames(lin), each = nrow(lin)), levels = colnames(lin)),
  pct     = as.numeric(as.matrix(lin))
)
print(ggplot(lin_long, aes(panel, cluster, fill = pct)) +
        geom_tile(colour = "grey90") +
        geom_text(aes(label = ifelse(is.na(pct), "", sprintf("%.0f", pct))), size = 2.6) +
        scale_fill_gradient(low = "white", high = "#2166AC", na.value = "grey95",
                            limits = c(0, 100), name = "% cells") +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = sprintf("%s res %s - lineage detection rate per cluster", SUBSET, RES),
             subtitle = sprintf("native lineage for this subset: %s", NATIVE),
             x = NULL, y = "cluster"))

for (cc in score_cols) {
  print(FeaturePlot(sub, features = cc, reduction = "umap", raster = TRUE,
                    pt.size = 0.3, order = TRUE) +
          scale_colour_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B",
                                 midpoint = 0) +
          ggtitle(sprintf("%s - %s programme score", SUBSET, sub("^score_", "", cc))))
  print(VlnPlot(sub, features = cc, pt.size = 0) + NoLegend() +
          ggtitle(sprintf("%s - %s by cluster", SUBSET, sub("^score_", "", cc))))
}

for (nm in names(lineage)) {
  g <- intersect(lineage[[nm]], rownames(sub))
  if (length(g) == 0) next
  print(DotPlot(sub, features = g) + RotatedAxis() +
          ggtitle(sprintf("%s res %s - %s panel", SUBSET, RES, nm)))
}

topg <- unique(unlist(lapply(split(mk, mk$cluster), function(d) head(d$gene, 3))))
topg <- head(intersect(topg, rownames(sub)), 60)
if (length(topg)) {
  print(DotPlot(sub, features = topg) + RotatedAxis() +
          theme(axis.text.x = element_text(size = 6)) +
          ggtitle(sprintf("%s res %s - top 3 markers per cluster", SUBSET, RES)))
}

comp_long <- data.frame(
  cluster = factor(rep(rownames(comp), times = ncol(comp)), levels = rownames(comp)),
  dataset = rep(colnames(comp), each = nrow(comp)),
  frac    = as.numeric(comp)
)
print(ggplot(comp_long, aes(cluster, frac, fill = dataset)) +
        geom_col() + theme_minimal(base_size = 11) +
        geom_hline(yintercept = PURITY_FLAG, linetype = 2) +
        labs(title = sprintf("%s res %s - dataset composition per cluster", SUBSET, RES),
             subtitle = sprintf("dashed line = %.0f%% purity flag", 100 * PURITY_FLAG),
             y = "fraction of cluster", x = "cluster"))

dev.off()

# ============================ 7. save + provenance ===========================
saveRDS(sub, IN, compress = FALSE)   # re-save with module scores + cluster col
ts("re-wrote ", IN, " with programme scores attached")

writeLines(c(
  sprintf("PROVENANCE - Step 4 markers - subset '%s', resolution %s", SUBSET, RES),
  sprintf("generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("script   : rebuild_04_markers.R"),
  sprintf("input    : %s", IN),
  "",
  "METHOD",
  sprintf("  FindAllMarkers: wilcox, only.pos=TRUE, min.pct=%.2f, logfc.threshold=%.2f",
          MIN_PCT, LOGFC),
  sprintf("  max.cells.per.ident=%d -- clusters larger than this are SUBSAMPLED.", MAX_CELLS_PER_IDENT),
  "    This bounds runtime and has little effect on ranking at these cluster",
  "    sizes, but it means p-values are not comparable across clusters of very",
  "    different size. Rank by avg_log2FC and pct.1-pct.2, not by p_val_adj.",
  sprintf("  seed=%d", SEED),
  "",
  "FLAG THRESHOLDS (prompts to inspect, not verdicts)",
  sprintf("  dataset purity > %.0f%%", 100 * PURITY_FLAG),
  sprintf("  foreign lineage detected in > %d%% of cells while native (%s) < %d%%",
          FOREIGN_PCT_FLAG, NATIVE, NATIVE_PCT_FLAG),
  "  median nFeature < 1000 (low complexity / ambient)",
  "  median nCount > 2.5x subset median (doublet candidate)",
  "  median percent.mt > 10 (dying) or < 0.5 (broken MT column, see GSE154778)",
  "",
  "CONSTRAINTS CARRIED FORWARD",
  "  - apCAF requires MHC-II AND collagen in the same cells. MHC-II-high with",
  "    low collagen is myeloid contamination.",
  "  - No gamma-delta calls in CRA001160 or GSE155698 (TR constant regions absent).",
  "  - No cross-dataset isotype work (IG constant regions absent in the same two).",
  "  - No compositional claim may rest on GSE211644 (sorted, CD8-biased).",
  "  - Masked HVG genes are still present and plottable; the mask governed",
  "    feature selection only.",
  ""
), F_PROV)
ts("wrote ", F_PROV)

ts("done")
print(sessionInfo())
