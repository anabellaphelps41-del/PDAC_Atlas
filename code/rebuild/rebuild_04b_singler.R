#!/usr/bin/env Rscript
# =============================================================================
# rebuild_04b_singler.R   --   reference-based lineage call + the lineage gate
#
# Usage:  Rscript rebuild_04b_singler.R <subset> [resolution]
#         <subset> in {TNK, stromal, Bcell, myeloid, endothelial}   [resolution] default 0.6
#         SingleR runs on the default assay (decontX in v4, RNA in v3).
#
# Input : subset_<subset>_harmony.qs2   (already carries the resolution columns)
#         ref/celldex_HPCA.rds, ref/celldex_BlueprintEncode.rds
#         (cache them on a LOGIN node first: Rscript ~/setup_singler_refs.R)
# Output: step4b_<subset>_res<r>_singler_percell.tsv
#         step4b_<subset>_res<r>_gate.tsv     <<< the table you act on
#         step4b_<subset>_res<r>_singler.pdf
#         PROVENANCE_step4b_<subset>_res<r>.txt
#         re-saves subset_<subset>_harmony.qs2 with singler_* columns attached
#
# -----------------------------------------------------------------------------
# WHY THIS EXISTS
#
# Decision 9 in handoff v4 section 5 specifies a lineage gate: score apCAF only
# in clusters that PASS the fibroblast panel and FAIL myeloid, endothelial and
# epithelial. The positive half of that gate cannot work on this data.
#
# Measured on step4_stromal_res0.6_decision.tsv: the fibroblast panel
# (COL1A1, COL1A2, COL3A1, DCN, LUM) is detected in >=60% of cells in 17 of 20
# clusters -- 84.5% in the macrophage cluster, 76.7% in the T cell cluster,
# 96.8% in the ductal epithelium cluster, and 80.6% in the low-complexity
# artifact. Ambient collagen is atlas-wide, and only 2 of the 6 source datasets
# had ambient RNA removed (GSE205013 and GSE212966, both SoupX; the source
# methods audit of 2026-08-29 has the table).
#
# The downstream consequence is visible in rebuild_04_markers.R lines 217-223:
# FOREIGN fires only when a foreign lineage exceeds FOREIGN_PCT_FLAG *and* the
# native lineage falls below NATIVE_PCT_FLAG. Because fibroblast almost never
# falls below 40, the native condition almost never fails. Nine clusters at
# res 0.6 carry a foreign lineage above 40%; exactly one (cluster 12,
# endothelial, fibroblast 28.1%) was flagged.
#
# So: keep the marker panels for the NEGATIVE half of the gate, where they work,
# and replace the POSITIVE half with a reference-based call that does not
# depend on collagen at all.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix); library(ggplot2)
})
for (p in c("SingleR", "celldex", "SummarizedExperiment", "BiocParallel")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("'%s' not installed. On a LOGIN node: Rscript ~/setup_singler_refs.R", p),
         call. = FALSE)
}

ts <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                paste0(..., collapse = "")))

# =============================== parameters ==================================
args   <- commandArgs(trailingOnly = TRUE)
SUBSET <- if (length(args) >= 1) args[1] else "stromal"
RES    <- if (length(args) >= 2) as.numeric(args[2]) else 0.6
source("/ix/tbruno/anp795/PDAC_Atlas/code/rebuild/io.R")   # read_obj(), save_obj(), dirs
stopifnot(SUBSET %in% SUBSETS)

DIR     <- V4_DIR
REF_DIR <- "/ix/tbruno/anp795/PDAC_Atlas/ref"   # shared, not per rebuild
IN      <- file.path(DIR, sprintf("subset_%s_harmony.qs2", SUBSET))
tag     <- sprintf("%s_res%s", SUBSET, RES)
F_CELL  <- file.path(DIR, sprintf("step4b_%s_singler_percell.tsv", tag))
F_GATE  <- file.path(DIR, sprintf("step4b_%s_gate.tsv", tag))
F_PDF   <- file.path(DIR, sprintf("step4b_%s_singler.pdf", tag))
F_PROV  <- file.path(DIR, sprintf("PROVENANCE_step4b_%s.txt", tag))

SINGLER_MIN_FRAC <- 0.50   # cluster's top lineage must reach this to count
FOREIGN_PCT_FLAG <- 40     # matches rebuild_04_markers.R
NTHREADS         <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
SEED             <- 42
set.seed(SEED)

# NATIVE: the SingleR lineages (map_label vocabulary below) that belong in the
# subset. The positive half of the gate asks that >= SINGLER_MIN_FRAC of a
# cluster's cells carry one of them, so an NK cluster passes in TNK and a
# plasma cluster in Bcell (v3 used one lineage and sent both to REASSIGN).
# stromal: HPCA cannot tell pericytes from fibroblasts (both mostly land on
# Tissue_stem_cells), so this gate is a stromal-LINEAGE check. Which stromal
# clusters are fibroblasts, and so get CAF scores, is decided in code/calls.
NATIVE <- list(
  TNK         = c("Tcell", "NK"),
  stromal     = c("fibroblast", "mural", "neural"),
  Bcell       = c("Bcell", "plasma"),
  myeloid     = "myeloid",
  endothelial = "endothelial"
)[[SUBSET]]
NATIVE_LABEL <- paste(NATIVE, collapse = "/")
# Marker panels whose PRESENCE disqualifies a cluster.
# myeloid leaves out fibroblast on purpose: macrophages carry collagen reads
# (84.5% fibroblast panel in the v3 stromal macrophage cluster), so the panel
# would fail real macrophage clusters.
DISQUALIFY <- list(
  TNK         = c("myeloid", "endothelial", "epithelial", "fibroblast"),
  stromal     = c("myeloid", "endothelial", "epithelial", "Tcell", "Bcell"),
  Bcell       = c("myeloid", "endothelial", "epithelial", "fibroblast"),
  myeloid     = c("Tcell", "Bcell", "endothelial", "epithelial"),
  endothelial = c("myeloid", "Tcell", "Bcell", "epithelial", "fibroblast")
)[[SUBSET]]

# ===================== marker panels (verbatim from rebuild_04_markers.R) =====
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
  acinar      = c("PRSS1", "CTRB1", "CPA1", "REG1A"),
  MHCII       = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD74"),
  neutrophil  = c("CSF3R", "FCGR3B", "CXCR2", "CXCR1"),
  mast        = c("TPSAB1", "TPSB2", "CPA3", "MS4A2"),
  pDC         = c("LILRA4", "CLEC4C", "LRRC26", "PTCRA"),
  schwann     = c("SOX10", "PLP1", "S100B", "MPZ", "CDH19")
)

# ===================== reference label -> our vocabulary ======================
# Applied to label.main from either reference. Order matters: first match wins.
#
# HPCA has no CAF profile. Tumour fibroblasts land on its other mesenchymal
# references -- Tissue_stem_cells (bone-marrow MSC-like stromal cells),
# Chondrocytes, Osteoblasts -- and Blueprint's Mesangial cells are the same
# lineage. The first version mapped those to 'stemlike'/'chondro', so 0 of 20
# stromal clusters could pass a gate that asks for 'fibroblast', and the two
# references "agreed" on only 22.8% of cells while Blueprint called 17/20
# clusters fibroblast. They are a lineage call here: fibroblast/mesenchymal.
# The original reference label is kept in singler_<ref>_raw for review.
map_label <- function(x) {
  x0 <- tolower(x)
  out <- rep(NA_character_, length(x0))
  rule <- list(
    fibroblast  = "fibroblast|^msc$|mesenchymal stem|tissue_stem|chondrocyte|osteoblast|mesangial",
    mural       = "smooth_muscle|smooth muscle|pericyte",
    endothelial = "endothelial",
    neural      = "neuron|neuroepithelial|astrocyte|schwann|glia|oligodendro",
    epithelial  = "epithelial|keratinocyte|hepatocyte",
    # myelocyte/CMP/GMP: HPCA has no mast-cell profile and calls mast cells
    # (and some neutrophils) Pro-Myelocyte, Myelocyte, CMP or GMP
    myeloid     = "macrophage|monocyte|dc$|^dc|dendritic|neutrophil|myeloid|eosinophil|basophil|mast|myelocyte|^cmp$|^gmp$|granulocyte",
    Tcell       = "t_cell|t-cell|t cells|cd4\\+|cd8\\+|tregs",
    NK          = "nk_cell|nk cell|nk cells",
    Bcell       = "b_cell|b-cell|b cells|pre-b|pro-b",
    plasma      = "plasma",
    RBC         = "erythro",
    adipocyte   = "adipocyte",
    stemlike    = "^ips|embryonic"
  )
  for (nm in names(rule)) {
    hit <- is.na(out) & grepl(rule[[nm]], x0)
    out[hit] <- nm
  }
  out[is.na(out)] <- "other"
  out
}

# ============================ 1. load ========================================
ts("subset = ", SUBSET, " | resolution = ", RES, " | threads = ", NTHREADS)
ts("loading ", IN)
sub <- read_obj(IN)
ts("loaded: ", format(nrow(sub), big.mark = ","), " genes x ",
   format(ncol(sub), big.mark = ","), " cells")

res_col <- res_col_of(sub, RES)
ts("assay: ", DefaultAssay(sub), " | cluster column: ", res_col)
Idents(sub) <- res_col
cl <- Idents(sub)
ts("clusters: ", nlevels(cl))

# SingleR wants log-normalised expression. The data layer is
# LogNormalize(scale.factor=1e4), carried unchanged from step 1
# (rebuild_02_subset.R lines 14-16). Correct input, no re-normalisation.
expr <- GetAssayData(sub, assay = DefaultAssay(sub), layer = "data")
ts("expression matrix: ", paste(dim(expr), collapse = " x "))

# ============================ 2. SingleR, per cell ===========================
bp <- BiocParallel::MulticoreParam(workers = NTHREADS)
refs <- c(HPCA = "celldex_HPCA", Blueprint = "celldex_BlueprintEncode")
preds <- list()

for (rn in names(refs)) {
  rf <- file.path(REF_DIR, paste0(refs[[rn]], ".rds"))
  if (!file.exists(rf))
    stop(sprintf("reference not cached: %s\n  On a LOGIN node: Rscript ~/setup_singler_refs.R", rf),
         call. = FALSE)
  ts("loading reference ", rn, " from ", rf)
  ref <- readRDS(rf)
  common <- intersect(rownames(expr), rownames(ref))
  ts("  genes shared with the atlas: ", format(length(common), big.mark = ","),
     " of ", format(nrow(ref), big.mark = ","), " in the reference")
  if (length(common) < 2000)
    warning(sprintf("only %d shared genes with %s -- labels will be unreliable",
                    length(common), rn))

  ts("  SingleR (label.main, per cell) ...")
  p <- SingleR::SingleR(test = expr, ref = ref, labels = ref$label.main,
                        de.method = "classic", BPPARAM = bp)
  preds[[rn]] <- p
  ts("  done. pruned (low confidence): ",
     sum(is.na(p$pruned.labels)), " of ", nrow(p),
     sprintf(" (%.1f%%)", 100 * mean(is.na(p$pruned.labels))))

  raw <- p$pruned.labels; raw[is.na(raw)] <- "unassigned"
  sub[[paste0("singler_", rn, "_raw")]] <- raw
  sub[[paste0("singler_", rn)]]         <- ifelse(raw == "unassigned", "unassigned", map_label(raw))
  rm(ref); gc(verbose = FALSE)
}

# report anything the mapping did not recognise, so it can be fixed rather than
# silently absorbed into "other"
for (rn in names(refs)) {
  raw <- sub[[paste0("singler_", rn, "_raw")]][, 1]
  mp  <- sub[[paste0("singler_", rn)]][, 1]
  unk <- sort(table(raw[mp == "other"]), decreasing = TRUE)
  if (length(unk)) {
    cat(sprintf("\n--- %s labels that fell through to 'other' ---\n", rn))
    print(head(unk, 20))
  }
}

agree <- sub[["singler_HPCA"]][, 1] == sub[["singler_Blueprint"]][, 1]
sub$singler_agree <- agree
ts("HPCA / Blueprint agreement on mapped lineage: ",
   sprintf("%.1f%% of cells", 100 * mean(agree)))

write.table(data.frame(cell = colnames(sub), cluster = as.character(cl),
                       dataset_id = as.character(sub$dataset_id),
                       HPCA_raw      = sub[["singler_HPCA_raw"]][, 1],
                       HPCA          = sub[["singler_HPCA"]][, 1],
                       Blueprint_raw = sub[["singler_Blueprint_raw"]][, 1],
                       Blueprint     = sub[["singler_Blueprint"]][, 1],
                       agree = agree),
            F_CELL, sep = "\t", quote = FALSE, row.names = FALSE)
ts("wrote ", F_CELL)

# ============================ 3. per-cluster aggregation =====================
idx <- split(seq_len(ncol(sub)), cl)

top_label <- function(v, ii) {
  t <- sort(table(v[ii]), decreasing = TRUE)
  c(label = names(t)[1], frac = as.numeric(t[1]) / length(ii))
}
hpca <- t(vapply(idx, function(ii) top_label(sub[["singler_HPCA"]][, 1], ii), character(2)))
blue <- t(vapply(idx, function(ii) top_label(sub[["singler_Blueprint"]][, 1], ii), character(2)))

# marker panels, same computation as rebuild_04_markers.R lines 151-157
pct_panel <- function(genes) {
  genes <- intersect(genes, rownames(sub))
  if (!length(genes)) return(setNames(rep(NA_real_, length(idx)), names(idx)))
  d <- expr[genes, , drop = FALSE]
  vapply(idx, function(ii) mean(Matrix::rowMeans(d[, ii, drop = FALSE] > 0)) * 100, numeric(1))
}
lin <- as.data.frame(lapply(lineage, pct_panel)); rownames(lin) <- names(idx)

# ============================ 4. THE GATE ====================================
# positive : SingleR (HPCA) top mapped lineage == NATIVE, at >= SINGLER_MIN_FRAC
# negative : no disqualifying marker panel above FOREIGN_PCT_FLAG
# The marker panels keep the negative half because that is where they work --
# CD68/LYZ/C1QA, PECAM1/VWF/CDH5 and EPCAM/KRT19 are not ambient in this tissue.
disq <- intersect(DISQUALIFY, names(lineage))
native_frac <- vapply(idx, function(ii) mean(sub[["singler_HPCA"]][ii, 1] %in% NATIVE), numeric(1))
show_panels <- intersect(unique(c(NATIVE, disq, "mural", "MHCII", "acinar")), names(lineage))
gate <- data.frame(
  cluster        = names(idx),
  n_cells        = vapply(idx, length, integer(1)),
  pct_of_subset  = round(100 * vapply(idx, length, integer(1)) / ncol(sub), 2),
  singler_HPCA   = hpca[, "label"],
  HPCA_frac      = round(as.numeric(hpca[, "frac"]), 3),
  singler_BP     = blue[, "label"],
  BP_frac        = round(as.numeric(blue[, "frac"]), 3),
  refs_agree     = hpca[, "label"] == blue[, "label"],
  HPCA_native_frac = round(native_frac, 3),
  round(lin[, show_panels, drop = FALSE], 1),
  row.names = NULL, check.names = FALSE
)

gate$positive  <- gate$HPCA_native_frac >= SINGLER_MIN_FRAC
gate$worst_foreign <- vapply(seq_len(nrow(gate)), function(i) {
  v <- as.numeric(lin[i, disq]); names(v) <- disq
  v <- v[!is.na(v)]
  if (!length(v) || max(v) <= FOREIGN_PCT_FLAG) "" else
    sprintf("%s=%.0f%%", names(which.max(v)), max(v))
}, character(1))
gate$negative  <- !nzchar(gate$worst_foreign)
gate$gate_pass <- gate$positive & gate$negative

gate$verdict <- ifelse(
  gate$gate_pass, sprintf("PASS - %s lineage", NATIVE_LABEL),
  ifelse(!gate$positive & !(gate$singler_HPCA %in% NATIVE),
         sprintf("REASSIGN - SingleR calls this %s (%.0f%%)", gate$singler_HPCA, 100 * gate$HPCA_frac),
         ifelse(!gate$negative,
                sprintf("FAIL - foreign marker signal %s", gate$worst_foreign),
                sprintf("FAIL - SingleR %s only %.0f%% (below %.0f%%)",
                        NATIVE_LABEL, 100 * gate$HPCA_native_frac, 100 * SINGLER_MIN_FRAC))))

write.table(gate, F_GATE, sep = "\t", quote = FALSE, row.names = FALSE)
ts("wrote ", F_GATE)

cat("\n================ LINEAGE GATE ================\n")
print(gate[, c("cluster","n_cells","pct_of_subset","singler_HPCA","HPCA_frac",
               "singler_BP","refs_agree","gate_pass","verdict")], row.names = FALSE)
cat(sprintf("\n  %d of %d clusters pass the gate (%s cells, %.1f%% of the subset)\n",
            sum(gate$gate_pass), nrow(gate),
            format(sum(gate$n_cells[gate$gate_pass]), big.mark = ","),
            100 * sum(gate$n_cells[gate$gate_pass]) / ncol(sub)))
cat("  Programme scores (myCAF / iCAF / apCAF / stellate) are interpretable\n")
cat("  ONLY on the passing clusters. On the rest they measure ambient RNA and\n")
cat("  MHC-II, which is how cluster 11 (macrophage) scored highest for apCAF.\n\n")
cat("  Where SingleR and the marker panels disagree, look before you decide.\n")
cat("  SingleR references contain no CAF states, so a fibroblast call is a\n")
cat("  LINEAGE call, never a myCAF/iCAF/apCAF call.\n\n")

# ============================ 5. plots =======================================
ts("writing ", F_PDF)
pdf(F_PDF, width = 12, height = 8.5)
for (rn in names(refs)) {
  print(DimPlot(sub, reduction = "umap", group.by = paste0("singler_", rn),
                label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.3) +
          ggtitle(sprintf("%s res %s - SingleR (%s), mapped lineage", SUBSET, RES, rn)))
}
print(DimPlot(sub, reduction = "umap", group.by = "singler_agree",
              raster = TRUE, pt.size = 0.3) +
        ggtitle(sprintf("%s res %s - HPCA / Blueprint agreement (%.1f%% of cells)",
                        SUBSET, RES, 100 * mean(agree))))

ct <- prop.table(table(cl, sub[["singler_HPCA"]][, 1]), 1)
dfl <- data.frame(cluster = factor(rep(rownames(ct), times = ncol(ct)), levels = rownames(ct)),
                  lineage = rep(colnames(ct), each = nrow(ct)),
                  frac    = as.numeric(ct))
print(ggplot(dfl, aes(cluster, frac, fill = lineage)) + geom_col() +
        geom_hline(yintercept = SINGLER_MIN_FRAC, linetype = 2) +
        theme_minimal(base_size = 11) +
        labs(title = sprintf("%s res %s - SingleR (HPCA) lineage composition per cluster", SUBSET, RES),
             subtitle = sprintf("dashed line = %.0f%% gate threshold", 100 * SINGLER_MIN_FRAC),
             y = "fraction of cluster", x = "cluster"))

gl <- data.frame(cluster = factor(gate$cluster, levels = gate$cluster),
                 pass = gate$gate_pass, n = gate$n_cells)
print(ggplot(gl, aes(cluster, n, fill = pass)) + geom_col() +
        scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#2166AC"),
                          name = "gate passed") +
        theme_minimal(base_size = 11) +
        labs(title = sprintf("%s res %s - clusters passing the lineage gate", SUBSET, RES),
             subtitle = sprintf("positive: SingleR %s >= %.0f%%   negative: no %s panel > %d%%",
                                NATIVE_LABEL, 100 * SINGLER_MIN_FRAC,
                                paste(disq, collapse = "/"), FOREIGN_PCT_FLAG),
             y = "cells", x = "cluster"))
dev.off()

# ============================ 6. save + provenance ===========================
save_obj(sub, IN)
ts("re-wrote ", IN, " with singler_* columns attached")

writeLines(c(
  sprintf("PROVENANCE - Step 4b SingleR lineage gate - '%s', resolution %s", SUBSET, RES),
  sprintf("generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "script   : rebuild_04b_singler.R",
  sprintf("input    : %s", IN),
  "",
  "WHY",
  "  Decision 9 (handoff v4 section 5) gates programme scores on the fibroblast",
  "  marker panel. Measured on step4_stromal_res0.6_decision.tsv, that panel is",
  "  detected in >=60% of cells in 17 of 20 clusters, including the macrophage",
  "  cluster (84.5%), the T cell cluster (76.7%), the ductal epithelium cluster",
  "  (96.8%) and the low-complexity artifact (80.6%). Ambient collagen is",
  "  atlas-wide: only GSE205013 and GSE212966 had SoupX applied by their authors.",
  "  In rebuild_04_markers.R the FOREIGN flag requires native < 40%, which almost",
  "  never happens, so 9 clusters carried a foreign lineage above 40% and only",
  "  one was flagged.",
  "",
  "METHOD",
  "  SingleR per cell against two celldex references, label.main, de.method",
  "  'classic', on the LogNormalize data layer (correct input; not re-normalised).",
  "    HumanPrimaryCellAtlas - non-immune stroma coverage",
  "    BlueprintEncode       - myeloid/lymphoid calibration",
  "  Raw labels are mapped to the rebuild_04_markers.R lineage vocabulary; labels",
  "  falling through to 'other' are printed in the log rather than hidden.",
  sprintf("  HPCA/Blueprint agreement: %.1f%% of cells", 100 * mean(agree)),
  "",
  "THE GATE",
  sprintf("  positive : SingleR (HPCA) lineage in {%s} for >= %.0f%% of the cluster's cells",
          paste(NATIVE, collapse = ", "), 100 * SINGLER_MIN_FRAC),
  sprintf("  assay    : %s, cluster column %s", DefaultAssay(sub), res_col),
  sprintf("  negative : no panel in {%s} detected above %d%%",
          paste(disq, collapse = ", "), FOREIGN_PCT_FLAG),
  sprintf("  result   : %d of %d clusters pass (%s cells, %.1f%% of subset)",
          sum(gate$gate_pass), nrow(gate),
          format(sum(gate$n_cells[gate$gate_pass]), big.mark = ","),
          100 * sum(gate$n_cells[gate$gate_pass]) / ncol(sub)),
  "",
  "LIMITS - READ BEFORE USING",
  "  - SingleR references contain no CAF states. A 'fibroblast' call is a",
  "    LINEAGE call. myCAF / iCAF / apCAF remain your judgement, gated on this.",
  "  - The references are bulk-sorted primary cells. Rare or tissue-specific",
  "    states (Schwann, mesothelium, chondrocyte-like) may map to 'other'.",
  "    'other' is not evidence of anything; read the raw labels.",
  "  - Ambient RNA affects SingleR too, just far less than a 5-gene collagen",
  "    panel. This reduces the problem; it does not remove it.",
  "  - Agreement between the two references is the confidence measure. A cluster",
  "    where they disagree needs looking at, whatever the gate says.",
  ""
), F_PROV)
ts("wrote ", F_PROV)
ts("done")
print(sessionInfo())
