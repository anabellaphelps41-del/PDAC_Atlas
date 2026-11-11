#!/usr/bin/env Rscript
# =============================================================================
# rebuild_05c_atlas_umap.R   --   one embedding of every kept cell (display only)
#
# Usage:  Rscript rebuild_05c_atlas_umap.R [all|embed|plot]        (default all)
#         embed: build the embedding; plot: redraw the figures from it.
#
# Input : <V4_DIR>/subset_<subset>_final.qs2, the five step-05 objects (read only)
#         <V4_DIR>/exclude_batchgenes.rds, batchgene_detection_<subset>.tsv (02b)
# Output: <V4_DIR>/atlas_v4_embedding.qs2   list: meta, pca, harmony, umap,
#                                           umap_unintegrated, hvg, masked, params, diag
#         <V4_DIR>/atlas_v4_umap.tsv.gz     per cell: both UMAPs + the reviewed labels
#         <V4_DIR>/figures_step5/atlas/{0_by_subset,0_all_clusters,1_compartments,2_lineages,3_batch,4_qc}.png
#         <V4_DIR>/figures_step5/by_cluster/<subset>.png   the subset's own UMAP, coloured by cluster
#         <V4_DIR>/PROVENANCE_step5c_atlas.txt
#
# Recipe = step 03 (rebuild_03_pipeline.R), run once on all five subsets:
#   decontX counts of the final objects, merged -> LogNormalize (1e4) ->
#   vst HVGs from a 5,000 pool, step 03 mask applied, top 2,000 -> ScaleData ->
#   PCA 50 -> Harmony on dataset_id (dims 1:50) -> UMAP on harmony 1:30, seed 42.
#   Also a UMAP of the uncorrected PCA and step 03's kNN diagnostics.
# Batch-gene mask: 02b writes one list per subset. Here the union of the five,
#   minus every gene 02b protected as a lineage marker in any subset
#   (GRADED_PROTECTED), so no compartment loses its markers.
# Display only: nothing is clustered and no label is changed. Labels come from
#   the reviewed subset calls (lineage_v4, celltype_v4, celltype_pool_v4).
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table)
  library(ggplot2); library(patchwork); library(shadowtext)
})
source("/ix/tbruno/anp795/PDAC_Atlas/code/rebuild/io.R")

MODE <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(MODE)) MODE <- "all"
if (!MODE %in% c("all", "embed", "plot")) stop("usage: rebuild_05c_atlas_umap.R [all|embed|plot]")
for (p in c("harmony", "RANN", "ragg", "pals", "ggrepel"))
  if (!requireNamespace(p, quietly = TRUE)) stop("package '", p, "' is not installed (install on a login node)")

DIR  <- V4_DIR
EMB  <- file.path(DIR, "atlas_v4_embedding.qs2")
TSV  <- file.path(DIR, "atlas_v4_umap.tsv.gz")
FIGD <- file.path(DIR, "figures_step5", "atlas")
PROV <- file.path(DIR, "PROVENANCE_step5c_atlas.txt")

N_HVG <- 2000; N_HVG_POOL <- 5000; N_PCS <- 50; N_DIMS <- 30; SEED <- 42; K_MIX <- 30
HARMONY_KEY <- "dataset_id"
NTHREADS <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
NAME <- c(TNK = "T/NK", stromal = "Stromal", Bcell = "B cell", myeloid = "Myeloid", endothelial = "Endothelial")

ts  <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))
fmt <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
show <- function(x) capture.output(print(x, row.names = FALSE))
set.seed(SEED)

# =============================== embed ========================================
if (MODE %in% c("all", "embed")) {
  t0 <- Sys.time()
  META_COLS <- c("dataset_id", "sample_key", "patient_key", "tissue_class", "cluster_v4", "subcluster_v4",
                 "lineage_v4", "celltype_v4", "celltype_pool_v4", "cell_flag", "exclude_from_de_lr",
                 "score_dissociation", "decontX_contamination", "nCount_RNA", "nFeature_RNA")
  mats <- list(); metas <- list()
  for (s in SUBSETS) {
    f <- file.path(DIR, sprintf("subset_%s_final.qs2", s))
    o <- read_obj(f)
    mats[[s]] <- LayerData(o, assay = "decontX", layer = "counts")
    md <- o@meta.data
    mt <- as.data.table(md[, intersect(META_COLS, names(md)), drop = FALSE], keep.rownames = "cell")
    for (cc in setdiff(META_COLS, names(mt))) mt[, (cc) := NA]
    mt[, subset := s]
    metas[[s]] <- mt
    ts(s, ": ", fmt(ncol(mats[[s]])), " cells, ", fmt(nrow(mats[[s]])), " genes (", f, ")")
    rm(o, md); invisible(gc())
  }
  genes <- Reduce(union, lapply(mats, rownames))
  if (!all(vapply(mats, function(m) identical(rownames(m), genes), TRUE))) {
    ts("gene sets differ between subsets; padding to the union of ", fmt(length(genes)), " genes with zeros")
    mats <- lapply(mats, function(m) {
      miss <- setdiff(genes, rownames(m))
      if (length(miss)) m <- rbind(m, Matrix(0, length(miss), ncol(m), sparse = TRUE, dimnames = list(miss, colnames(m))))
      m[genes, , drop = FALSE]
    })
  }
  M <- do.call(cbind, unname(mats)); rm(mats); invisible(gc())
  meta <- rbindlist(metas, use.names = TRUE, fill = TRUE); rm(metas)
  if (anyDuplicated(meta$cell)) stop("a cell is in more than one final object: ", meta$cell[anyDuplicated(meta$cell)])
  stopifnot(identical(meta$cell, colnames(M)))
  mdf <- as.data.frame(meta[, !"cell"]); rownames(mdf) <- meta$cell
  atl <- CreateSeuratObject(counts = M, assay = "decontX", meta.data = mdf, min.cells = 0, min.features = 0)
  rm(M, mdf); invisible(gc())
  ts("merged: ", fmt(nrow(atl)), " genes x ", fmt(ncol(atl)), " cells")
  ds_tab <- sort(table(atl$dataset_id), decreasing = TRUE)
  print(ds_tab)

  atl <- NormalizeData(atl, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

  # ---- HVG mask: step 03 section 2, with the batch-gene union ----------------
  excl_all <- readRDS(file.path(DIR, "exclude_batchgenes.rds"))
  protected <- unique(unlist(lapply(SUBSETS, function(s) {
    f <- file.path(DIR, sprintf("batchgene_detection_%s.tsv", s))
    if (!file.exists(f)) return(character(0))
    x <- fread(f, sep = "\t", select = c("gene", "class"), colClasses = "character")
    x[class == "GRADED_PROTECTED", gene]
  })))
  batch_excl <- setdiff(unique(unlist(excl_all)), protected)
  # keep in sync with rebuild_03_pipeline.R section 2
  mask_rx <- paste0("^IG[HKL][VJC]|^IGH[GAMDE]|^IGLL|^JCHAIN|",
                    "^TR[ABGD][VJC]|^MT-|^MTRNR|^RP[SL]|^HB[ABDGQZ][0-9]?$")
  ambient <- c("MALAT1", "NEAT1", "XIST",
               "PRSS1", "PRSS2", "PRSS3", "CTRB1", "CTRB2", "CELA2A", "CELA2B", "CELA3A", "CELA3B",
               "CPA1", "CPA2", "CPB1", "PNLIP", "PNLIPRP1", "CLPS", "AMY2A", "AMY2B", "CTRC", "CTRL",
               "REG1A", "REG1B", "REG3A", "REG3G", "SYCN", "CEL", "PLA2G1B", "SPINK1", "PPY")
  stress <- c("FOS", "JUN", "JUNB", "EGR1", "HSPA1A", "HSPA1B", "DNAJB1",
              "FOSB", "JUND", "ATF3", "IER2", "IER3", "SOCS3", "ZFP36",
              "HSPB1", "HSPH1", "HSPA8", "DUSP1", "KLF6", "PPP1R15A",
              "HSPA6", "DNAJA4", "DNAJB4", "BAG3", "ZFAND2A")
  stress_rx <- "^MT1[A-Z]$|^MT2A$"
  g <- rownames(atl)
  m_rx      <- grep(mask_rx, g, value = TRUE)
  m_ambient <- intersect(ambient, g)
  m_stress  <- union(intersect(stress, g), grep(stress_rx, g, value = TRUE))
  m_batch   <- intersect(batch_excl, g)
  masked    <- unique(c(m_rx, m_ambient, m_stress, m_batch))
  ts("mask: regex ", length(m_rx), ", ambient ", length(m_ambient), ", stress ", length(m_stress),
     ", batch ", length(m_batch), " (union of 5 subsets minus ", length(protected), " protected) -> ", length(masked), " genes")

  atl <- FindVariableFeatures(atl, selection.method = "vst", nfeatures = N_HVG_POOL, verbose = FALSE)
  hvg_pool <- VariableFeatures(atl)
  hvg <- head(setdiff(hvg_pool, masked), N_HVG)
  top_masked <- head(intersect(hvg_pool[1:100], masked), 40)
  ts("HVG pool ", length(hvg_pool), " -> ", sum(hvg_pool %in% masked), " masked -> ", length(hvg), " used")
  VariableFeatures(atl) <- hvg

  atl <- ScaleData(atl, features = hvg, verbose = FALSE)
  atl <- RunPCA(atl, features = hvg, npcs = N_PCS, seed.use = SEED, verbose = FALSE)
  ts("PCA done; UMAP of the uncorrected PCA")
  atl <- RunUMAP(atl, reduction = "pca", dims = 1:N_DIMS, reduction.name = "umap_unintegrated",
                 reduction.key = "umapun_", seed.use = SEED, verbose = FALSE)
  ts("RunHarmony(group.by.vars = '", HARMONY_KEY, "', dims 1:", N_PCS, ")")
  atl <- harmony::RunHarmony(atl, group.by.vars = HARMONY_KEY, reduction.use = "pca", dims.use = 1:N_PCS,
                             reduction.save = "harmony", ncores = NTHREADS, plot_convergence = FALSE, verbose = FALSE)
  ts("UMAP of harmony dims 1:", N_DIMS)
  atl <- RunUMAP(atl, reduction = "harmony", dims = 1:N_DIMS, reduction.name = "umap",
                 reduction.key = "umap_", seed.use = SEED, verbose = FALSE)

  # ---- kNN diagnostics (step 03), plus compartment / lineage neighbourhoods ---
  ts("kNN diagnostics (k = ", K_MIX, ")")
  knn_of  <- function(emb, k = K_MIX) RANN::nn2(emb, k = k + 1)$nn.idx[, -1, drop = FALSE]
  nb_frac <- function(nn, lab, same) {
    l <- as.character(lab); m <- matrix(l[nn], nrow = nrow(nn))
    mean(rowMeans(if (same) m == l else m != l))
  }
  mdx <- atl@meta.data
  nn_p <- knn_of(Embeddings(atl, "pca")[, 1:N_DIMS])
  nn_h <- knn_of(Embeddings(atl, "harmony")[, 1:N_DIMS])
  diag <- data.table(
    metric = c("neighbours from another dataset", "neighbours with the same tissue class",
               "neighbours from the same compartment", "neighbours with the same lineage"),
    before_harmony = c(nb_frac(nn_p, mdx$dataset_id, FALSE), nb_frac(nn_p, mdx$tissue_class, TRUE),
                       nb_frac(nn_p, mdx$subset, TRUE), nb_frac(nn_p, mdx$lineage_v4, TRUE)),
    after_harmony  = c(nb_frac(nn_h, mdx$dataset_id, FALSE), nb_frac(nn_h, mdx$tissue_class, TRUE),
                       nb_frac(nn_h, mdx$subset, TRUE), nb_frac(nn_h, mdx$lineage_v4, TRUE)),
    reference = c(1 - sum(prop.table(table(mdx$dataset_id))^2), sum(prop.table(table(mdx$tissue_class))^2),
                  sum(prop.table(table(mdx$subset))^2), sum(prop.table(table(mdx$lineage_v4))^2)))
  diag[, (2:4) := lapply(.SD, round, 3), .SDcols = 2:4]
  rm(nn_p, nn_h); invisible(gc())
  print(diag)

  cells <- colnames(atl)
  res <- list(
    meta = as.data.table(atl@meta.data, keep.rownames = "cell"),
    pca = Embeddings(atl, "pca")[cells, ], harmony = Embeddings(atl, "harmony")[cells, ],
    umap = Embeddings(atl, "umap")[cells, ], umap_unintegrated = Embeddings(atl, "umap_unintegrated")[cells, ],
    hvg = hvg, masked = list(regex = m_rx, ambient = m_ambient, stress = m_stress, batch = m_batch, protected = protected),
    params = list(N_HVG = N_HVG, N_HVG_POOL = N_HVG_POOL, N_PCS = N_PCS, N_DIMS = N_DIMS, SEED = SEED,
                  HARMONY_KEY = HARMONY_KEY, K_MIX = K_MIX),
    diag = diag, created = Sys.time())
  save_obj(res, EMB)
  ts("wrote ", EMB)
  out <- res$meta[, .(cell, subset, dataset_id, sample_key, tissue_class, lineage_v4, celltype_v4,
                      celltype_pool_v4, cell_flag, exclude_from_de_lr)]
  out[, `:=`(umap_1 = res$umap[, 1], umap_2 = res$umap[, 2],
             umap_unintegrated_1 = res$umap_unintegrated[, 1], umap_unintegrated_2 = res$umap_unintegrated[, 2])]
  fwrite(out, TSV, sep = "\t")
  ts("wrote ", TSV)

  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  writeLines(c(
    "PROVENANCE - Step 5c atlas-level embedding (display only)",
    sprintf("generated: %s (%.1f min)", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), mins),
    "script   : rebuild_05c_atlas_umap.R",
    sprintf("inputs   : %s", paste(sprintf("subset_%s_final.qs2", SUBSETS), collapse = ", ")),
    sprintf("cells    : %s from %d samples, %d datasets; %s genes", fmt(ncol(atl)), uniqueN(mdx$sample_key),
            uniqueN(mdx$dataset_id), fmt(nrow(atl))),
    "", "CELLS PER SUBSET", capture.output(print(table(mdx$subset))),
    "", "CELLS PER DATASET", capture.output(print(ds_tab)),
    "", "RECIPE (step 03, rebuild_03_pipeline.R)",
    "  decontX counts -> LogNormalize (1e4) -> vst HVGs, pool 5,000, mask, top 2,000 -> ScaleData",
    sprintf("  -> PCA %d -> Harmony on %s (dims 1:%d) -> UMAP on harmony 1:%d, seed %d", N_PCS, HARMONY_KEY, N_PCS, N_DIMS, SEED),
    "  also: UMAP on the uncorrected PCA (umap_unintegrated). Nothing is clustered.",
    "", "HVG MASK",
    sprintf("  regex families (IG/TR/MT/RP/HB) : %d", length(m_rx)),
    sprintf("  ambient pancreatic + nuisance   : %d", length(m_ambient)),
    sprintf("  dissociation stress (wide set)  : %d", length(m_stress)),
    sprintf("  batch genes                     : %d (union of the 5 per-subset 02b lists, minus %d GRADED_PROTECTED: %s)",
            length(m_batch), length(protected), paste(protected, collapse = ", ")),
    sprintf("  unique masked                   : %d", length(masked)),
    sprintf("  HVG pool %d -> %d masked -> %d used", length(hvg_pool), sum(hvg_pool %in% masked), length(hvg)),
    sprintf("  masked genes in the vst top 100: %s", paste(top_masked, collapse = ", ")),
    "", sprintf("kNN DIAGNOSTICS (k = %d, dims 1:%d; reference = expected under random neighbours)", K_MIX, N_DIMS),
    show(diag),
    "", "OUTPUTS", paste0("  ", c(EMB, TSV), collapse = "\n"),
    "  The labels are the reviewed subset calls; this embedding is for display and must not be",
    "  used to re-cluster or relabel cells."
  ), PROV)
  ts("wrote ", PROV)
  rm(atl); invisible(gc())
}

# =============================== plot =========================================
if (MODE %in% c("all", "plot")) {
  res <- read_obj(EMB)
  dir.create(FIGD, recursive = TRUE, showWarnings = FALSE)
  SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; ZERO <- "#f0efec"
  HI <- "#2a78d6"; BGPT <- "#dcdbd4"; ALLPT <- "#b5b4ad"
  SEQ <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
  DPI <- 300

  m <- res$meta
  d <- m[, .(cell, subset, dataset = dataset_id, sample = sample_key, tissue = tissue_class,
             lineage = lineage_v4, id = as.character(subcluster_v4), celltype = celltype_v4,
             excl = as.logical(exclude_from_de_lr) %in% TRUE,
             dissoc = as.numeric(score_dissociation), contam = as.numeric(decontX_contamination),
             nfeat = as.numeric(nFeature_RNA))]
  d[, `:=`(x = res$umap[cell, 1], y = res$umap[cell, 2],
           ux = res$umap_unintegrated[cell, 1], uy = res$umap_unintegrated[cell, 2])]
  N <- nrow(d)
  ts("plotting ", fmt(N), " cells")
  psize <- function(scale = 1) max(0.03, min(0.6, 0.5 * sqrt(1e4 / N))) * scale

  th_umap <- theme_void(base_size = 8) + theme(
    plot.background  = element_rect(fill = SURF, colour = NA),
    panel.background = element_rect(fill = SURF, colour = NA),
    plot.title    = element_text(colour = INK, size = 8, face = "bold", lineheight = 1.05, margin = margin(b = 1.5)),
    plot.subtitle = element_text(colour = INK2, size = 6.5, lineheight = 1.05, margin = margin(b = 2)),
    plot.margin = margin(4, 5, 4, 5), legend.position = "none")
  th_page <- theme(plot.background = element_rect(fill = SURF, colour = NA),
                   plot.title = element_text(colour = INK, size = 12, face = "bold"),
                   plot.caption = element_text(colour = MUTED, size = 7, hjust = 0))
  blank <- function() ggplot() + theme_void() + theme(plot.background = element_rect(fill = SURF, colour = NA))
  aspect <- function(xc, yc) diff(range(d[[yc]])) / diff(range(d[[xc]]))
  PW <- 3.2                                             # panel width, inches
  ph <- function(xc = "x", yc = "y") (PW - 0.2) * aspect(xc, yc) + 0.45
  save_png <- function(p, f, w, h) {
    ggsave(file.path(FIGD, f), p, width = w, height = h, dpi = DPI, device = ragg::agg_png, bg = SURF)
    ts("wrote ", f, sprintf(" (%.1f x %.1f in)", w, h)); f
  }
  hl_panel <- function(sel, title, subtitle, xc = "x", yc = "y") {
    a <- d[!sel, .(px = get(xc), py = get(yc))]; b <- d[sel, .(px = get(xc), py = get(yc))]
    hs <- min(0.45, psize(0.6) * max(1, sqrt(5000 / max(1, nrow(b)))))   # small groups drawn larger, or they vanish
    ggplot() +
      geom_point(data = a, aes(px, py), colour = BGPT, size = psize(0.6), shape = 16, stroke = 0) +
      geom_point(data = b, aes(px, py), colour = HI, size = hs, shape = 16, stroke = 0) +
      coord_equal() + labs(title = title, subtitle = subtitle) + th_umap
  }
  value_panel <- function(v, title, subtitle) {
    dd <- d[!is.na(get(v))]
    q <- quantile(dd[[v]], c(0.01, 0.99))
    dd[, val := pmin(pmax(get(v), q[1]), q[2])]
    setorder(dd, val)
    ggplot(dd, aes(x, y, colour = val)) +
      geom_point(size = psize(0.6), shape = 16, stroke = 0) +
      scale_colour_gradientn(colours = c(ZERO, SEQ), name = NULL, breaks = scales::breaks_pretty(3)) +
      coord_equal() + labs(title = title, subtitle = subtitle) + th_umap +
      theme(legend.position = "right", legend.key.width = unit(5, "pt"), legend.key.height = unit(16, "pt"),
            legend.text = element_text(colour = INK2, size = 6))
  }
  NOTE <- "Atlas-level embedding of the step-05 final cells (decontX counts; step 03 recipe; Harmony on dataset_id). Display only: labels are the reviewed subset calls."
  files <- character(0)

  # ---- 1. compartments ---------------------------------------------------------
  cn <- d[, .(n = .N, s = uniqueN(sample)), by = subset][order(-n)]
  big <- d[, .N, by = .(subset, lineage)][order(-N)][, .SD[1], by = subset]
  lab <- d[big[, .(subset, lineage)], on = .(subset, lineage)][, .(x = median(x), y = median(y)), by = subset]
  lab[, txt := NAME[subset]]
  p_all <- ggplot(d, aes(x, y)) +
    geom_point(colour = ALLPT, size = psize(0.6), shape = 16, stroke = 0) +
    geom_shadowtext(data = lab, aes(x, y, label = txt), colour = INK, bg.colour = SURF, bg.r = 0.15,
                    size = 3.1, fontface = "bold") +
    coord_equal() + th_umap +
    labs(title = "All kept cells", subtitle = sprintf("%s cells · %d samples · %d datasets", fmt(N), uniqueN(d$sample), uniqueN(d$dataset)))
  p_c <- c(list(p_all), lapply(seq_len(nrow(cn)), function(i)
    hl_panel(d$subset == cn$subset[i], NAME[[cn$subset[i]]],
             sprintf("%s cells (%.0f%%) · %d samples", fmt(cn$n[i]), 100 * cn$n[i] / N, cn$s[i]))))
  p1 <- wrap_plots(p_c, ncol = 3) +
    plot_annotation(title = sprintf("PDAC atlas v4: %s kept cells in five compartments", fmt(N)),
                    caption = NOTE, theme = th_page)
  files <- c(files, save_png(p1, "1_compartments.png", 3 * PW + 0.4, 2 * ph() + 0.9))

  # ---- 0. atlas coloured by subset (names on the plot as well) -------------------
  SUBCOL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")   # reference palette slots 1-5
  scol <- setNames(SUBCOL[seq_len(nrow(cn))], cn$subset)
  set.seed(SEED)
  p0 <- ggplot(d[sample(.N)], aes(x, y, colour = subset)) +
    geom_point(size = psize(2), shape = 16, stroke = 0) +
    geom_shadowtext(data = lab, aes(x, y, label = txt), inherit.aes = FALSE, colour = INK, bg.colour = SURF,
                    bg.r = 0.15, size = 4.2, fontface = "bold") +
    scale_colour_manual(values = scol, breaks = cn$subset, name = NULL,
                        labels = sprintf("%s · %s cells (%.0f%%)", NAME[cn$subset], fmt(cn$n), 100 * cn$n / N)) +
    guides(colour = guide_legend(override.aes = list(size = 3.2))) +
    coord_equal() + th_umap +
    theme(legend.position = "right", legend.text = element_text(colour = INK, size = 8.5)) +
    labs(title = sprintf("PDAC atlas v4: %s kept cells by subset", fmt(N)), caption = NOTE) +
    theme(plot.title = element_text(colour = INK, size = 12, face = "bold"),
          plot.caption = element_text(colour = MUTED, size = 7, hjust = 0))
  files <- c(files, save_png(p0, "0_by_subset.png", 10, 7.4 * aspect("x", "y") + 1.1))

  # ---- by cluster: each subset's own UMAP coloured by cluster --------------------
  # Up to 20 clusters per subset: Kelly's max-contrast colours (white, black and
  # grey removed) plus three. Cluster ids are printed on the plot, because 20 hues
  # cannot all be told apart by colour alone.
  CLCOL <- c(unname(pals::kelly(22))[-c(1, 2, 9)], "#1B9E9E", "#7F3C8D", "#11A579")
  id_order <- function(x) { x <- unique(x); x[order(suppressWarnings(as.integer(sub("[a-z]+$", "", x))), sub("^[0-9]+", "", x))] }
  CLD <- file.path(DIR, "figures_step5", "by_cluster")
  dir.create(CLD, recursive = TRUE, showWarnings = FALSE)
  rep_lab <- function(data, xc, yc) ggrepel::geom_text_repel(
    data = data, aes(.data[[xc]], .data[[yc]], label = id), inherit.aes = FALSE, colour = INK,
    bg.colour = SURF, bg.r = 0.15, size = 3, fontface = "bold", segment.colour = NA,
    box.padding = 0.05, point.padding = 0, max.overlaps = Inf, seed = SEED)
  for (s in SUBSETS) {
    o  <- read_obj(file.path(DIR, sprintf("subset_%s_final.qs2", s)))
    su <- Embeddings(o, "umap")[, 1:2]; rm(o); invisible(gc())
    ds <- d[subset == s]
    if (!all(ds$cell %in% rownames(su))) stop(s, ": cells of the embedding are missing from the final object")
    ds[, `:=`(sx = su[cell, 1], sy = su[cell, 2])]
    ids <- id_order(ds$id)
    if (length(ids) > length(CLCOL)) stop(s, ": ", length(ids), " clusters but only ", length(CLCOL), " colours")
    cc <- setNames(CLCOL[seq_along(ids)], ids)
    g <- ds[, .(n = .N, celltype = celltype[1], sx = median(sx), sy = median(sy), x = median(x), y = median(y)), by = id][match(ids, id)]
    leg <- sprintf("%s  %s (%s)", g$id, g$celltype, fmt(g$n))
    ds[, id := factor(id, levels = ids)]
    set.seed(SEED); ds <- ds[sample(.N)]
    pss <- 2 * max(0.03, min(0.45, 0.5 * sqrt(1e4 / nrow(ds))))        # larger dots so the clusters look filled
    p <- ggplot(ds, aes(sx, sy, colour = id)) +
      geom_point(size = pss, shape = 16, stroke = 0) + rep_lab(g, "sx", "sy") +
      scale_colour_manual(values = cc, breaks = ids, labels = leg, name = NULL, drop = FALSE) +
      guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2.8))) +
      coord_equal() + th_umap +
      theme(legend.position = "right", legend.text = element_text(colour = INK, size = 7.5),
            legend.key.size = unit(10, "pt"), legend.background = element_rect(fill = SURF, colour = NA),
            plot.title = element_text(colour = INK, size = 12, face = "bold"),
            plot.caption = element_text(colour = MUTED, size = 7, hjust = 0), plot.caption.position = "plot") +
      labs(title = sprintf("%s: %s kept cells in %d clusters", NAME[[s]], fmt(nrow(ds)), length(ids)),
           caption = "UMAP: v4 run-2 Harmony embedding of this subset (batch = dataset_id), where the clusters were made; kept cells after step 05, not re-embedded.")
    asp <- diff(range(ds$sy)) / diff(range(ds$sx))
    f <- file.path(CLD, sprintf("%s.png", s))
    hh <- max(8.2 * asp + 1.1, 0.2 * length(ids) + 1.8)
    ggsave(f, p, width = 14, height = hh, dpi = DPI, device = ragg::agg_png, bg = SURF)
    ts("wrote ", f, sprintf(" (14 x %.1f in)", hh))
  }

  # ---- 0. every cluster of every subset on the atlas UMAP -------------------------
  # 79 groups; ids are the subset letter (T, M, S, E, B) + cluster, printed on the plot
  # and in the key. Colours: clusters that border each other on this UMAP get the
  # most different colours. Bordering = at least 1% of either cluster's 15 nearest
  # UMAP neighbours fall in the other. A greedy pass (most-connected cluster first)
  # gives each cluster the palette colour farthest in CIE Lab from its coloured
  # neighbours; a colour repeats only between clusters that do not touch.
  PREFIX <- c(TNK = "T", myeloid = "M", stromal = "S", endothelial = "E", Bcell = "B")
  PAL <- c(unname(pals::kelly(22))[-c(1, 2, 9)],
           "#1B9E9E", "#7F3C8D", "#11A579", "#3969AC", "#E73F74", "#80BA5A", "#CF1C90")
  PLAB <- grDevices::convertColor(t(grDevices::col2rgb(PAL)) / 255, from = "sRGB", to = "Lab")
  gall <- d[, .(n = .N, celltype = celltype[1], x = median(x), y = median(y)), by = .(subset, id)]
  gall[, rank := match(id, id_order(id)), by = subset]
  gall[, gid := paste0(PREFIX[subset], id)]
  d[, gid := paste0(PREFIX[subset], id)]
  gi  <- match(d$gid, gall$gid)
  knn <- RANN::nn2(as.matrix(d[, .(x, y)]), k = 16)$nn.idx[, -1, drop = FALSE]
  prs <- data.table(a = rep(gi, ncol(knn)), b = gi[as.vector(knn)])[a != b, .N, by = .(a, b)]
  rm(knn); invisible(gc())
  prs[, w := N / (15 * gall$n[a])]
  adj <- prs[w >= 0.002, .(a, b)]
  # also neighbours: clusters whose centres lie close (small scattered clusters rarely share kNN)
  cd  <- as.matrix(dist(gall[, .(x, y)]))
  rad <- 0.06 * sqrt(diff(range(d$x))^2 + diff(range(d$y))^2)
  cl  <- which(cd < rad & upper.tri(cd), arr.ind = TRUE)
  adj <- unique(rbind(adj, data.table(a = cl[, 1], b = cl[, 2])))
  adj <- unique(rbind(adj, adj[, .(a = b, b = a)]))
  nbr <- split(adj$b, factor(adj$a, levels = seq_len(nrow(gall))))
  deg <- lengths(nbr)
  col_i <- rep(NA_integer_, nrow(gall)); used <- integer(length(PAL))
  for (i in order(-deg, -gall$n)) {
    js <- nbr[[i]]; js <- js[!is.na(col_i[js])]
    sc <- vapply(seq_along(PAL), function(k) {
      if (!length(js)) return(Inf)
      min(sqrt(rowSums((PLAB[col_i[js], , drop = FALSE] - matrix(PLAB[k, ], length(js), 3, byrow = TRUE))^2)))
    }, 0)
    # near-best colours (>= 80% of the best distance), preferring one not yet used in this
    # subset, then the least used overall, then the most distant
    cand  <- if (is.infinite(max(sc))) seq_along(PAL) else which(sc >= 0.8 * max(sc))
    fresh <- setdiff(cand, col_i[gall$subset == gall$subset[i]])
    if (length(fresh)) cand <- fresh
    cand  <- cand[used[cand] == min(used[cand])]
    best  <- cand[which.max(sc[cand])]
    col_i[i] <- best; used[best] <- used[best] + 1L
  }
  allcol <- setNames(PAL[col_i], gall$gid)
  de_adj <- sqrt(rowSums((PLAB[col_i[adj$a], , drop = FALSE] - PLAB[col_i[adj$b], , drop = FALSE])^2))
  rep_in_sub <- sum(gall[, .(r = .N - uniqueN(col_i[.I])), by = subset]$r)
  ts("all-clusters colours: ", nrow(adj) / 2, " neighbouring pairs, ", length(unique(col_i)), " of ", length(PAL),
     " colours used; smallest Lab distance between neighbours ", round(min(de_adj), 1),
     " (median ", round(median(de_adj), 1), "); colour repeats within a subset: ", rep_in_sub)
  set.seed(SEED)
  p_all <- ggplot(d[sample(.N)], aes(x, y, colour = gid)) +
    geom_point(size = psize(2), shape = 16, stroke = 0) +
    ggrepel::geom_text_repel(data = gall, aes(x, y, label = gid), inherit.aes = FALSE, colour = INK,
                             bg.colour = SURF, bg.r = 0.15, size = 2.5, fontface = "bold",
                             segment.colour = INK2, segment.size = 0.2, min.segment.length = 0.4,
                             box.padding = 0.12, point.padding = 0, max.overlaps = Inf, seed = SEED) +
    scale_colour_manual(values = allcol, guide = "none") +
    coord_equal() + th_umap
  kcols <- list(c("TNK", "Bcell", "endothelial"), c("myeloid", "stromal"))
  krows <- list()
  for (ci in seq_along(kcols)) {
    yy <- 0
    for (s in kcols[[ci]]) {
      gg <- gall[subset == s][order(rank)]
      krows[[length(krows) + 1]] <- data.table(col = ci, y = yy, head = TRUE, gid = NA_character_,
                                               txt = sprintf("%s · %s cells", NAME[[s]], fmt(sum(gg$n))))
      krows[[length(krows) + 1]] <- data.table(col = ci, y = yy - seq_len(nrow(gg)), head = FALSE, gid = gg$gid,
                                               txt = sprintf("%s  %s (%s)", gg$gid, gg$celltype, fmt(gg$n)))
      yy <- yy - nrow(gg) - 1.6
    }
  }
  key <- rbindlist(krows)
  p_key <- ggplot() +
    geom_point(data = key[head == FALSE], aes((col - 1), y, colour = gid), size = 2.3, shape = 16) +
    geom_text(data = key[head == FALSE], aes((col - 1) + 0.025, y, label = txt), hjust = 0, size = 2.45, colour = INK) +
    geom_text(data = key[head == TRUE], aes((col - 1) - 0.012, y, label = txt), hjust = 0, size = 2.9,
              fontface = "bold", colour = INK) +
    scale_colour_manual(values = allcol, guide = "none") +
    scale_x_continuous(limits = c(-0.02, 2), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(add = 0.8)) +
    coord_cartesian(clip = "off") + theme_void() +
    theme(plot.background = element_rect(fill = SURF, colour = NA))
  UW <- 11
  p_ac <- (p_all | p_key) + plot_layout(widths = c(UW, 9.4)) +
    plot_annotation(title = sprintf("PDAC atlas v4: all %d clusters of the five subsets on the atlas UMAP", nrow(gall)),
                    caption = NOTE, theme = th_page)
  files <- c(files, save_png(p_ac, "0_all_clusters.png", UW + 9.8,
                             max(UW * aspect("x", "y") + 1.0, 0.27 * max(-min(key$y)) + 1.4)))

  # ---- 2. lineages -------------------------------------------------------------
  ln <- d[, .(n = .N, subset = subset[1]), by = lineage]
  ln[, so := match(subset, cn$subset)]; setorder(ln, so, -n)
  p_l <- lapply(seq_len(nrow(ln)), function(i)
    hl_panel(d$lineage == ln$lineage[i], ln$lineage[i],
             sprintf("%s · %s cells (%.1f%%)", NAME[[ln$subset[i]]], fmt(ln$n[i]), 100 * ln$n[i] / N)))
  nc <- 5
  p2 <- wrap_plots(p_l, ncol = nc) +
    plot_annotation(title = "PDAC atlas v4: reviewed lineages on the atlas UMAP", caption = NOTE, theme = th_page)
  files <- c(files, save_png(p2, "2_lineages.png", nc * PW + 0.4, ceiling(nrow(ln) / nc) * ph() + 0.9))

  # ---- 3. batch: datasets before / after Harmony, tissue class -----------------
  dn <- d[, .(n = .N, s = uniqueN(sample)), by = dataset][order(-n)]
  tn <- d[, .(n = .N, s = uniqueN(sample)), by = tissue][order(-n)]
  nb <- max(nrow(dn), nrow(tn))
  row_b <- lapply(seq_len(nrow(dn)), function(i)
    hl_panel(d$dataset == dn$dataset[i], dn$dataset[i], "before Harmony (uncorrected PCA)", "ux", "uy"))
  row_a <- lapply(seq_len(nrow(dn)), function(i)
    hl_panel(d$dataset == dn$dataset[i], dn$dataset[i],
             sprintf("after Harmony · %s cells (%.0f%%) · %d samples", fmt(dn$n[i]), 100 * dn$n[i] / N, dn$s[i])))
  row_t <- lapply(seq_len(nrow(tn)), function(i)
    hl_panel(d$tissue == tn$tissue[i], gsub("_", " ", tn$tissue[i]),
             sprintf("%s cells (%.0f%%) · %d samples", fmt(tn$n[i]), 100 * tn$n[i] / N, tn$s[i])))
  pad <- function(l) c(l, lapply(seq_len(nb - length(l)), function(i) blank()))
  p3 <- wrap_plots(c(pad(row_b), pad(row_a), pad(row_t)), ncol = nb,
                   heights = c(ph("ux", "uy"), ph(), ph())) +
    plot_annotation(title = "PDAC atlas v4: datasets before and after Harmony, and tissue class",
                    caption = NOTE, theme = th_page)
  files <- c(files, save_png(p3, "3_batch.png", nb * PW + 0.4, ph("ux", "uy") + 2 * ph() + 0.95))

  # ---- 4. QC: DE/LR exclusion, dissociation, contamination, depth ---------------
  p_q <- list(hl_panel(d$excl, "Kept out of DE / ligand-receptor",
                       sprintf("%s cells (%.1f%%)", fmt(sum(d$excl)), 100 * mean(d$excl))))
  if (any(!is.na(d$dissoc))) p_q <- c(p_q, list(value_panel("dissoc", "Dissociation score", "van den Brink core set; clipped at 1st / 99th percentile")))
  if (any(!is.na(d$contam))) p_q <- c(p_q, list(value_panel("contam", "DecontX contamination", "estimated ambient fraction per cell; clipped")))
  if (any(!is.na(d$nfeat))) { d[, lnf := log10(nfeat)]; p_q <- c(p_q, list(value_panel("lnf", "Genes detected (log10)", "RNA assay nFeature; clipped"))) }
  p4 <- wrap_plots(p_q, ncol = 2) +
    plot_annotation(title = "PDAC atlas v4: quality covariates on the atlas UMAP", caption = NOTE, theme = th_page)
  files <- c(files, save_png(p4, "4_qc.png", 2 * (PW + 0.5) + 0.4, ceiling(length(p_q) / 2) * ph() + 0.9))

  ts("figures: ", paste(files, collapse = ", "))
}
ts("done")
