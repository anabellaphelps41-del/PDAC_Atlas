#!/usr/bin/env Rscript
# =============================================================================
# rebuild_03_pipeline.R   --   Step 3 of the PDAC atlas rebuild
#
# Per-subset: HVG (masked) -> ScaleData -> PCA -> Harmony(dataset_id)
#             -> FindNeighbors -> FindClusters resolution sweep -> UMAP
#
# Usage:  Rscript rebuild_03_pipeline.R <subset>
#         <subset> in {TNK, stromal, Bcell, myeloid, endothelial}
#
# Runs on the default assay: decontX in v4 (rebuild_02_subset.R), RNA in v3.
# Cluster columns are named <assay>_snn_res.<r>.
#
# BATCH KEY: Harmony runs on PDAC_HARMONY_KEY (default dataset_id) and, for
# comparison only, once more on the other key (sample_key). The comparison
# reports batch mixing, how often a cell's neighbours share its tissue_class
# (a drop means real tumour/normal differences are being corrected away), and
# the adjusted Rand index between the two res-0.6 clusterings.
#
# Input : <DIR>/subset_<subset>.qs2
#         (built by rebuild_02_subset.R from atlas_harmonized_counts_v3.qs2)
# Output: subset_<subset>_harmony.qs2
#         step3_<subset>_qc.pdf
#         step3_<subset>_clustercomp.tsv
#         PROVENANCE_step3_<subset>.txt
#
# The HVG mask governs FEATURE SELECTION ONLY. Masked genes remain in the
# object and remain plottable as markers. This is deliberate -- IG/TR genes
# and ambient pancreatic enzymes must not drive the embedding, but JCHAIN,
# TRAC etc. are still legitimate things to plot at Step 4.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

ts <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                paste0(..., collapse = "")))

# --- fail fast on missing deps (compute nodes have no outbound internet) -----
for (p in c("harmony", "RANN")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop(sprintf(paste0("Package '%s' is not installed.\n",
                        "  Install it on a LOGIN node (compute nodes have no internet):\n",
                        "    module load r/4.6.0\n",
                        "    R -e 'install.packages(\"%s\", repos=\"https://cloud.r-project.org\")'"),
                 p, p), call. = FALSE)
  }
}

# =============================== parameters ==================================
args   <- commandArgs(trailingOnly = TRUE)
SUBSET <- if (length(args) >= 1) args[1] else "stromal"
source("/ix/tbruno/anp795/PDAC_Atlas/code/rebuild/io.R")   # read_obj(), save_obj(), dirs
stopifnot(SUBSET %in% SUBSETS)

DIR  <- V4_DIR
IN   <- file.path(DIR, sprintf("subset_%s.qs2", SUBSET))
OUT  <- file.path(DIR, sprintf("subset_%s_harmony.qs2", SUBSET))
PDFF <- file.path(DIR, sprintf("step3_%s_qc.pdf", SUBSET))
COMP <- file.path(DIR, sprintf("step3_%s_clustercomp.tsv", SUBSET))
PROV <- file.path(DIR, sprintf("PROVENANCE_step3_%s.txt", SUBSET))

MIN_CELLS_PER_BATCH <- 150     # Harmony fits a centroid per batch; below this it is noise
N_HVG               <- 2000    # HVGs kept AFTER masking
N_HVG_POOL          <- 5000    # vst pool to draw from, so masking does not starve the set
N_PCS               <- 50      # PCs computed
N_DIMS              <- 30      # PCs handed to Harmony/neighbours/UMAP
RESOLUTIONS         <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6)
SEED                <- 42
NTHREADS            <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))  # Harmony BLAS threads
K_MIX               <- 30      # k for the batch-mixing diagnostic
HARMONY_KEY         <- Sys.getenv("PDAC_HARMONY_KEY", "dataset_id")
ALT_KEY             <- if (HARMONY_KEY == "dataset_id") "sample_key" else "dataset_id"
stopifnot(HARMONY_KEY %in% c("dataset_id", "sample_key"))

set.seed(SEED)

ts("subset = ", SUBSET)
ts("params: N_HVG=", N_HVG, " pool=", N_HVG_POOL, " nPCs=", N_PCS,
   " dims=", N_DIMS, " seed=", SEED)

# ============================ 1. load + drop tiny batches ====================
ts("loading ", IN)
sub <- read_obj(IN)
ts("loaded: ", format(nrow(sub), big.mark = ","), " genes x ",
   format(ncol(sub), big.mark = ","), " cells")

stopifnot(all(c("dataset_id", "sample_key", "tissue_class") %in% colnames(sub@meta.data)))
ASSAY      <- DefaultAssay(sub)
RES_PREFIX <- paste0(ASSAY, "_snn_res.")
ts("default assay: ", ASSAY, " | Harmony key: ", HARMONY_KEY, " (comparison: ", ALT_KEY, ")")

tab     <- sort(table(sub$dataset_id), decreasing = TRUE)
ts("cells per dataset_id:")
print(tab)

keep_ds <- names(tab)[tab >= MIN_CELLS_PER_BATCH]
drop_ds <- setdiff(names(tab), keep_ds)

if (length(drop_ds) > 0) {
  n_drop <- sum(tab[drop_ds])
  ts("dropping ", length(drop_ds), " dataset(s) below ", MIN_CELLS_PER_BATCH,
     " cells: ", paste(sprintf("%s=%d", drop_ds, tab[drop_ds]), collapse = "  "))
  ts("cells dropped: ", n_drop, " (", sprintf("%.3f%%", 100 * n_drop / ncol(sub)), " of subset)")
  cells_keep <- colnames(sub)[sub$dataset_id %in% keep_ds]
  sub <- subset(sub, cells = cells_keep)
  sub$dataset_id <- droplevels(factor(sub$dataset_id))
  ts("after drop: ", format(ncol(sub), big.mark = ","), " cells, ",
     length(keep_ds), " batches")
} else {
  ts("no dataset below the ", MIN_CELLS_PER_BATCH, "-cell floor; nothing dropped")
  sub$dataset_id <- droplevels(factor(sub$dataset_id))
}

# ============================ 2. build the HVG mask ==========================
# 2a. structural/graded batch genes from rebuild_02b_batchgenes.R
excl_all <- readRDS(file.path(DIR, "exclude_batchgenes.rds"))
if (!SUBSET %in% names(excl_all)) {
  stop(sprintf("subset '%s' not found in exclude_batchgenes.rds (names: %s)",
               SUBSET, paste(names(excl_all), collapse = ", ")))
}
batch_excl <- as.character(excl_all[[SUBSET]])

# 2b. regex families
#     ^TR[ABGD][VJC] is the one people forget and it matters most here: TCR
#     V/J segments are clonally distributed, so left in they produce
#     patient-specific clusters that look like real T cell subtypes.
mask_rx <- paste0("^IG[HKL][VJC]|^IGH[GAMDE]|^IGLL|^JCHAIN|",
                  "^TR[ABGD][VJC]|^MT-|^MTRNR|^RP[SL]|^HB[ABDGQZ][0-9]?$")

# 2c. ambient pancreatic enzymes + the usual high-expression nuisances
ambient <- c("MALAT1", "NEAT1", "XIST",
             "PRSS1", "PRSS2", "PRSS3", "CTRB1", "CTRB2",
             "CELA2A", "CELA2B", "CELA3A", "CELA3B",
             "CPA1", "CPA2", "CPB1", "PNLIP", "PNLIPRP1", "CLPS",
             "AMY2A", "AMY2B", "CTRC", "CTRL",
             "REG1A", "REG1B", "REG3A", "REG3G",
             "SYCN", "CEL", "PLA2G1B", "SPINK1", "PPY")

# 2d. dissociation-stress genes. Six datasets, six dissociation protocols;
#     this is a plausible batch axis and no Step 4 call depends on these.
#     Still plottable -- if a cluster turns out to be stress-defined you will
#     see it in the marker step.
stress <- c("FOS", "JUN", "JUNB", "EGR1", "HSPA1A", "HSPA1B", "DNAJB1")
# Wider set, enabled 2026-11-09 (v4 run 2) by the user's decision. With only
# the 7 genes above masked, HSPA6, DNAJA4, HSPB1 and the metallothioneins still
# shaped v4 run-1 clusters: TNK c8, myeloid c2 and B c4 dissolved into
# ordinary cell types once these were masked too (cluster stability 0.83 ->
# 0.26, 0.92 -> 0.55, 0.77 -> 0.34). Sources: the van den Brink et al. 2017
# genes in the old comment here, the heat-shock genes that defined those
# clusters, and the metallothioneins (MT1*, MT2A; TNK c13 was MT-defined).
stress <- c(stress, "FOSB", "JUND", "ATF3", "IER2", "IER3", "SOCS3", "ZFP36",
            "HSPB1", "HSPH1", "HSPA8", "DUSP1", "KLF6", "PPP1R15A",
            "HSPA6", "DNAJA4", "DNAJB4", "BAG3", "ZFAND2A")
stress_rx <- "^MT1[A-Z]$|^MT2A$"

g          <- rownames(sub)
m_rx       <- grep(mask_rx, g, value = TRUE)
m_ambient  <- intersect(ambient, g)
m_stress   <- union(intersect(stress, g), grep(stress_rx, g, value = TRUE))
m_batch    <- intersect(batch_excl, g)
masked     <- unique(c(m_rx, m_ambient, m_stress, m_batch))

ts("HVG mask composition (genes present in this subset):")
cat(sprintf("    regex families (IG/TR/MT/RP/HB) : %5d\n", length(m_rx)))
cat(sprintf("    ambient pancreatic + nuisance   : %5d\n", length(m_ambient)))
cat(sprintf("    dissociation stress             : %5d\n", length(m_stress)))
cat(sprintf("    batch genes (structural+graded) : %5d\n", length(m_batch)))
cat(sprintf("    ------------------------------------------\n"))
cat(sprintf("    unique masked total             : %5d  (%.1f%% of %d genes)\n",
            length(masked), 100 * length(masked) / length(g), length(g)))

# ============================ 3. HVG selection ===============================
ts("FindVariableFeatures(vst, pool = ", N_HVG_POOL, ")")
sub <- FindVariableFeatures(sub, selection.method = "vst",
                            nfeatures = N_HVG_POOL, verbose = FALSE)
hvg_pool <- VariableFeatures(sub)            # vst-ranked, best first
hvg      <- head(setdiff(hvg_pool, masked), N_HVG)   # setdiff preserves rank

n_removed <- sum(hvg_pool %in% masked)
ts("pool ", length(hvg_pool), " -> masked out ", n_removed,
   " (", sprintf("%.1f%%", 100 * n_removed / length(hvg_pool)), ") -> kept top ", length(hvg))

if (length(hvg) < N_HVG) {
  warning(sprintf("only %d HVGs survived masking (wanted %d); raise N_HVG_POOL",
                  length(hvg), N_HVG))
}

ts("masked genes that were in the vst top 100 (these are the ones that would ",
   "have driven the embedding):")
print(head(intersect(hvg_pool[1:100], masked), 40))

VariableFeatures(sub) <- hvg

# ============================ 4. scale / PCA =================================
ts("ScaleData on ", length(hvg), " features")
sub <- ScaleData(sub, features = hvg, verbose = FALSE)

ts("RunPCA(npcs = ", N_PCS, ")")
sub <- RunPCA(sub, features = hvg, npcs = N_PCS, seed.use = SEED, verbose = FALSE)

ts("PC stdev (1:", N_PCS, "):")
print(round(Stdev(sub, reduction = "pca"), 3))

# unintegrated UMAP, kept for the before/after comparison
ts("RunUMAP on raw PCA (unintegrated reference embedding)")
sub <- RunUMAP(sub, reduction = "pca", dims = 1:N_DIMS,
               reduction.name = "umap_unintegrated", reduction.key = "umapun_",
               seed.use = SEED, verbose = FALSE)

# ============================ 5. Harmony =====================================
ts("RunHarmony(group.by.vars = '", HARMONY_KEY, "', dims 1:", N_PCS, ")")
sub <- harmony::RunHarmony(sub,
                           group.by.vars = HARMONY_KEY,
                           reduction.use = "pca",
                           dims.use      = 1:N_PCS,
                           reduction.save = "harmony",
                           ncores        = NTHREADS,
                           plot_convergence = FALSE,
                           verbose = TRUE)
ts("harmony embedding: ", paste(dim(Embeddings(sub, "harmony")), collapse = " x "))

ts("RunHarmony(group.by.vars = '", ALT_KEY, "') for comparison only")
sub <- harmony::RunHarmony(sub, group.by.vars = ALT_KEY, reduction.use = "pca",
                           dims.use = 1:N_PCS, reduction.save = "harmony_alt",
                           ncores = NTHREADS, plot_convergence = FALSE, verbose = FALSE)

# ---- diagnostics on the kNN graph of each embedding
#   mixing : fraction of a cell's neighbours from ANOTHER dataset
#   tissue : fraction of a cell's neighbours with the SAME tissue_class
knn_of  <- function(emb, k = K_MIX) RANN::nn2(emb, k = k + 1)$nn.idx[, -1, drop = FALSE]
nb_frac <- function(nn, lab, same) {
  l <- as.character(lab); m <- matrix(l[nn], nrow = nrow(nn))
  mean(rowMeans(if (same) m == l else m != l))
}
nn_pca  <- knn_of(Embeddings(sub, "pca")[, 1:N_DIMS])
nn_harm <- knn_of(Embeddings(sub, "harmony")[, 1:N_DIMS])
nn_alt  <- knn_of(Embeddings(sub, "harmony_alt")[, 1:N_DIMS])
p_ds     <- as.numeric(prop.table(table(sub$dataset_id)))
mix_ideal<- 1 - sum(p_ds^2)      # expected score under perfect mixing
mix_pca  <- nb_frac(nn_pca,  sub$dataset_id, same = FALSE)
mix_harm <- nb_frac(nn_harm, sub$dataset_id, same = FALSE)
mix_alt  <- nb_frac(nn_alt,  sub$dataset_id, same = FALSE)
tis_rand <- sum(as.numeric(prop.table(table(sub$tissue_class)))^2)   # random neighbours
tis_pca  <- nb_frac(nn_pca,  sub$tissue_class, same = TRUE)
tis_harm <- nb_frac(nn_harm, sub$tissue_class, same = TRUE)
tis_alt  <- nb_frac(nn_alt,  sub$tissue_class, same = TRUE)
rm(nn_pca, nn_harm, nn_alt)

ts("kNN batch mixing (k = ", K_MIX, "), fraction of neighbours from another dataset:")
cat(sprintf("    before Harmony (PCA)  : %.3f\n", mix_pca))
cat(sprintf("    after  Harmony        : %.3f\n", mix_harm))
cat(sprintf("    perfect-mixing target : %.3f\n", mix_ideal))
cat(sprintf("    recovered             : %.1f%% of the gap\n",
            100 * (mix_harm - mix_pca) / max(mix_ideal - mix_pca, 1e-9)))
cat("    READ: if 'after' is not well above 'before', Harmony under-corrected.\n")
cat("          Options then: more dims, or switch to scVI/scANVI (GPU).\n")
cat("          If 'after' overshoots the target, you may be over-correcting\n")
cat("          and erasing real biology -- check that known subtypes survive.\n")
cat(sprintf("    with Harmony on %-11s: %.3f\n", ALT_KEY, mix_alt))
ts("tissue_class neighbourhood (fraction of neighbours with the same tissue_class):")
cat(sprintf("    PCA %.3f | Harmony %s %.3f | Harmony %s %.3f | random %.3f\n",
            tis_pca, HARMONY_KEY, tis_harm, ALT_KEY, tis_alt, tis_rand))
cat("    READ: the further an embedding falls from PCA toward random, the more\n")
cat("          tumour vs normal difference it has removed.\n")

# ============================ 6. cluster sweep + UMAP ========================
ts("FindNeighbors on harmony dims 1:", N_DIMS)
sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:N_DIMS, verbose = FALSE)

for (r in RESOLUTIONS) {
  sub <- FindClusters(sub, resolution = r, algorithm = 1,
                      random.seed = SEED, verbose = FALSE)
  ts("  resolution ", r, " -> ",
     nlevels(sub[[paste0(RES_PREFIX, r)]][, 1]), " clusters")
}

# the comparison embedding, clustered at 0.6 only
ari <- function(a, b) {
  tab <- table(a, b); c2 <- function(x) sum(choose(x, 2))
  e <- c2(rowSums(tab)) * c2(colSums(tab)) / choose(sum(tab), 2)
  (c2(tab) - e) / ((c2(rowSums(tab)) + c2(colSums(tab))) / 2 - e)
}
main_idents <- Idents(sub); main_sc <- sub$seurat_clusters   # FindClusters below resets both
sub <- FindNeighbors(sub, reduction = "harmony_alt", dims = 1:N_DIMS,
                     graph.name = c("alt_nn", "alt_snn"), verbose = FALSE)
sub <- FindClusters(sub, graph.name = "alt_snn", resolution = 0.6, algorithm = 1,
                    cluster.name = "alt_res.0.6", random.seed = SEED, verbose = FALSE)
Idents(sub) <- main_idents; sub$seurat_clusters <- main_sc
ari_06 <- ari(sub[[paste0(RES_PREFIX, 0.6)]][, 1], sub$alt_res.0.6)
ts("res 0.6: ", nlevels(sub[[paste0(RES_PREFIX, 0.6)]][, 1]), " clusters (", HARMONY_KEY, ") vs ",
   nlevels(sub$alt_res.0.6), " (", ALT_KEY, "); adjusted Rand index ", sprintf("%.3f", ari_06))

ts("RunUMAP on harmony dims 1:", N_DIMS)
sub <- RunUMAP(sub, reduction = "harmony", dims = 1:N_DIMS,
               reduction.name = "umap", reduction.key = "umap_",
               seed.use = SEED, verbose = FALSE)
sub <- RunUMAP(sub, reduction = "harmony_alt", dims = 1:N_DIMS,
               reduction.name = "umap_alt", reduction.key = "umapalt_",
               seed.use = SEED, verbose = FALSE)

# ============================ 7. cluster composition =========================
# The question this answers: did any cluster survive Harmony as a dataset
# cluster? A cluster that is >80% one dataset is a batch artefact until proven
# otherwise -- unless that dataset is genuinely a distinct biological sample.
comp_rows <- list()
ts("per-resolution cluster purity (clusters dominated by a single dataset):")
for (r in RESOLUTIONS) {
  col  <- paste0(RES_PREFIX, r)
  cl   <- sub[[col]][, 1]
  ct   <- table(cl, sub$dataset_id)
  frac <- prop.table(ct, 1)
  mx   <- apply(frac, 1, max)
  dom  <- colnames(frac)[apply(frac, 1, which.max)]
  n    <- as.numeric(table(cl))
  flag <- mx > 0.80
  cat(sprintf("    res %-4s : %2d clusters | %d flagged >80%% single-dataset | %s cells in them\n",
              r, nlevels(cl), sum(flag), format(sum(n[flag]), big.mark = ",")))
  if (any(flag)) {
    cat(sprintf("        cluster %-3s n=%-7s %.0f%% %s\n",
                rownames(frac)[flag], format(n[flag], big.mark = ","),
                100 * mx[flag], dom[flag]))
  }
  comp_rows[[as.character(r)]] <- data.frame(
    resolution = r, cluster = rownames(frac), n_cells = n,
    max_dataset_frac = round(as.numeric(mx), 4), dominant_dataset = dom,
    round(as.data.frame.matrix(frac), 4), check.names = FALSE
  )
}
comp <- do.call(rbind, comp_rows)
write.table(comp, COMP, sep = "\t", quote = FALSE, row.names = FALSE)
ts("wrote ", COMP)

# ============================ 8. QC plots ====================================
marker_panels <- list(
  stromal = c("COL1A1", "DCN", "LUM", "ACTA2", "TAGLN", "POSTN",
              "IL6", "CXCL12", "CFD", "HLA-DRA", "CD74", "RGS5"),
  TNK     = c("CD3D", "CD3E", "CD4", "CD8A", "IL7R", "CCR7",
              "FOXP3", "CTLA4", "GZMK", "GNLY", "KLRD1", "MKI67"),
  Bcell   = c("MS4A1", "CD79A", "BANK1", "TCL1A", "CD27", "AIM2",
              "MZB1", "XBP1", "TNFRSF17", "JCHAIN", "BCL6", "AICDA"),
  myeloid = c("CD68", "C1QA", "APOE", "SPP1", "CD14", "FCN1",
              "CSF3R", "FCGR3B", "CLEC9A", "CD1C", "LAMP3", "LILRA4", "TPSAB1", "MKI67"),
  endothelial = c("PECAM1", "VWF", "CLDN5", "GJA5", "HEY1", "ACKR1",
                  "SELP", "CA4", "RGCC", "ESM1", "PROX1", "CCL21", "MKI67")
)
mk <- intersect(marker_panels[[SUBSET]], rownames(sub))

ts("writing ", PDFF)
pdf(PDFF, width = 11, height = 8.5)

print(ElbowPlot(sub, ndims = N_PCS) +
        ggtitle(sprintf("%s - PC stdev (using %d dims)", SUBSET, N_DIMS)) +
        geom_vline(xintercept = N_DIMS, linetype = 2, colour = "red"))

print(DimPlot(sub, reduction = "umap_unintegrated", group.by = "dataset_id",
              shuffle = TRUE, raster = TRUE, pt.size = 0.3) +
        ggtitle(sprintf("%s - BEFORE Harmony (mixing %.3f)", SUBSET, mix_pca)))

print(DimPlot(sub, reduction = "umap", group.by = "dataset_id",
              shuffle = TRUE, raster = TRUE, pt.size = 0.3) +
        ggtitle(sprintf("%s - AFTER Harmony (mixing %.3f, target %.3f)",
                        SUBSET, mix_harm, mix_ideal)))

print(DimPlot(sub, reduction = "umap", group.by = "dataset_id",
              split.by = "dataset_id", ncol = 3, raster = TRUE, pt.size = 0.3) +
        NoLegend() + ggtitle(sprintf("%s - after Harmony, split by dataset", SUBSET)))

print(DimPlot(sub, reduction = "umap_alt", group.by = "dataset_id",
              shuffle = TRUE, raster = TRUE, pt.size = 0.3) +
        ggtitle(sprintf("%s - COMPARISON: Harmony on %s (mixing %.3f)", SUBSET, ALT_KEY, mix_alt)))
print(DimPlot(sub, reduction = "umap_alt", group.by = paste0(RES_PREFIX, 0.6),
              label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.3) + NoLegend() +
        ggtitle(sprintf("%s - Harmony on %s, coloured by the %s res-0.6 clusters (ARI %.3f)",
                        SUBSET, ALT_KEY, HARMONY_KEY, ari_06)))
print(DimPlot(sub, reduction = "umap", group.by = "tissue_class",
              shuffle = TRUE, raster = TRUE, pt.size = 0.3) +
        ggtitle(sprintf("%s - tissue_class (Harmony on %s)", SUBSET, HARMONY_KEY)))

if ("celltype_level3" %in% colnames(sub@meta.data)) {
  print(DimPlot(sub, reduction = "umap", group.by = "celltype_level3",
                label = TRUE, repel = TRUE, raster = TRUE, pt.size = 0.3) +
          ggtitle(sprintf("%s - inherited celltype_level3 (context only, NOT trusted)", SUBSET)))
}

for (r in RESOLUTIONS) {
  col <- paste0(RES_PREFIX, r)
  print(DimPlot(sub, reduction = "umap", group.by = col, label = TRUE,
                repel = TRUE, raster = TRUE, pt.size = 0.3) + NoLegend() +
          ggtitle(sprintf("%s - resolution %s (%d clusters)",
                          SUBSET, r, nlevels(sub[[col]][, 1]))))
}

if (length(mk) > 0) {
  print(FeaturePlot(sub, features = mk, reduction = "umap", ncol = 4,
                    raster = TRUE, pt.size = 0.3, order = TRUE) &
          theme(plot.title = element_text(size = 9)))
  print(DotPlot(sub, features = mk, group.by = paste0(RES_PREFIX, 0.6)) +
          RotatedAxis() + ggtitle(sprintf("%s - markers at res 0.6", SUBSET)))
}

print(VlnPlot(sub, features = intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                          "decontX_contamination"),
                                        colnames(sub@meta.data)),
              group.by = paste0(RES_PREFIX, 0.6), pt.size = 0, ncol = 1) & NoLegend())

dev.off()

# ============================ 9. save + provenance ===========================
ts("writing ", OUT, " (qs2)")
save_obj(sub, OUT)
ts("wrote ", sprintf("%.2f GB", file.size(OUT) / 1e9))

prov <- c(
  sprintf("PROVENANCE - Step 3 pipeline - subset '%s'", SUBSET),
  sprintf("generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("script   : rebuild_03_pipeline.R"),
  sprintf("input    : %s", IN),
  sprintf("output   : %s", OUT),
  "",
  "ANNOTATION PROVENANCE",
  "  subset defined by celltype_level2 plus the reviewed res-0.6 routing",
  "  (see PROVENANCE_02_subset.txt).",
  sprintf("  assay used: %s (cluster columns %s<r>)", ASSAY, RES_PREFIX),
  "  celltype_level3 is carried for context only and is NOT trusted -- the",
  "  hierarchical pass relabelled 98.2% of the T/NK compartment as unresolved.",
  "",
  "BATCHES",
  sprintf("  min cells per batch : %d", MIN_CELLS_PER_BATCH),
  sprintf("  kept    : %s", paste(keep_ds, collapse = ", ")),
  sprintf("  dropped : %s", if (length(drop_ds)) paste(drop_ds, collapse = ", ") else "(none)"),
  sprintf("  cells   : %s", format(ncol(sub), big.mark = ",")),
  "",
  "HVG MASK (feature selection only -- masked genes remain in the object and",
  "are still plottable as markers at Step 4)",
  sprintf("  regex   : %s", mask_rx),
  sprintf("  ambient : %s", paste(ambient, collapse = ", ")),
  sprintf("  stress  : %s + %s (%d genes present: %s)", paste(stress, collapse = ", "), stress_rx,
          length(m_stress), paste(sort(m_stress), collapse = ", ")),
  sprintf("  batch   : %d genes from exclude_batchgenes.rds[['%s']]", length(m_batch), SUBSET),
  sprintf("  masked total %d of %d genes", length(masked), length(g)),
  sprintf("  vst pool %d -> %d masked out -> %d HVGs used",
          length(hvg_pool), n_removed, length(hvg)),
  "",
  "EMBEDDING",
  sprintf("  PCs computed %d, dims used 1:%d, seed %d", N_PCS, N_DIMS, SEED),
  sprintf("  Harmony on %s", HARMONY_KEY),
  sprintf("  kNN batch mixing by dataset (k=%d): before %.3f | after %.3f | target %.3f",
          K_MIX, mix_pca, mix_harm, mix_ideal),
  sprintf("  same-tissue_class neighbours: PCA %.3f | after %.3f | random %.3f",
          tis_pca, tis_harm, tis_rand),
  sprintf("  COMPARISON, Harmony on %s: mixing %.3f | same-tissue %.3f | ARI vs main at res 0.6 %.3f",
          ALT_KEY, mix_alt, tis_alt, ari_06),
  "    (reduction harmony_alt, umap_alt, column alt_res.0.6; not used downstream)",
  sprintf("  resolutions swept: %s", paste(RESOLUTIONS, collapse = ", ")),
  "",
  "UPSTREAM EXCLUSIONS (rebuild_02_subset.R, PROVENANCE_02_subset.txt)",
  "  - manifest: peripheral blood, liver/other metastases, non-PDAC tumours, and",
  "    all of GSE211644 (sorted, CD8-biased) are not in this object.",
  "  - scDblFinder doublets (per sample, whole atlas) removed.",
  "  - 'Mixed/possible doublet' (celltype_level2) excluded by construction.",
  "  - v4: per dataset x compartment QC floor and reviewed 'drop' clusters removed;",
  "    misrouted cells moved between subsets (rebuild_02a_routing.R).",
  "",
  "KNOWN CONSTRAINTS",
  "  - CRA001160 and GSE155698 share an older reference build with the IG_C and",
  "    TR_C biotypes filtered out. Consequences: no gamma-delta calls in those",
  "    two datasets; no cross-dataset isotype work; CD24 unusable atlas-wide.",
  "  - GSE154778 has 45% percent.mt == 0 atlas-wide, but 0-5% inside these",
  "    compartments. Never apply a uniform percent.mt filter -- doing so",
  "    manufactures a batch effect out of a QC step.",
  "  - Ambient RNA was removed by the authors for GSE205013 and GSE212966 only;",
  "    v4 runs DecontX on every dataset (rebuild_01d_decontx.R).",
  ""
)
writeLines(prov, PROV)
ts("wrote ", PROV)

ts("done")
print(sessionInfo())
