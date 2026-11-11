#!/usr/bin/env Rscript
# =============================================================================
# rebuild_05b_figures.R   --   figures of one final subset (step 05 output)
#
# Usage:  Rscript rebuild_05b_figures.R <subset>
#
# Input : <V4_DIR>/subset_<subset>_final.qs2          (read only)
# Output: <V4_DIR>/figures_step5/<subset>/
#           1_overview.png     UMAP coloured by lineage, cluster ids on the plot, key at the right
#           2_celltypes.png    one panel per celltype_v4: n cells, pool, % kept out of DE/LR
#           3_batch.png        one panel per dataset and per tissue class
#           4_flags.png        cells kept out of DE/LR, one panel per cell_flag, dissociation score
#           5_markers.png      dot plot: lineage core genes + every gene named in the labels
#           6_composition.png  share of each sample's cells per celltype_pool_v4, by tissue class
#           PROVENANCE_step5b_<subset>.txt
#
# The UMAP is the run-2 embedding ("umap", Harmony on dataset_id) restricted to
# the kept cells. Step 05 drops cells without re-embedding, so the gaps where
# dropped clusters sat are expected.
#
# Colour. Identity never rests on colour alone. At most three lineage hues (the
# palette slots that stay apart on a scatter for every pair, incl. colour-vision
# deficiency); a fourth lineage and beyond fold into the third hue, and the
# cluster ids printed on the plot carry identity. Anything with more groups is
# drawn as small multiples: one highlighted group per panel on all cells.
# Magnitudes use one blue ramp, light to dark.
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table)
  library(ggplot2); library(patchwork); library(ggrepel)
})
source("/ix/tbruno/anp795/PDAC_Atlas/code/rebuild/io.R")

a <- commandArgs(trailingOnly = TRUE)
if (!length(a)) stop("usage: rebuild_05b_figures.R <subset>")
SUBSET <- a[1]
if (!SUBSET %in% SUBSETS) stop("unknown subset ", SUBSET, "; one of ", paste(SUBSETS, collapse = ", "))

IN   <- file.path(V4_DIR, sprintf("subset_%s_final.qs2", SUBSET))
OUTD <- file.path(V4_DIR, "figures_step5", SUBSET)
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
PROV <- file.path(OUTD, sprintf("PROVENANCE_step5b_%s.txt", SUBSET))
MIN_SAMPLE_CELLS <- 20     # composition: samples with fewer kept cells are not shown
MAX_GENES        <- 80     # dot plot
DPI              <- 300
COLW             <- 2.35   # inches per column in the small-multiple grids
TW               <- 28     # characters per line: panel titles (7 pt bold) ...
SW               <- 40     # ... and panel subtitles (6 pt), so they stay inside a COLW column
NAME <- c(TNK = "T/NK", stromal = "Stromal", Bcell = "B cell", myeloid = "Myeloid",
          endothelial = "Endothelial")[[SUBSET]]

# ---- palette (light mode of the reference data-viz palette) ------------------
SURF <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"
GRID <- "#e1e0d9"; ZERO <- "#f0efec"
HI   <- "#2a78d6"                              # highlighted cells / single series
BGPT <- "#dcdbd4"                              # all other cells
CAT3 <- c("#2a78d6", "#eb6834", "#1baf7a")     # slots 1-3, validated all-pairs for scatter
SEQ  <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")

# Lineage core genes for the dot plot; the genes named in the labels are added.
CORE <- list(
  TNK         = c("CD3E", "CD4", "CD8A", "CD8B", "FOXP3", "CCR7", "LEF1", "IL7R", "GZMK", "GZMB",
                  "NKG7", "GNLY", "KLRF1", "FCGR3A", "NCAM1", "KIT", "RORC", "MKI67"),
  stromal     = c("PDGFRA", "DCN", "LUM", "COL1A1", "PDGFRB", "RGS5", "ACTA2", "MYH11", "S100B",
                  "PLP1", "MKI67"),
  Bcell       = c("MS4A1", "CD79A", "CD19", "TCL1A", "IGHD", "CD27", "MZB1", "JCHAIN", "XBP1",
                  "SDC1", "IGHG1", "IGHA1", "MKI67"),
  myeloid     = c("CD68", "CD14", "FCGR3A", "LYZ", "C1QA", "APOE", "CSF3R", "FCGR3B", "CXCR2",
                  "TPSAB1", "KIT", "CLEC9A", "CD1C", "LAMP3", "LILRA4", "MKI67"),
  endothelial = c("PECAM1", "CDH5", "VWF", "PLVAP", "CA4", "ACKR1", "SELE", "GJA5", "SEMA3G",
                  "PROX1", "MKI67"))

ts  <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))
fmt <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
wrap <- function(x, w) vapply(x, function(s) paste(strwrap(s, w), collapse = "\n"), "")
pad_lines <- function(x) {                     # same line count, so panel titles align
  k <- lengths(strsplit(x, "\n", fixed = TRUE)); m <- max(k)
  paste0(x, strrep("\n ", m - k))
}
nlines <- function(x) max(lengths(strsplit(x, "\n", fixed = TRUE)))
save_png <- function(p, f, w, h) {
  ggsave(file.path(OUTD, f), p, width = w, height = h, dpi = DPI, device = ragg::agg_png, bg = SURF)
  ts("wrote ", f, sprintf(" (%.1f x %.1f in)", w, h))
  f
}

# ---- data -----------------------------------------------------------------------
ts("step 05b ", SUBSET, ": reading ", IN)
o  <- read_obj(IN)
md <- o@meta.data
need <- c("cluster_v4", "subcluster_v4", "lineage_v4", "celltype_v4", "celltype_pool_v4",
          "cell_flag", "exclude_from_de_lr", "sample_key", "dataset_id", "tissue_class")
if (length(miss <- setdiff(need, names(md)))) stop("final object lacks ", paste(miss, collapse = ", "), "; rerun step 05")
if (!"umap" %in% Reductions(o)) stop("final object has no 'umap' reduction")
um <- Embeddings(o, "umap")
d <- data.table(
  cell = rownames(md), x = um[, 1], y = um[, 2],
  id = as.character(md$subcluster_v4), lineage = as.character(md$lineage_v4),
  celltype = as.character(md$celltype_v4), pool = as.character(md$celltype_pool_v4),
  flag = ifelse(is.na(md$cell_flag), "", as.character(md$cell_flag)),
  excl = as.logical(md$exclude_from_de_lr) %in% TRUE,
  sample = as.character(md$sample_key), dataset = as.character(md$dataset_id),
  tissue = as.character(md$tissue_class),
  dissoc = if ("score_dissociation" %in% names(md)) md$score_dissociation else NA_real_)
N <- nrow(d)
if (nrow(bad <- d[, .(k = uniqueN(celltype)), by = id][k > 1])) stop("ids with more than one label: ", paste(bad$id, collapse = ", "))
ts(fmt(N), " cells, ", uniqueN(d$id), " cell types, ", uniqueN(d$sample), " samples")

id_rank <- function(v) order(as.integer(sub("[a-z]+$", "", v)), sub("^[0-9]+", "", v))
grp <- d[, .(n = .N, lineage = lineage[1], celltype = celltype[1], pool = pool[1],
             pct_excl = 100 * mean(excl), x = median(x), y = median(y)), by = id]
grp <- grp[id_rank(id)]

# Lineage hues: up to three; beyond that the third hue covers the smaller lineages.
lin_n <- d[, .N, by = lineage][order(-N)]
lin_map <- if (nrow(lin_n) <= 3) setNames(lin_n$lineage, lin_n$lineage) else
  setNames(c(lin_n$lineage[1:2], rep(paste(lin_n$lineage[-(1:2)], collapse = " / "), nrow(lin_n) - 2)), lin_n$lineage)
lg_cols <- setNames(CAT3[seq_along(unique(lin_map))], unique(lin_map))   # largest lineage first
d[, lgroup := factor(lin_map[lineage], levels = names(lg_cols))]

XL <- range(d$x); YL <- range(d$y)
ASP <- diff(YL) / diff(XL)
row_h <- function(nt, ns) (COLW - 0.15) * ASP + 0.105 * nt + 0.09 * ns + 0.2   # one grid row, inches
psize <- function(scale = 1) max(0.04, min(0.6, 0.5 * sqrt(1e4 / N))) * scale
th_umap <- theme_void(base_size = 8) + theme(
  plot.background  = element_rect(fill = SURF, colour = NA),
  panel.background = element_rect(fill = SURF, colour = NA),
  plot.title    = element_text(colour = INK, size = 7, face = "bold", lineheight = 1.05, margin = margin(b = 1.5)),
  plot.subtitle = element_text(colour = INK2, size = 6, lineheight = 1.05, margin = margin(b = 2)),
  plot.margin = margin(4, 5, 4, 5), legend.position = "none")
th_page <- theme(plot.background = element_rect(fill = SURF, colour = NA),
                 plot.title = element_text(colour = INK, size = 11, face = "bold"),
                 plot.caption = element_text(colour = MUTED, size = 6.5, hjust = 0))
umap_coord <- coord_equal(xlim = XL, ylim = YL, expand = TRUE)

hl_panel <- function(sel, title, subtitle, ps) {
  ggplot() +
    geom_point(data = d[!sel], aes(x, y), colour = BGPT, size = ps, shape = 16, stroke = 0) +
    geom_point(data = d[sel],  aes(x, y), colour = HI,   size = ps, shape = 16, stroke = 0) +
    umap_coord + labs(title = title, subtitle = subtitle) + th_umap
}
grid_dims <- function(k) { nc <- if (k > 12) 5 else if (k > 6) 4 else max(k, 1); c(nc, ceiling(k / nc)) }
blank <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = SURF, colour = NA))
# specs: list of list(title, sub, sel = logical) or list(title, sub, plot = ggplot without labels).
# Titles and subtitles are wrapped and padded together, so every panel in a grid lines up.
make_panels <- function(specs) {
  tt <- pad_lines(wrap(vapply(specs, `[[`, "", "title"), TW))
  ss <- pad_lines(wrap(vapply(specs, `[[`, "", "sub"), SW))
  pl <- lapply(seq_along(specs), function(i) {
    s <- specs[[i]]
    if (is.null(s$sel)) s$plot + labs(title = tt[i], subtitle = ss[i]) else hl_panel(s$sel, tt[i], ss[i], psize(0.6))
  })
  list(plots = pl, h = row_h(nlines(tt), nlines(ss)))
}
files <- character(0)

# ---- 1. overview ------------------------------------------------------------------
set.seed(1)
p_ov <- ggplot(d[base::sample(.N)], aes(x, y, colour = lgroup)) +   # random draw order: no lineage hides another
  geom_point(size = psize(1.15), shape = 16, stroke = 0) +
  geom_text_repel(data = grp, aes(x, y, label = id), inherit.aes = FALSE, colour = INK,
                  bg.color = SURF, bg.r = 0.18, size = 2.7, fontface = "bold", seed = 1,
                  box.padding = 0.15, point.size = NA, min.segment.length = Inf, max.overlaps = Inf) +
  scale_colour_manual(values = lg_cols, breaks = names(lg_cols), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.6))) +
  umap_coord + th_umap +
  theme(legend.position = if (length(lg_cols) > 1) "inside" else "none",
        legend.position.inside = c(0.01, 0.99), legend.justification = c(0, 1),
        legend.text = element_text(colour = INK2, size = 7), legend.key.size = unit(9, "pt"),
        legend.background = element_rect(fill = SURF, colour = NA))
key <- copy(grp)[, row := .I][, txt := sprintf("%s   %s   (%s)", id, celltype, fmt(n))]
p_key <- ggplot(key, aes(y = -row)) +
  geom_point(aes(x = 0, colour = lin_map[lineage]), size = 1.8, shape = 16) +
  geom_text(aes(x = 0.035, label = txt), hjust = 0, size = 2.25, colour = INK) +
  scale_colour_manual(values = lg_cols, guide = "none") +
  scale_x_continuous(limits = c(-0.02, 1), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(add = 0.9)) +
  theme_void() + theme(plot.background = element_rect(fill = SURF, colour = NA))
p1 <- (p_ov | p_key) + plot_layout(widths = c(1.45, 1)) +
  plot_annotation(
    title = sprintf("%s: %s kept cells in %d cell types", NAME, fmt(N), nrow(grp)),
    theme = th_page)
files <- c(files, save_png(p1, "1_overview.png", 12, 6.4))

# ---- 2. one panel per cell type ---------------------------------------------------------
m2 <- make_panels(lapply(seq_len(nrow(grp)), function(i) list(
  sel = d$id == grp$id[i],
  title = sprintf("%s  %s", grp$id[i], grp$celltype[i]),
  sub = paste0(fmt(grp$n[i]), " cells",
               if (grp$pool[i] != grp$celltype[i]) paste0(" · pool: ", grp$pool[i]) else "",
               if (grp$pct_excl[i] > 0) sprintf(" · %.0f%% out of DE/LR", grp$pct_excl[i]) else ""))))
gd <- grid_dims(length(m2$plots))
p2 <- wrap_plots(m2$plots, ncol = gd[1]) +
  plot_annotation(title = sprintf("%s: each cell type on the UMAP", NAME), theme = th_page)
files <- c(files, save_png(p2, "2_celltypes.png", COLW * gd[1], gd[2] * m2$h + 0.55))

# ---- 3. batch: dataset and tissue class -----------------------------------------------
spec_by <- function(col, lab) {
  tb <- d[, .(n = .N, s = uniqueN(sample)), by = c(col)][order(-n)]
  lapply(seq_len(nrow(tb)), function(i) list(
    sel = d[[col]] == tb[[col]][i], title = paste0(lab, ": ", gsub("_", " ", tb[[col]][i])),
    sub = sprintf("%s cells (%.0f%%) · %d samples", fmt(tb$n[i]), 100 * tb$n[i] / N, tb$s[i])))
}
s_ds <- spec_by("dataset", "dataset"); s_tc <- spec_by("tissue", "tissue")
m3 <- make_panels(c(s_ds, s_tc))
p_ds <- m3$plots[seq_along(s_ds)]; p_tc <- m3$plots[-seq_along(s_ds)]
nc <- max(length(p_ds), length(p_tc))
fill_row <- function(l) c(l, rep(list(blank), nc - length(l)))
p3 <- wrap_plots(c(fill_row(p_ds), fill_row(p_tc)), ncol = nc) +
  plot_annotation(title = sprintf("%s: cells by dataset (top) and tissue class (bottom)", NAME), theme = th_page)
files <- c(files, save_png(p3, "3_batch.png", COLW * nc, 2 * m3$h + 0.55))

# ---- 4. flags ---------------------------------------------------------------------------
fl <- d[nzchar(flag), .(n = .N, ex = 100 * mean(excl)), by = flag][order(-n)]
s_fl <- c(list(list(sel = d$excl, title = "Kept out of DE / ligand-receptor",
                    sub = sprintf("%s cells (%.1f%%)", fmt(sum(d$excl)), 100 * mean(d$excl)))),
          lapply(seq_len(nrow(fl)), function(i) list(
            sel = d$flag == fl$flag[i], title = paste0("flag: ", fl$flag[i]),
            sub = sprintf("%s cells (%.1f%%); %s out of DE/LR", fmt(fl$n[i]), 100 * fl$n[i] / N,
                          if (fl$ex[i] >= 99.95) "all" else sprintf("%.0f%%", fl$ex[i])))))
if (any(!is.na(d$dissoc))) {
  q <- quantile(d$dissoc, c(0.01, 0.99), na.rm = TRUE)
  dd <- d[!is.na(dissoc)][, v := pmin(pmax(dissoc, q[1]), q[2])][order(v)]
  s_fl <- c(s_fl, list(list(
    title = "Dissociation score", sub = "van den Brink core set; clipped at the 1st / 99th percentile",
    plot = ggplot(dd, aes(x, y, colour = v)) + geom_point(size = psize(0.6), shape = 16, stroke = 0) +
      scale_colour_gradientn(colours = c(ZERO, SEQ), name = NULL, breaks = scales::breaks_pretty(3)) +
      umap_coord + th_umap +
      theme(legend.position = "right", legend.key.width = unit(5, "pt"), legend.key.height = unit(14, "pt"),
            legend.text = element_text(colour = INK2, size = 5.5)))))
}
m4 <- make_panels(s_fl)
gd <- grid_dims(length(m4$plots))
p4 <- wrap_plots(m4$plots, ncol = gd[1]) +
  plot_annotation(title = sprintf("%s: flagged cells", NAME), theme = th_page)
files <- c(files, save_png(p4, "4_flags.png", COLW * gd[1], gd[2] * m4$h + 0.55))

# ---- 5. marker dot plot -------------------------------------------------------------
lab_genes <- unique(unlist(regmatches(grp$celltype, gregexpr("[A-Za-z0-9.-]+(?=\\+)", grp$celltype, perl = TRUE))))
genes_all <- unique(c(CORE[[SUBSET]], lab_genes))
core_missing <- setdiff(CORE[[SUBSET]], rownames(o))
lab_not_gene <- setdiff(setdiff(lab_genes, CORE[[SUBSET]]), rownames(o))   # e.g. "ISG", "IgG": label words, not genes
genes_missing <- c(core_missing, lab_not_gene)
genes <- head(intersect(genes_all, rownames(o)), MAX_GENES)
X  <- LayerData(o, assay = "decontX", layer = "data")[genes, d$cell, drop = FALSE]
gf <- factor(d$id, levels = grp$id)
Fm <- t(Matrix::fac2sparse(gf))                         # cells x groups
nn <- as.numeric(table(gf))
avg <- sweep(as.matrix(X %*% Fm), 2, nn, "/")
pct <- sweep(as.matrix((X > 0) %*% Fm), 2, nn, "/") * 100
rel <- avg / pmax(apply(avg, 1, max), 1e-9)
top <- apply(rel, 1, which.max)
genes <- genes[order(top, -pct[cbind(seq_along(top), top)])]       # roughly diagonal
ylab <- setNames(sprintf("%s  %s", grp$id, ifelse(nchar(grp$celltype) > 58, paste0(substr(grp$celltype, 1, 57), "…"), grp$celltype)), grp$id)
dp <- data.table(gene = rep(rownames(rel), ncol(rel)), id = rep(colnames(rel), each = nrow(rel)),
                 rel = as.vector(rel), pct = as.vector(pct))
dp[, gene := factor(gene, levels = genes)]
dp[, lab := factor(ylab[id], levels = rev(ylab))]
W5 <- max(7, 3.6 + 0.135 * length(genes))
wrapW <- function(s) paste(strwrap(s, floor((W5 - 0.4) / 0.062)), collapse = "\n")   # ~0.062 in per 8-pt character
cap5 <- paste(c(if (length(core_missing)) paste("Core genes not in the object:", paste(core_missing, collapse = ", ")),
                if (length(lab_not_gene)) paste("Label terms that are not gene names (not plotted):", paste(lab_not_gene, collapse = ", "))),
              collapse = ". ")
p5 <- ggplot(dp[pct > 0], aes(gene, lab)) +
  geom_point(aes(size = pct, colour = rel), shape = 16) +
  scale_size_area(max_size = 3.4, limits = c(0, 100), breaks = c(25, 50, 75, 100), name = "% of cells\nexpressing") +
  scale_colour_gradientn(colours = SEQ, limits = c(0, 1), breaks = c(0, 0.5, 1),
                         name = "mean expression,\nshare of the\ngene's maximum") +
  labs(x = NULL, y = NULL, title = sprintf("%s: marker genes by cell type", NAME),
       caption = if (nzchar(cap5)) wrapW(cap5) else NULL) +
  theme_minimal(base_size = 7) +
  theme(plot.background = element_rect(fill = SURF, colour = NA),
        plot.title.position = "plot", plot.caption.position = "plot",
        panel.grid.major = element_line(colour = GRID, linewidth = 0.2), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic", colour = INK2, size = 6),
        axis.text.y = element_text(colour = INK, size = 6.5),
        legend.title = element_text(colour = INK2, size = 6.5), legend.text = element_text(colour = INK2, size = 6),
        plot.title = element_text(colour = INK, size = 11, face = "bold"),
        plot.caption = element_text(colour = MUTED, size = 6.5, hjust = 0))
files <- c(files, save_png(p5, "5_markers.png", W5, 1.4 + 0.2 * nrow(grp)))

# ---- 6. composition -------------------------------------------------------------------
TIS <- c("primary_tumor", "adjacent_normal", "normal_pancreas")
sn <- d[, .(n = .N), by = .(sample, dataset, tissue)]
shown <- sn[n >= MIN_SAMPLE_CELLS]
shown[, tissue := factor(tissue, levels = c(intersect(TIS, unique(tissue)), setdiff(unique(tissue), TIS)))]
setorder(shown, tissue, dataset, sample)
pool_lv <- unique(grp$pool)
cm <- d[sample %in% shown$sample, .N, by = .(sample, pool)][, share := N / sum(N), by = sample]
cg <- CJ(sample = shown$sample, pool = pool_lv)[cm, on = .(sample, pool), share := i.share]
cg <- shown[, .(sample, tissue)][cg, on = "sample"]
cg[, `:=`(sample = factor(sample, levels = shown$sample), pool = factor(pool, levels = rev(pool_lv)))]
shown[, sample := factor(sample, levels = shown$sample)]
# Tissue class is a coloured strip with a legend, not facet labels: a facet holding two
# samples is too narrow for its name. The facets still separate the groups.
# Palette slots 2-4 (orange, aqua, yellow; validated as adjacent pairs), not blue: blue is the heatmap's magnitude ramp.
TISC    <- c("#eb6834", "#1baf7a", "#eda100")
tis_lv  <- levels(shown$tissue)
tis_col <- setNames(TISC[seq_along(tis_lv)], tis_lv)
tis_n   <- shown[, .N, by = tissue]
tis_lab <- setNames(sprintf("%s (%d)", gsub("_", " ", tis_lv), tis_n$N[match(tis_lv, tis_n$tissue)]), tis_lv)
th_comp <- theme_minimal(base_size = 7) +
  theme(plot.background = element_rect(fill = SURF, colour = NA), strip.text = element_blank(),
        legend.title = element_text(colour = INK2, size = 6.5), legend.text = element_text(colour = INK2, size = 6))
p_bar <- ggplot(shown, aes(sample, n)) + geom_col(fill = HI, width = 0.72) +
  scale_y_log10(labels = scales::label_comma(), breaks = c(20, 100, 1000, 10000)) +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "kept cells\n(log scale)") + th_comp +
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = GRID, linewidth = 0.2),
        axis.text.y = element_text(colour = INK2, size = 6), axis.title.y = element_text(colour = INK2, size = 6.5))
p_tis <- ggplot(shown, aes(sample, 1, fill = tissue)) + geom_tile(colour = SURF, linewidth = 0.25) +
  scale_fill_manual(values = tis_col, breaks = tis_lv, labels = tis_lab, name = "tissue class (samples)") +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL) + th_comp +
  theme(axis.text = element_blank(), panel.grid = element_blank())
p_heat <- ggplot(cg, aes(sample, pool, fill = share)) +
  geom_tile(colour = SURF, linewidth = 0.25) +
  scale_fill_gradientn(colours = SEQ, na.value = ZERO, transform = "sqrt", limits = c(0, max(cg$share, na.rm = TRUE)),
                       labels = scales::label_percent(accuracy = 1), name = "share of the\nsample's cells\n(sqrt scale)") +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL) + th_comp +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, colour = INK2, size = 4.6),
        axis.text.y = element_text(colour = INK, size = 6.5))
n_hidden <- nrow(sn) - nrow(shown)
cap <- sprintf("%d samples with at least %d kept cells; %d with fewer are not shown. Light grey tile: no cells of that type.",
               nrow(shown), MIN_SAMPLE_CELLS, n_hidden)
if (SUBSET == "myeloid") cap <- paste(cap, "GSE205013:P26 is left out of myeloid composition comparisons (decision 2026-11-09).")
p6 <- (p_bar / p_tis / p_heat) +
  plot_layout(heights = c(1, 0.1, max(3, 0.32 * length(pool_lv))), guides = "collect") +
  plot_annotation(title = sprintf("%s: composition per sample", NAME), caption = cap, theme = th_page)
files <- c(files, save_png(p6, "6_composition.png", max(8, 3.8 + 0.105 * nrow(shown)), 2.6 + 0.2 * length(pool_lv) + 0.9))

# ---- provenance ------------------------------------------------------------------------
writeLines(c(
  sprintf("PROVENANCE - Step 5b figures - %s", SUBSET),
  sprintf("generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "script   : rebuild_05b_figures.R",
  sprintf("object   : %s (read only)", IN),
  sprintf("cells    : %s kept, %d cell types (celltype_v4), %d pooled labels, %d samples",
          fmt(N), nrow(grp), length(pool_lv), uniqueN(d$sample)),
  "embedding: umap (run-2 Harmony on dataset_id), kept cells only, not recomputed",
  sprintf("dot plot : %d genes; not in the object: %s", length(genes),
          if (length(genes_missing)) paste(genes_missing, collapse = ", ") else "none"),
  sprintf("composition: %d samples shown (>= %d cells), %d not shown", nrow(shown), MIN_SAMPLE_CELLS, n_hidden),
  "", "files:", paste0("  ", files)
), PROV)
ts("wrote ", PROV)
