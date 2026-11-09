#!/usr/bin/env Rscript
# =============================================================================
# rebuild_04c_subcluster.R   --   split one cluster into subclusters, in a new column
#
# Usage:  Rscript rebuild_04c_subcluster.R <subset> <res> <cluster> <markersA> <markersB> [nameA] [nameB]
#         markers are comma-separated, e.g.  XCR1,CLEC9A  LAMP3,CCR7,CCL22
#
# Writes column  subcluster_res<res>  into subset_<subset>_harmony.qs2 (re-saves
# it, like 04/04b): every cell keeps its cluster id, except cells of <cluster>,
# which become <cluster>a / <cluster>b. Also writes
#   step4c_<subset>_res<res>_c<cluster>_subcluster.tsv   (per cell)
#   PROVENANCE_step4c_<subset>_res<res>_c<cluster>.txt
#
# Why marker-defined, not graph-defined: the use case (v4, 2026-11-09) is
# myeloid cDC1 vs mregDC. They sit in one cluster up to resolution 1.6 but are
# two separate groups by marker (in run 1, 169 cDC1 and 144 mregDC cells, only
# 2 expressing both), so more resolution does not split them. Each cell is
# scored on the two panels (fraction of panel genes detected, decontX
# counts); a cell goes to the panel it scores higher on, and is
# 'u' (unassigned) when it detects neither panel or ties. Seurat's
# FindSubCluster result on the same cells is reported alongside as a check.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table) })
source("/ix/tbruno/anp795/PDAC_Atlas/code/rebuild/io.R")

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 5) stop("usage: <subset> <res> <cluster> <markersA> <markersB> [nameA] [nameB]")
SUBSET <- a[1]; RES <- a[2]; CL <- a[3]
MA <- strsplit(a[4], ",")[[1]]; MB <- strsplit(a[5], ",")[[1]]
NA_ <- if (length(a) >= 6) a[6] else "A"; NB_ <- if (length(a) >= 7) a[7] else "B"
stopifnot(SUBSET %in% SUBSETS)
DIR <- V4_DIR
IN  <- file.path(DIR, sprintf("subset_%s_harmony.qs2", SUBSET))
OUT <- file.path(DIR, sprintf("step4c_%s_res%s_c%s_subcluster.tsv", SUBSET, RES, CL))
PROV <- file.path(DIR, sprintf("PROVENANCE_step4c_%s_res%s_c%s.txt", SUBSET, RES, CL))
ts <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))

sub <- read_obj(IN)
rc <- res_col_of(sub, RES)
cl <- as.character(sub[[rc]][, 1])
ii <- which(cl == CL)
if (!length(ii)) stop("cluster ", CL, " not found in ", rc)
ts(SUBSET, " ", rc, " cluster ", CL, ": ", length(ii), " cells")

cnt <- LayerData(sub, assay = DefaultAssay(sub), layer = "counts")
MA <- intersect(MA, rownames(cnt)); MB <- intersect(MB, rownames(cnt))
fa <- Matrix::colMeans(cnt[MA, ii, drop = FALSE] > 0)
fb <- Matrix::colMeans(cnt[MB, ii, drop = FALSE] > 0)
call <- ifelse(fa > fb, "a", ifelse(fb > fa, "b", "u"))
call[fa == 0 & fb == 0] <- "u"

# graph-based check
Idents(sub) <- rc
gname <- paste0(DefaultAssay(sub), "_snn")
fs <- tryCatch({
  s2 <- FindSubCluster(sub, cluster = CL, graph.name = gname, subcluster.name = "fsc", resolution = 0.5)
  as.character(s2$fsc[ii])
}, error = function(e) rep(NA_character_, length(ii)))

lab <- cl; lab[ii] <- paste0(CL, call)
newcol <- paste0("subcluster_res", RES)
sub[[newcol]] <- lab
save_obj(sub, IN)
ts("re-wrote ", IN, " with column ", newcol)

det <- function(g, idx) round(100 * mean(cnt[g, ii[idx]] > 0), 1)
tab <- rbindlist(lapply(c("a", "b", "u"), function(k) {
  idx <- which(call == k)
  if (!length(idx)) return(NULL)
  data.table(sub = paste0(CL, k), name = c(a = NA_, b = NB_, u = "unassigned")[[k]], n = length(idx),
             t(setNames(vapply(c(MA, MB), function(g) det(g, idx), numeric(1)), c(MA, MB))))
}), fill = TRUE)
both <- sum(fa > 0 & fb > 0)
print(tab); print(table(marker_call = call, FindSubCluster = fs, useNA = "ifany"))

fwrite(data.table(cell = colnames(sub)[ii], cluster = CL, subcluster = paste0(CL, call),
                  frac_A = round(fa, 3), frac_B = round(fb, 3), findsubcluster = fs), OUT, sep = "\t")
show <- function(x) capture.output(print(x, row.names = FALSE))
writeLines(c(
  sprintf("PROVENANCE - Step 4c subcluster - %s %s cluster %s", SUBSET, rc, CL),
  sprintf("generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "script   : rebuild_04c_subcluster.R",
  sprintf("panel a (%s): %s", NA_, paste(MA, collapse = ", ")),
  sprintf("panel b (%s): %s", NB_, paste(MB, collapse = ", ")),
  "rule     : a cell goes to the panel with the higher fraction of genes detected (decontX counts); u = neither or tie",
  sprintf("cells detecting both panels: %d of %d", both, length(ii)),
  sprintf("new column: %s (other clusters unchanged)", newcol), "",
  "% of cells detecting each gene, by subcluster:", show(tab), "",
  "agreement with Seurat FindSubCluster (resolution 0.5):", capture.output(print(table(marker_call = call, FindSubCluster = fs, useNA = "ifany")))
), PROV)
ts("wrote ", OUT, " and ", PROV)
