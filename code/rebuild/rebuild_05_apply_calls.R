#!/usr/bin/env Rscript
# =============================================================================
# rebuild_05_apply_calls.R   --   apply the reviewed cluster calls to one subset
#
# Usage:  Rscript rebuild_05_apply_calls.R <subset> [res]          (res default 0.6)
#         PDAC_DRAFT_CALLS=1 Rscript rebuild_05_apply_calls.R <subset>   (draft run)
#
# Input : <DIR>/subset_<subset>_harmony.qs2                    (read only, never re-saved)
#         code/calls/cluster_calls_<subset>_v4_res<res>.tsv    (one row per cluster)
#         code/calls/v4_run2_decisions.tsv                     (whole-cluster drops / flags)
#         code/calls/subcluster_labels_v4_res<res>.tsv         (names for subcluster_res<res>)
# Output: <DIR>/subset_<subset>_final.qs2                      (kept cells only)
#         <DIR>/step5_<subset>_res<res>_cells.tsv              (every cell, kept and dropped)
#         <DIR>/step5_<subset>_res<res>_summary.tsv            (per cluster / subcluster)
#         <DIR>/PROVENANCE_step5_<subset>_res<res>.txt
# (DIR = rebuild_v4; a draft run writes all four to rebuild_v4/step5_draft/ instead)
#
# What it does
#   - Stops unless the call table lists exactly the object's clusters with the
#     same n_cells. A mismatch means 03 was rerun and renumbered the clusters,
#     so the calls no longer describe these cells.
#   - The call for a cluster is reviewer_call, or the draft `call` where the
#     reviewer left it empty. A normal run stops if any cluster is unreviewed;
#     PDAC_DRAFT_CALLS=1 accepts draft calls and writes to step5_draft/, so a
#     draft result can never be mistaken for the reviewed one.
#   - keep / drop only. A 'route' call stops the run: moving cells to another
#     subset means re-embedding them there, i.e. a rerun from 02a.
#   - v4_run2_decisions.tsv must agree with the call table (same n_cells and
#     action). Its cell_flag goes to every cell of the cluster that has no
#     flag yet, and its exclude_from_de_lr to every cell of the cluster.
#     Per-cell flags from the v4 run-1 decisions (set in 02) are kept as they are.
#   - Labels: celltype_v4 / lineage_v4 from the call table, replaced by the
#     subcluster label where the cluster was split (rebuild_04c).
#   - The harmony object is only read. Its embeddings and graphs are kept for
#     the kept cells; nothing is re-clustered.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(data.table) })
source("/ix/tbruno/anp795/PDAC_Atlas/code/rebuild/io.R")

a <- commandArgs(trailingOnly = TRUE)
if (!length(a)) stop("usage: rebuild_05_apply_calls.R <subset> [res]")
SUBSET <- a[1]
RES    <- if (length(a) >= 2) a[2] else "0.6"
if (!SUBSET %in% SUBSETS) stop("unknown subset ", SUBSET, "; one of ", paste(SUBSETS, collapse = ", "))
DRAFT  <- Sys.getenv("PDAC_DRAFT_CALLS", "0") == "1"

DIR   <- V4_DIR
OUTD  <- if (DRAFT) file.path(DIR, "step5_draft") else DIR
dir.create(OUTD, showWarnings = FALSE)
IN        <- file.path(DIR, sprintf("subset_%s_harmony.qs2", SUBSET))
CALLS     <- file.path(CALLS_DIR, sprintf("cluster_calls_%s_v4_res%s.tsv", SUBSET, RES))
DEC_F     <- file.path(CALLS_DIR, "v4_run2_decisions.tsv")
SUB_F     <- file.path(CALLS_DIR, sprintf("subcluster_labels_v4_res%s.tsv", RES))
OBJ_OUT   <- file.path(OUTD, sprintf("subset_%s_final.qs2", SUBSET))
CELLS_OUT <- file.path(OUTD, sprintf("step5_%s_res%s_cells.tsv", SUBSET, RES))
SUM_OUT   <- file.path(OUTD, sprintf("step5_%s_res%s_summary.tsv", SUBSET, RES))
PROV      <- file.path(OUTD, sprintf("PROVENANCE_step5_%s_res%s.txt", SUBSET, RES))

msg  <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))
show <- function(x) capture.output(print(x, row.names = FALSE))
die  <- function(x, ...) { if (!is.null(x)) print(x, row.names = FALSE); stop(..., call. = FALSE) }
read_tsv <- function(f) fread(f, sep = "\t", colClasses = "character", na.strings = NULL)

msg("step 05 ", SUBSET, " res ", RES, if (DRAFT) "  (DRAFT: unreviewed calls allowed)" else "")

# ---- call table --------------------------------------------------------------
if (!file.exists(CALLS)) stop("no call table: ", CALLS)
tab  <- read_tsv(CALLS)
need <- c("cluster", "n_cells", "call", "lineage", "celltype_label", "reviewer_call")
if (length(miss <- setdiff(need, names(tab)))) stop(CALLS, ": missing columns ", paste(miss, collapse = ", "))
if (anyDuplicated(tab$cluster)) stop(CALLS, ": a cluster is listed twice")
tab[, n_cells := as.integer(n_cells)]
tab[, reviewed := nzchar(reviewer_call)]
tab[, final_call := ifelse(reviewed, reviewer_call, call)]

if (!DRAFT && any(!tab$reviewed))
  stop(sprintf("%d of %d clusters in %s have no reviewer_call (clusters %s). Review them, or run with PDAC_DRAFT_CALLS=1 for a draft (written to %s).",
               sum(!tab$reviewed), nrow(tab), basename(CALLS), paste(tab$cluster[!tab$reviewed], collapse = ", "),
               file.path(DIR, "step5_draft")), call. = FALSE)
bad <- tab[!final_call %in% c("keep", "drop", "route")]
if (nrow(bad)) die(bad[, .(cluster, call, reviewer_call)], CALLS, ": calls must be keep, drop or route")
bad <- tab[final_call == "route"]
if (nrow(bad)) die(bad[, .(cluster, final_call, route_to = if ("route_to" %in% names(bad)) route_to else NA)],
                   "route calls need a rerun from 02a (the cells must be re-embedded in their new subset); step 05 only keeps or drops")

# ---- object, and the consistency check --------------------------------------
sub <- read_obj(IN)
rc  <- res_col_of(sub, RES)
md  <- sub@meta.data
cl  <- as.character(md[[rc]])
msg("read ", IN, ": ", ncol(sub), " cells, ", length(unique(cl)), " clusters in ", rc)

obs <- data.table(cluster = names(table(cl)), n_obj = as.integer(table(cl)))
chk <- merge(obs, tab[, .(cluster, n_cells)], by = "cluster", all = TRUE)
bad <- chk[is.na(n_obj) | is.na(n_cells) | n_obj != n_cells]
if (nrow(bad)) die(bad, basename(CALLS), " does not match ", rc, " of ", basename(IN),
                   ". Was 03 rerun? Then the clusters are renumbered and the table must be redrafted.")

# ---- whole-cluster decisions (v4 run 2) -------------------------------------
dec <- if (file.exists(DEC_F)) read_tsv(DEC_F) else data.table()
if (nrow(dec)) dec <- dec[subset == SUBSET]
if (nrow(dec)) {
  if (RES != "0.6") stop(basename(DEC_F), " is keyed to the res-0.6 clusters; it cannot be applied at res ", RES)
  if (anyDuplicated(dec$run2_cluster)) stop(DEC_F, ": a cluster is listed twice for ", SUBSET)
  dec[, n_cells := as.integer(n_cells)]
  dec[, exclude_from_de_lr := as.logical(exclude_from_de_lr)]
  if (anyNA(dec$exclude_from_de_lr)) stop(DEC_F, ": exclude_from_de_lr must be TRUE or FALSE")
  m   <- merge(dec, tab[, .(cluster, n_table = n_cells, final_call)], by.x = "run2_cluster", by.y = "cluster", all.x = TRUE)
  bad <- m[is.na(n_table) | n_table != n_cells | action != final_call]
  if (nrow(bad)) die(bad[, .(run2_cluster, n_cells, n_table, action, final_call)],
                     basename(DEC_F), " disagrees with ", basename(CALLS), " (cluster, n_cells or action)")
}

# ---- per-cell calls and labels ----------------------------------------------
ti       <- match(cl, tab$cluster)
action   <- tab$final_call[ti]
src      <- ifelse(tab$reviewed[ti], "reviewed", "draft")
lineage  <- tab$lineage[ti]
celltype <- tab$celltype_label[ti]

scol  <- paste0("subcluster_res", RES)
subcl <- if (scol %in% names(md)) as.character(md[[scol]]) else cl
split <- subcl != cl
sl    <- if (file.exists(SUB_F)) read_tsv(SUB_F) else data.table()
if (nrow(sl)) sl <- sl[subset == SUBSET]
if (any(split) || nrow(sl)) {
  if (!nrow(sl)) stop(scol, " splits clusters ", paste(unique(cl[split]), collapse = ", "),
                      " but ", basename(SUB_F), " has no labels for ", SUBSET)
  sl[, n_cells := as.integer(n_cells)]
  obs2 <- data.table(subcluster = subcl[split])[, .(n_obj = .N), by = subcluster]
  chk2 <- merge(obs2, sl[, .(subcluster, n_cells)], by = "subcluster", all = TRUE)
  bad  <- chk2[is.na(n_obj) | is.na(n_cells) | n_obj != n_cells]
  if (nrow(bad)) die(bad, basename(SUB_F), " does not match ", scol, " (rerun of 03 or 04c?)")
  si  <- match(subcl, sl$subcluster)
  has <- !is.na(si)
  lineage[has]  <- sl$lineage[si[has]]
  celltype[has] <- sl$celltype_label[si[has]]
}

flag <- if ("cell_flag" %in% names(md)) as.character(md$cell_flag) else rep("", nrow(md))
flag[is.na(flag)] <- ""
excl <- if ("exclude_from_de_lr" %in% names(md)) as.logical(md$exclude_from_de_lr) else rep(FALSE, nrow(md))
excl[is.na(excl)] <- FALSE
fsrc <- ifelse(nzchar(flag), "v4 run 1 (per cell)", "")
n_before <- sum(nzchar(flag)); x_before <- sum(excl)
if (nrow(dec)) for (k in seq_len(nrow(dec))) {
  ii <- cl == dec$run2_cluster[k]
  if (nzchar(dec$cell_flag[k])) {
    fill <- ii & !nzchar(flag)
    flag[fill] <- dec$cell_flag[k]
    fsrc[fill] <- sprintf("v4 run 2 (cluster %s)", dec$run2_cluster[k])
  }
  if (dec$exclude_from_de_lr[k]) excl[ii] <- TRUE
}

why <- ifelse(action == "drop", sprintf("v4 run-2 call: drop (cluster %s, %s)", cl, src), "")
cells <- data.table(
  cell = rownames(md),
  sample_key = if ("sample_key" %in% names(md)) md$sample_key else NA_character_,
  dataset_id = if ("dataset_id" %in% names(md)) md$dataset_id else NA_character_,
  cluster = cl, subcluster = subcl, action, call_source = src,
  lineage_v4 = lineage, celltype_v4 = celltype,
  cell_flag = flag, flag_source = fsrc, exclude_from_de_lr = excl, drop_reason = why)
fwrite(cells, CELLS_OUT, sep = "\t")

summ <- cells[, .(n_cells = .N, n_samples = uniqueN(sample_key),
                  pct_flagged = round(100 * mean(nzchar(cell_flag)), 1),
                  pct_excluded_de_lr = round(100 * mean(exclude_from_de_lr), 1)),
              by = .(cluster, subcluster, action, call_source, lineage_v4, celltype_v4)]
summ <- summ[order(as.integer(cluster), subcluster)]
fwrite(summ, SUM_OUT, sep = "\t")
msg("wrote ", CELLS_OUT, " and ", SUM_OUT)

# ---- final object -------------------------------------------------------------
keep <- action == "keep"
sub$cluster_v4         <- cl
sub$subcluster_v4      <- subcl
sub$lineage_v4         <- lineage
sub$celltype_v4        <- celltype
sub$call_source        <- src
sub$cell_flag          <- flag
sub$flag_source        <- fsrc
sub$exclude_from_de_lr <- excl
fin <- subset(sub, cells = colnames(sub)[keep])
Idents(fin) <- "celltype_v4"
save_obj(fin, OBJ_OUT)
msg("wrote ", OBJ_OUT, ": ", ncol(fin), " of ", ncol(sub), " cells")

# ---- provenance ----------------------------------------------------------------
by_lin <- cells[action == "keep", .(cells = .N, flagged = sum(nzchar(cell_flag)), excluded_de_lr = sum(exclude_from_de_lr)),
                by = lineage_v4][order(-cells)]
by_flag <- cells[action == "keep" & nzchar(cell_flag), .(cells = .N), by = .(cell_flag, flag_source)][order(-cells)]
writeLines(c(
  sprintf("PROVENANCE - Step 5 apply cluster calls - %s %s%s", SUBSET, rc, if (DRAFT) "  ** DRAFT **" else ""),
  sprintf("generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "script   : rebuild_05_apply_calls.R",
  sprintf("object   : %s (read only)", IN),
  sprintf("calls    : %s (%d of %d clusters reviewed)", CALLS, sum(tab$reviewed), nrow(tab)),
  sprintf("decisions: %s (%d rows for %s)", DEC_F, nrow(dec), SUBSET),
  sprintf("subclusters: %s", if (any(split)) sprintf("%s (%s)", SUB_F, scol) else "none"),
  if (DRAFT) c("", "DRAFT: clusters without a reviewer_call use the draft `call`. These outputs are",
               "for looking only; the reviewed run writes to rebuild_v4/ itself.") else character(0),
  "",
  sprintf("cells: %d in, %d kept, %d dropped", nrow(cells), sum(keep), sum(!keep)),
  sprintf("flags: %d cells flagged before step 05, %d after; exclude_from_de_lr %d before, %d after (all cells)",
          n_before, sum(nzchar(flag)), x_before, sum(excl)),
  "", "PER CLUSTER", show(summ),
  "", "KEPT CELLS BY LINEAGE", show(by_lin),
  "", "KEPT CELLS BY FLAG", if (nrow(by_flag)) show(by_flag) else "  none",
  "", "New metadata columns in the final object: cluster_v4, subcluster_v4, lineage_v4,",
  "  celltype_v4 (also the Idents), call_source, flag_source; cell_flag and",
  "  exclude_from_de_lr are updated. Cells with exclude_from_de_lr stay out of DE",
  "  and ligand-receptor work."
), PROV)
msg("wrote ", PROV)
