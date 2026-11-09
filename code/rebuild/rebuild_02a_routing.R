#!/usr/bin/env Rscript
# =============================================================================
# rebuild_02a_routing.R   --   reviewed cluster calls -> per-cell routing table
#
# Input : code/calls/cluster_calls_<subset>_res0.6_reviewed.tsv   (TNK, stromal, Bcell)
#         <UP>/step4b_<subset>_res0.6_singler_percell.tsv         (cell -> v3 cluster)
# Output: <DIR>/cell_routing.tsv
#         <DIR>/cell_calls_v4run1.tsv   (second round, see the end of the script)
#         <DIR>/PROVENANCE_02a_routing.txt
# (UP = rebuild_v3, DIR = rebuild_v4; see io.R)
#
# One row per cell of the three v3 subsets, carrying its cluster's decision:
#   action    keep | drop | route
#   to_subset the subset the cell belongs to in v4 (empty for drop)
#   lineage   lineage_reviewed from the tables (draft lineage, with the
#             reviewer's relabels); rebuild_01d_decontx.R uses it as the
#             DecontX cell-type label
#   followup  checks still owed after re-clustering (recheck_after_decontx,
#             recheck_dissociation, recheck_after_reintegration, marker_check,
#             merge_into:<cluster>)
# The decisions live in the reviewed TSVs (columns action / route_to /
# followup); this script only applies them by barcode, so they survive
# re-clustering. rebuild_02_subset.R reads the result.
#
# It stops if a table's clusters or n_cells do not match the v3 per-cell
# assignments: that would mean the tables describe a different clustering.
# =============================================================================

suppressPackageStartupMessages(library(data.table))
source("/ix/tbruno/anp795/PDAC_Atlas/code/rebuild/io.R")   # dirs

UP  <- UP_DIR
DIR <- V4_DIR
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)
RES <- "0.6"
OUT  <- file.path(DIR, "cell_routing.tsv")
PROV <- file.path(DIR, "PROVENANCE_02a_routing.txt")

msg <- function(...) { cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = ""); flush.console() }

ACTIONS <- c("keep", "drop", "route")
out <- list(); summ <- list()

for (s in c("TNK", "stromal", "Bcell")) {
  ft <- file.path(CALLS_DIR, sprintf("cluster_calls_%s_res%s_reviewed.tsv", s, RES))
  fc <- file.path(UP, sprintf("step4b_%s_res%s_singler_percell.tsv", s, RES))
  calls <- fread(ft, colClasses = list(character = c("route_to", "followup", "celltype_label", "lineage_reviewed")))
  cells <- fread(fc, select = c("cell", "cluster", "dataset_id"))

  need <- c("cluster", "n_cells", "action", "route_to", "followup", "celltype_label", "reviewer_call",
            "lineage_reviewed")
  if (length(setdiff(need, names(calls)))) stop(ft, " lacks: ", paste(setdiff(need, names(calls)), collapse = ", "))
  if (!all(calls$action %in% ACTIONS)) stop(ft, ": unknown action(s) ", paste(setdiff(calls$action, ACTIONS), collapse = ", "))
  bad_rt <- calls[(action == "route") != nzchar(route_to)]
  if (nrow(bad_rt)) stop(ft, ": route_to must be set exactly when action == route (clusters ",
                         paste(bad_rt$cluster, collapse = ", "), ")")
  if (!all(calls$route_to[nzchar(calls$route_to)] %in% SUBSETS)) stop(ft, ": route_to outside ", paste(SUBSETS, collapse = "/"))

  n_obs <- cells[, .N, by = cluster]
  chk   <- merge(calls[, .(cluster, n_cells)], n_obs, by = "cluster", all = TRUE)
  if (anyNA(chk) || any(chk$n_cells != chk$N)) {
    print(chk[is.na(n_cells) | is.na(N) | n_cells != N])
    stop(s, ": call table does not match the v3 res-", RES, " clustering")
  }

  r <- merge(cells, calls[, .(cluster, action, route_to, followup, celltype_label, reviewer_call,
                              lineage = lineage_reviewed)],
             by = "cluster", sort = FALSE)
  r[, `:=`(from_subset = s,
           to_subset = fifelse(action == "drop", "", fifelse(action == "route", route_to, s)))]
  out[[s]] <- r[, .(cell, dataset_id, from_subset, v3_cluster = cluster, action, to_subset,
                    lineage, celltype_label, followup, reviewer_call)]
  summ[[s]] <- out[[s]][, .(cells = .N), by = .(from_subset, action, to_subset)][order(action, to_subset)]
  msg(s, ": ", nrow(r), " cells, ", nrow(calls), " clusters")
}

rt <- rbindlist(out)
if (anyDuplicated(rt$cell)) stop("a cell appears in more than one v3 subset")
fwrite(rt, OUT, sep = "\t")
msg("wrote ", OUT, " (", nrow(rt), " cells)")

st <- rbindlist(summ)
print(st)
into <- rt[action != "drop", .(cells = .N), by = .(to_subset)][order(to_subset)]

show <- function(x) capture.output(print(x, row.names = FALSE))
writeLines(c(
  "PROVENANCE - per-cell routing from the reviewed res-0.6 cluster calls",
  paste("generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "script   : rebuild_02a_routing.R",
  sprintf("calls    : %s/cluster_calls_{TNK,stromal,Bcell}_res%s_reviewed.tsv", CALLS_DIR, RES),
  sprintf("clusters : %s/step4b_<subset>_res%s_singler_percell.tsv", UP, RES),
  paste("output   :", OUT),
  "",
  "CELLS BY DECISION", show(st), "",
  "CELLS EACH v4 SUBSET RECEIVES FROM THE v3 SUBSETS (myeloid/endothelial also",
  "get their own level2 cells in rebuild_02_subset.R)", show(into), "",
  "Cells of the v3 subsets that are not in these tables (GSE154778 B cells",
  "dropped at step 3 for being < 150 cells) are reported by rebuild_02_subset.R.",
  "",
  "SECOND ROUND (v4 run 1): see PROVENANCE section appended below, if present."
), PROV)
msg("wrote ", PROV)

# =============================================================================
# Second round: decisions on the v4 run-1 res-0.6 clusters
#   code/calls/v4_run1_decisions.tsv  (subset, run1_cluster, action, cell_flag,
#                                      exclude_from_de_lr, label, reason)
#   <DIR>/run1_snapshot/step4b_<subset>_res0.6_singler_percell.tsv (cell -> run-1 cluster)
# -> <DIR>/cell_calls_v4run1.tsv, one row per cell of a decided cluster.
#    rebuild_02_subset.R drops action == "drop" and carries cell_flag and
#    exclude_from_de_lr into the objects, so a flag follows the cell into
#    whatever cluster it joins after re-clustering.
# =============================================================================
DEC_F  <- file.path(CALLS_DIR, "v4_run1_decisions.tsv")
SNAP   <- file.path(DIR, "run1_snapshot")
OUT2   <- file.path(DIR, "cell_calls_v4run1.tsv")
if (file.exists(DEC_F)) {
  dec <- fread(DEC_F, colClasses = list(character = c("cell_flag", "label", "reason")))
  if (!all(dec$action %in% c("keep", "drop"))) stop(DEC_F, ": action must be keep or drop")
  if (!all(dec$subset %in% SUBSETS)) stop(DEC_F, ": unknown subset")
  if (anyDuplicated(dec[, .(subset, run1_cluster)])) stop(DEC_F, ": a cluster is listed twice")
  dec[, exclude_from_de_lr := as.logical(exclude_from_de_lr)]
  per <- rbindlist(lapply(unique(dec$subset), function(s) {
    f <- file.path(SNAP, sprintf("step4b_%s_res%s_singler_percell.tsv", s, RES))
    if (!file.exists(f)) stop("run-1 snapshot missing: ", f)
    x <- fread(f, select = c("cell", "cluster", "dataset_id"))
    d <- dec[subset == s]
    miss <- setdiff(d$run1_cluster, unique(x$cluster))
    if (length(miss)) stop(s, ": run-1 cluster(s) not in the snapshot: ", paste(miss, collapse = ", "))
    m <- merge(x, d, by.x = "cluster", by.y = "run1_cluster")
    m[, .(cell, dataset_id, subset = s, run1_cluster = cluster, action, cell_flag, exclude_from_de_lr, label)]
  }))
  if (anyDuplicated(per$cell)) stop("a cell is in two decided run-1 clusters")
  fwrite(per, OUT2, sep = "\t")
  s2 <- per[, .(cells = .N), by = .(subset, run1_cluster, action, cell_flag)][order(subset, run1_cluster)]
  print(s2)
  msg("wrote ", OUT2, " (", nrow(per), " cells)")
  cat(c("", "SECOND ROUND - decisions on v4 run-1 clusters",
        paste("decisions:", DEC_F), paste("clusters :", SNAP), paste("output   :", OUT2), "",
        show(s2)), file = PROV, sep = "\n", append = TRUE)
} else {
  msg("no ", DEC_F, "; second round skipped")
  if (file.exists(OUT2)) file.remove(OUT2)
}
