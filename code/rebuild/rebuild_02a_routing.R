#!/usr/bin/env Rscript
# =============================================================================
# rebuild_02a_routing.R   --   reviewed cluster calls -> per-cell routing table
#
# Input : code/calls/cluster_calls_<subset>_res0.6_reviewed.tsv   (TNK, stromal, Bcell)
#         <UP>/step4b_<subset>_res0.6_singler_percell.tsv         (cell -> v3 cluster)
# Output: <DIR>/cell_routing.tsv
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
  "dropped at step 3 for being < 150 cells) are reported by rebuild_02_subset.R."
), PROV)
msg("wrote ", PROV)
