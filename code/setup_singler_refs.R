#!/usr/bin/env Rscript
# =============================================================================
# setup_singler_refs.R    ***  RUN THIS ON A LOGIN NODE, NOT VIA SBATCH  ***
#
# celldex fetches its reference datasets over the network from ExperimentHub.
# Compute nodes on HTC have no outbound internet (rebuild_03_pipeline.R lines
# 34-43 already document this for package installs; it applies to data too).
# A SingleR job submitted without this step will die on the first celldex call.
#
# This caches both references as plain .rds under DIR/ref/ so the compute node
# only ever does readRDS().
#
#   module load r/4.6.0
#   Rscript ~/setup_singler_refs.R
#
# Installs, if missing: SingleR, celldex, SummarizedExperiment (Bioconductor).
# Expect a few minutes on first run.
# =============================================================================

DIR     <- "/ix/tbruno/anp795/PDAC_Atlas"
REF_DIR <- "/ix/tbruno/anp795/PDAC_Atlas/ref"   # shared, not per rebuild
dir.create(REF_DIR, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

need <- c("SingleR", "celldex", "SummarizedExperiment", "BiocParallel")
miss <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
if (length(miss)) {
  msg("installing: ", paste(miss, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(miss, ask = FALSE, update = FALSE)
}
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE))
    stop("still missing after install: ", p, call. = FALSE)
}
msg("all packages available")

# ---------------------------------------------------------------------------
# Two references, deliberately. They disagree in useful ways and the gate uses
# agreement as its confidence measure.
#
#   HumanPrimaryCellAtlas - broad, has Fibroblasts, Smooth_muscle_cells,
#     Chondrocytes, Neurons/Astrocytes (for the Schwann cluster), Epithelial_cells.
#     Best coverage of the non-immune stroma.
#   BlueprintEncode - immune-heavy and better calibrated for myeloid and
#     lymphoid, which is exactly where the fibroblast panel is failing.
# ---------------------------------------------------------------------------
# celldex changed its API at 1.99: the named accessor functions were replaced by
# fetchReference(<name>, <version>). R 4.6.0 ships with a recent Bioconductor, so
# the legacy functions may be deprecated or gone. Try legacy first, fall back to
# fetchReference, and report which path was taken.
msg("celldex version: ", as.character(utils::packageVersion("celldex")))

fetch <- function(legacy_name, new_name, out_name) {
  f <- file.path(REF_DIR, paste0(out_name, ".rds"))
  if (file.exists(f)) { msg("already cached: ", f); return(invisible(NULL)) }

  ref <- NULL
  if (legacy_name %in% getNamespaceExports("celldex")) {
    msg("fetching ", out_name, " via celldex::", legacy_name, "() ...")
    ref <- tryCatch(getExportedValue("celldex", legacy_name)(),
                    error = function(e) { msg("  legacy call failed: ", conditionMessage(e)); NULL })
  }
  if (is.null(ref) && "fetchReference" %in% getNamespaceExports("celldex")) {
    vers <- tryCatch(celldex::surveyReferences()$version[
                       celldex::surveyReferences()$name == new_name],
                     error = function(e) character(0))
    vers <- if (length(vers)) max(vers) else "2024-02-26"
    msg("fetching ", out_name, " via celldex::fetchReference('", new_name, "', '", vers, "') ...")
    ref <- tryCatch(celldex::fetchReference(new_name, vers),
                    error = function(e) { msg("  fetchReference failed: ", conditionMessage(e)); NULL })
  }
  if (is.null(ref))
    stop(sprintf(paste0("could not obtain reference '%s'.\n",
                        "  Check celldex docs for this version, or tell me the celldex\n",
                        "  version and the output of celldex::surveyReferences()."),
                 out_name), call. = FALSE)

  if (!"label.main" %in% colnames(SummarizedExperiment::colData(ref)))
    stop(sprintf("reference '%s' has no label.main column; columns are: %s",
                 out_name, paste(colnames(SummarizedExperiment::colData(ref)), collapse = ", ")),
         call. = FALSE)

  msg("  ", nrow(ref), " genes x ", ncol(ref), " samples")
  msg("  main labels (", length(unique(ref$label.main)), "): ",
      paste(sort(unique(ref$label.main)), collapse = ", "))
  saveRDS(ref, f)
  msg("wrote ", f, "  (", sprintf("%.1f MB", file.size(f) / 1e6), ")")
}

fetch("HumanPrimaryCellAtlasData", "hpca",             "celldex_HPCA")
fetch("BlueprintEncodeData",       "blueprint_encode", "celldex_BlueprintEncode")

msg("done. Compute nodes can now readRDS these; no network needed.")
msg("Next: sbatch ~/rebuild_04b_singler.sbatch")
