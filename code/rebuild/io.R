# io.R -- object I/O shared by the rebuild scripts.
#
# /ix reads at ~27 MB/s, so file size is what loading costs. qs2 (zstd,
# multithreaded) stores the Seurat objects ~6x smaller than uncompressed .rds
# (stromal: 3.0 GB -> 0.47 GB) and loads them in seconds instead of minutes.
#
# read_obj() also accepts a legacy .rds with the same stem, so objects written
# before the switch still load. save_obj() writes to a temporary file and
# renames it into place, so a job that dies mid-write cannot corrupt the
# existing copy (steps 04/04b re-save their input).

IO_THREADS <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))

# Directories. Steps 00m/01/01c (sample manifest, harmonized counts, doublet
# calls) are unchanged since v3 and are read from rebuild_v3. Everything from
# 01d on (DecontX, routing, the five subsets and their steps 02b-04b) is
# written to rebuild_v4. rebuild_v3 stays frozen: code/calls refers to its
# res-0.6 cluster numbers.
ROOT_DIR <- "/ix/tbruno/anp795/PDAC_Atlas"
UP_DIR   <- Sys.getenv("PDAC_UPSTREAM_DIR", file.path(ROOT_DIR, "rebuild_v3"))
V4_DIR   <- Sys.getenv("PDAC_REBUILD_DIR",  file.path(ROOT_DIR, "rebuild_v4"))
CALLS_DIR <- file.path(ROOT_DIR, "code", "calls")
SUBSETS  <- c("TNK", "stromal", "Bcell", "myeloid", "endothelial")

# Cluster column for one resolution. FindClusters names it
# <assay>_snn_res.<r>: RNA_ in the v3 objects, decontX_ in v4.
res_col_of <- function(obj, res) {
  cand <- paste0(unique(c(SeuratObject::DefaultAssay(obj), "RNA")), "_snn_res.", res)
  hit  <- cand[cand %in% colnames(obj@meta.data)]
  if (!length(hit))
    stop(sprintf("no cluster column for resolution %s. Available: %s", res,
                 paste(grep("_snn_res\\.", colnames(obj@meta.data), value = TRUE),
                       collapse = ", ")), call. = FALSE)
  hit[1]
}

read_obj <- function(path) {
  if (!file.exists(path)) {
    legacy <- sub("\\.qs2$", ".rds", path)
    if (legacy == path || !file.exists(legacy)) stop("object not found: ", path)
    path <- legacy
  }
  if (grepl("\\.qs2$", path)) qs2::qs_read(path, nthreads = IO_THREADS) else readRDS(path)
}

save_obj <- function(x, path) {
  stopifnot(grepl("\\.qs2$", path))
  tmp <- paste0(path, ".tmp")
  qs2::qs_save(x, tmp, nthreads = IO_THREADS)
  if (!file.rename(tmp, path)) stop("could not move ", tmp, " to ", path)
  invisible(path)
}
