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
