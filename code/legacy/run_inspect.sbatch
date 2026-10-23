#!/bin/bash
#SBATCH --job-name=inspect_atlas
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=250g
#SBATCH --time=8:00:00
#SBATCH --output=logs/inspect_%j.out
#SBATCH --error=logs/inspect_%j.err

# NOTE: CRC runs a multi-cluster SLURM setup. Submit with an explicit cluster,
# e.g.  sbatch -M smp -p high-mem run_inspect.sbatch
# Check what exists and what the per-partition memory ceiling is BEFORE submitting:
#   sinfo -M all -o "%P %c %m %l"
#   crc-sinfo
# The htc partition may cap below 250 GB; smp's high-mem partition is the usual
# home for jobs this size.

set -euo pipefail

PROJ=/ix/tbruno/anp795/PDAC_Atlas
ATLAS="$PROJ/combined_9_datasets_hierarchical_annotation_seurat.rds"
OUTDIR="$PROJ/inspect"

mkdir -p "$OUTDIR" logs

module purge
module load r/4.6.0   # ships Seurat 5.5.0 / SeuratObject 5.4.0

# Fail fast if the batch environment can't see Seurat, rather than after a
# 40-minute readRDS. The RStudio Server app and a batch job do not always
# resolve the same library paths.
Rscript -e 'suppressPackageStartupMessages(library(Seurat)); cat("Seurat", as.character(packageVersion("Seurat")), "OK\n")' || {
  echo "ERROR: Seurat not available under module r/4.6.0 in batch context." >&2
  echo "Run in RStudio:  .libPaths()   and set R_LIBS_USER to match." >&2
  exit 1
}

echo "host:  $(hostname)"
echo "start: $(date)"
echo "mem allocated: ${SLURM_MEM_PER_NODE:-NA} MB"

Rscript "$PROJ/01_inspect_atlas.R" "$ATLAS" "$OUTDIR"

echo "end: $(date)"
