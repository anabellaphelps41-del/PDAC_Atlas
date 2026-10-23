#!/bin/bash
#SBATCH --job-name=subset_compartments
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=300g
#SBATCH --time=12:00:00
#SBATCH --output=logs/subset_%j.out
#SBATCH --error=logs/subset_%j.err

# Submit with an explicit cluster/partition, e.g.
#   sbatch -M smp -p high-mem run_subset.sbatch
# Needs MORE memory than the inspect job: the atlas and the subset coexist
# in RAM during subset(). Budget roughly 1.5x the inspect job's peak.

set -euo pipefail

PROJ=/ix/tbruno/anp795/PDAC_Atlas
ATLAS="$PROJ/combined_9_datasets_hierarchical_annotation_seurat.rds"
OUTDIR="$PROJ/compartments"
CFG="$PROJ/compartments.tsv"

mkdir -p "$OUTDIR" logs

module purge
module load r/4.6.0   # ships Seurat 5.5.0 / SeuratObject 5.4.0

Rscript -e 'suppressPackageStartupMessages(library(Seurat)); cat("Seurat", as.character(packageVersion("Seurat")), "OK\n")' || {
  echo "ERROR: Seurat not available under module r/4.6.0 in batch context." >&2
  exit 1
}

echo "host:  $(hostname)"
echo "start: $(date)"

Rscript "$PROJ/02_subset_compartments.R" "$ATLAS" "$OUTDIR" "$CFG"

echo "end: $(date)"

# After this finishes, check peak memory actually used so you can right-size
# future jobs instead of guessing:
#   sacct -M smp -j $SLURM_JOB_ID --format=JobID,MaxRSS,Elapsed,State
