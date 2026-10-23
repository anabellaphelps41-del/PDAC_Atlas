#!/bin/bash
#SBATCH --job-name=rb04b_singler_stromal
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=03:00:00
#SBATCH --output=/ihome/tbruno/anp795/logs/rebuild04b_singler_stromal_%j.out
#SBATCH --error=/ihome/tbruno/anp795/logs/rebuild04b_singler_stromal_%j.err
#SBATCH --mail-user=anp795@pitt.edu
#SBATCH --mail-type=END,FAIL

# PREREQUISITE, on a LOGIN node, once:
#     module load r/4.6.0
#     Rscript ~/setup_singler_refs.R
# celldex pulls its references over the network and compute nodes have no
# outbound internet. Without the cache this job dies on the first celldex call.
#
# Sizing: rebuild_04 stromal (41,597 cells) ran in 4:47 at 13.8 GiB. SingleR
# per cell against two references is the expensive part -- roughly 15-30 min on
# 8 cores at this cell count. 60G / 3h is generous. The script reads
# SLURM_CPUS_PER_TASK for BiocParallel, so keep cpus-per-task and threads in step.
#
# For T/NK later (132,899 cells) raise to --mem=120G --time=08:00:00 and change
# the argument to TNK.

module purge
module load r/4.6.0

Rscript /ihome/tbruno/anp795/rebuild_04b_singler.R stromal 0.6
