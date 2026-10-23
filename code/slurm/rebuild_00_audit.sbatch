#!/bin/bash
#SBATCH --job-name=rb00_audit
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=00:30:00
#SBATCH --output=/ihome/tbruno/anp795/logs/rebuild00_audit_%j.out
#SBATCH --error=/ihome/tbruno/anp795/logs/rebuild00_audit_%j.err
#SBATCH --mail-user=anp795@pitt.edu
#SBATCH --mail-type=END,FAIL

# Metadata only. The script drops the counts with rm()+gc() immediately after
# each readRDS, so peak memory is roughly one subset load: T/NK is 5.91 GB on
# disk. 40G is generous; 2 CPUs is all it needs. Small requests queue fastest,
# and this is the job standing between you and the Step 3 rerun.

module purge
module load r/4.6.0

Rscript /ihome/tbruno/anp795/rebuild_00_provenance_audit.R TNK stromal Bcell
