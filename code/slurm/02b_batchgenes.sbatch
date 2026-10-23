#!/bin/bash
#SBATCH --job-name=rb02b_batchgenes
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=01:30:00
#SBATCH --output=/ix/tbruno/anp795/PDAC_Atlas/logs/%x_%j.out
#SBATCH --error=/ix/tbruno/anp795/PDAC_Atlas/logs/%x_%j.err
#SBATCH --mail-user=anp795@pitt.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
module purge
module load r/4.6.0

Rscript /ix/tbruno/anp795/PDAC_Atlas/code/rebuild/rebuild_02b_batchgenes.R 
