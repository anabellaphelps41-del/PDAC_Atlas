#!/bin/bash
#SBATCH --job-name=rb01c_doublets
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --output=/ix/tbruno/anp795/PDAC_Atlas/logs/%x_%j.out
#SBATCH --error=/ix/tbruno/anp795/PDAC_Atlas/logs/%x_%j.err
#SBATCH --mail-user=anp795@pitt.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
module purge
module load r/4.6.0

Rscript /ix/tbruno/anp795/PDAC_Atlas/code/rebuild/rebuild_01c_doublets.R 
