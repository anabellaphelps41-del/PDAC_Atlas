#!/bin/bash
#SBATCH --job-name=rb00_audit
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=00:30:00
#SBATCH --output=/ix/tbruno/anp795/PDAC_Atlas/logs/%x_%j.out
#SBATCH --error=/ix/tbruno/anp795/PDAC_Atlas/logs/%x_%j.err
#SBATCH --mail-user=anp795@pitt.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
module purge
module load r/4.6.0

Rscript /ix/tbruno/anp795/PDAC_Atlas/code/rebuild/rebuild_00_provenance_audit.R TNK stromal Bcell
