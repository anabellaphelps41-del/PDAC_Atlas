#!/bin/bash
#SBATCH --job-name=rb4_01d_decontx
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

# One BLAS/OpenMP thread per forked worker. Without this, 16 workers each
# start a full set of BLAS threads (01c lost hours to that).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# A test run passes PDAC_DECONTX_SAMPLES (see the script header).
Rscript /ix/tbruno/anp795/PDAC_Atlas/code/rebuild/rebuild_01d_decontx.R
