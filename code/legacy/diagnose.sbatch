#!/bin/bash
#SBATCH --job-name=diagnose
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --time=01:00:00
#SBATCH --output=/ihome/tbruno/anp795/logs/diagnose_%j.out
#SBATCH --error=/ihome/tbruno/anp795/logs/diagnose_%j.err
#SBATCH --mail-user=anp795@pitt.edu
#SBATCH --mail-type=END,FAIL
module purge
module load r/4.6.0
Rscript /ihome/tbruno/anp795/diagnose.R
