#!/bin/bash
#SBATCH --job-name=rebuild03_stromal
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=04:00:00
#SBATCH --output=/ihome/tbruno/anp795/logs/rebuild03_stromal_%j.out
#SBATCH --error=/ihome/tbruno/anp795/logs/rebuild03_stromal_%j.err
#SBATCH --mail-user=anp795@pitt.edu
#SBATCH --mail-type=END,FAIL

module purge
module load r/4.6.0

Rscript /ihome/tbruno/anp795/rebuild_03_pipeline.R stromal
