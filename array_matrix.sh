#!/bin/bash

#SBATCH --job-name="matrix"
#SBATCH --mail-type=ALL
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=490G
#SBATCH --array=1

trans="$SLURM_ARRAY_TASK_ID"
Rscript NewWave_matrix.R "$trans"
