#!/bin/bash

#SBATCH --job-name="zinbwave"
#SBATCH --mail-type=ALL
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=480G
#SBATCH --array=1-10


trans="$SLURM_ARRAY_TASK_ID"
Rscript performance_zinbwave.R "$trans"