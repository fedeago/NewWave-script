#!/bin/bash

#SBATCH --job-name="std"
#SBATCH --mail-type=ALL
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --array=1-10


trans="$SLURM_ARRAY_TASK_ID"
Rscript performance_std.R "$trans"
