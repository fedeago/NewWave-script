#!/bin/bash

#SBATCH --job-name="common"
#SBATCH --mail-type=ALL
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=115G
#SBATCH --array=1-10

trans="$SLURM_ARRAY_TASK_ID"
Rscript performance_com_allmini.R "$trans"
