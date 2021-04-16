#!/bin/bash

#SBATCH --job-name="g_r_300"
#SBATCH --mail-type=ALL
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --array=1-10

trans="$SLURM_ARRAY_TASK_ID"
Rscript performance_gene_parmini_random.R "$trans" 310515
