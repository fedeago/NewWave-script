#!/bin/bash

#SBATCH --job-name="hdf5random"
#SBATCH --mail-type=ALL
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --array=1-5

trans="$SLURM_ARRAY_TASK_ID"
Rscript NewWave_hdf5_random.R "$trans"
