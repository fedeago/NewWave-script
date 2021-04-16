#!/bin/bash
#SBATCH --job-name="pcaTIME"
#SBATCH --mail-type=ALL
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=200G


Rscript pca_hdf5.R "$trans"
