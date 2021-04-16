#!/bin/bash

#SBATCH --job-name="common"
#SBATCH --mail-type=ALL
#SBATCH --time=00:50:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --exclude=xen3

trans=1
Rscript performance_com_allmini.R "$trans"
