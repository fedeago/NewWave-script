#!/bin/bash
#SBATCH --job-name="pcaRAM"
#SBATCH --mail-type=ALL
#SBATCH --time=07:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --nodelist=xen3

collectl -i1:20 -sZ -oT --procfilt C"R" > mem_pca_all.txt &
collectl_pid="$!"
Rscript /mnt/federico/BICCN_data/hdf5_example/memory_pca.R
pkill -TERM "$collectl_pid"
sleep 10

