#!/bin/bash
#SBATCH --job-name="matrixRAM"
#SBATCH --mail-type=ALL
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120G
#SBATCH --nodelist=xen3
collectl -i1:10 -sZ -oT --procfilt C"R" > mem_matrix_10k.txt &
collectl_pid="$!"
Rscript /mnt/federico/BICCN_data/hdf5_example/memory_matrix.R
pkill -TERM "$collectl_pid"
sleep 3
