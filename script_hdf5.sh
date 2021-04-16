#!/bin/bash
#SBATCH --job-name="hdf5RAM"
#SBATCH --mail-type=ALL
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=60G


collectl -i1:20 -sZ -oT --procfilt C"R" > mem_hdf5_10_40.txt &
collectl_pid="$!"
Rscript /mnt/federico/BICCN_data/hdf5_example/memory_hdf5.R
pkill -TERM "$collectl_pid"
sleep 3
