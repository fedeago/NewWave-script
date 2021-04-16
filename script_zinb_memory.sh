#!/bin/bash
#SBATCH --job-name="zinb_mem"
#SBATCH --mail-type=ALL
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=400G


collectl -i1:5 -sZ -oT --procfilt C"R" > mem_zinb_50k.txt &
collectl_pid="$!"
Rscript /mnt/temp1/BICCN_data/memory/memory_zinb.R
pkill -TERM "$collectl_pid"
sleep 3
