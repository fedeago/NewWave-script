#!/bin/bash
#SBATCH --job-name="newRAM"
#SBATCH --mail-type=ALL
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=400G

collectl -i1:5 -sZ -oT --procfilt C"R" > mem_new_50k.txt &
collectl_pid="$!"
Rscript /mnt/temp1/BICCN_data/memory/memory_new.R
pkill -TERM "$collectl_pid"
sleep 3
