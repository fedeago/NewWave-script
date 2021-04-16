#!/bin/bash
#SBATCH --job-name="nb_mem"
#SBATCH --mail-type=ALL
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=450G

collectl -i1:5 -sZ -oT --procfilt C"R" > mem_nb_200k.txt &
collectl_pid="$!"
Rscript /mnt/temp1/BICCN_data/memory/memory_nb.R
pkill -TERM "$collectl_pid"
sleep 3
