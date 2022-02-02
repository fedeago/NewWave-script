collectl -i1:10 -sZ -oT --procfilt C"R" > NewWave_mem_10.txt &
collectl_pid="$!"
Rscript 04_MemoryNewWave.R
pkill -TERM "$collectl_pid"
sleep 3
