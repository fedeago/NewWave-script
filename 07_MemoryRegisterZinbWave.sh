collectl -i1:10 -sZ -oT --procfilt C"R" > zinbwave_mem_10.txt &
collectl_pid="$!"
Rscript /path/to/file/06_MemoryZinbWave.R
pkill -TERM "$collectl_pid"
sleep 3
