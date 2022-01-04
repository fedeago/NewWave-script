collectl -i1:20 -sZ -oT --procfilt C"R" > mem_matrix_10.txt &
collectl_pid="$!"
Rscript /path/to/file/13_MemoryMatrix.R
pkill -TERM "$collectl_pid"
sleep 3
