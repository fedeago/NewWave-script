collectl -i1:20 -sZ -oT --procfilt C"R" > mem_pca_10.txt &
collectl_pid="$!"
Rscript /path/to/file/15_MemoryPCA.R
pkill -TERM "$collectl_pid"
sleep 10

