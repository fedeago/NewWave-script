collectl -i1:10 -sZ -oT --procfilt C"R" > mem_hdf5_10_10.txt &
collectl_pid="$!"
Rscript /path/to/file/09_MemoryHDF5.R
pkill -TERM "$collectl_pid"
sleep 3
