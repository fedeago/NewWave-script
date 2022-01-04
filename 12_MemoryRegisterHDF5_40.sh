collectl -i1:20 -sZ -oT --procfilt C"R" > mem_hdf5_10_40.txt &
collectl_pid="$!"
Rscript /path/to/file/11_MemoryHDF5_40.R
pkill -TERM "$collectl_pid"
sleep 3
