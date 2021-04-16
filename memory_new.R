library(NewWave)
library(SingleCellExperiment)

#310515
load("/mnt/temp1/BICCN_data/Campionamento/BICCN_data.Rdata")
load("/mnt/temp1/BICCN_data/hvg_1e+05.Rdata")

n_cell <- 50000
set.seed(1234)
dati <-all_data[hvg,sample(1:ncol(all_data),n_cell)]
rm(all_data)
gc()

res_zinbwave <- newWave(dati,
                         X = "~batch", K=10, children= 10)

