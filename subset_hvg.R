library(scran)
library(scater)
library(batchelor)
num_cell <- 100000
load("/mnt/temp1/BICCN_data/Campionamento/BICCN_data.Rdata")

set.seed(1234)
cell <- sample(1:ncol(all_data), num_cell)


dati <- all_data[,cell]
rm(all_data)
gc()
dati <- computeLibraryFactors(dati)
dati <- logNormCounts(dati)
logcounts(dati) <-as.matrix( assays(fastMNN(dati, batch = dati$batch, d=10), "reconstructed")[[1]])+1e6
dec <- modelGeneVar(dati)

hvg <- getTopHVGs(dec, n=1000)

name <- paste0("hvg_",num_cell,".Rdata")

save(hvg, file = name)
