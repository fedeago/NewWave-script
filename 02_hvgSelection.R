library(scran)
library(scater)
library(batchelor)
library(scran)
library(here)

num_cell <- 100000
load(here("BICCN_data.Rdata"))

set.seed(1234)
cell <- sample(1:ncol(all_data), num_cell)


dati <- all_data[,cell]

dati <- computeLibraryFactors(dati)
dati <- logNormCounts(dati)
logcounts(dati) <-as.matrix(assays(fastMNN(dati, batch = dati$batch, d=10), "reconstructed", withDimnames = F)[[1]])+1e6
dec <- modelGeneVar(dati)

hvg <- getTopHVGs(dec, n=1000)

name <- paste0("hvg.Rdata")

save(hvg, file = name)

all_data = all_data[hvg,]

save(all_data,file=here("BICCN_hvg.Rdata"))