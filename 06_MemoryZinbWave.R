library(zinbwave)
library(SingleCellExperiment)
library(BiocParallel)

load("/path/to/file/BICCN_hvg.Rdata")

n_cell <- 100000
set.seed(1234)
dati <-all_data[,sample(1:ncol(all_data),n_cell)]

res_zinbwave <- zinbFit(dati, X = "~batch", K=10, BPPARAM = MulticoreParam(10), zeroinflation=T)

