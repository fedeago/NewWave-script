library(NewWave)
library(SingleCellExperiment)
library(TENxBrainData)
library(BiocSingular)

tenx <- TENxBrainData()
n_cell <- 10000
load(here("10x_hvg.Rdata"))

set.seed(1234)
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]

res <- BiocSingular::runPCA(counts(dati), rank = 10, BPPARAM = BiocParallel::MulticoreParam(10), BSPARAM=ExactParam())
                                            
