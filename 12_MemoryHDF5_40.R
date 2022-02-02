library(NewWave)
library(SingleCellExperiment)
library(TENxBrainData)


tenx <- TENxBrainData()
n_cell <- 10000
load(here("10x_hvg.Rdata"))

set.seed(1234)
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]

res_newWave_genewise_allminibatch <- newFit(dati, K=10,
                                          commondispersion = F,
                                         children = 40,
                                 n_gene_par = 100, n_cell_par = ncol(dati)/10)
