
cache_dir <- file.path(tempdir(), "ExperimentHub")
dir.create(cache_dir)
ExperimentHub::setExperimentHubOption("CACHE", cache_dir)
system(paste("cp -L /mnt/federico/sassari/*", cache_dir))

library(NewWave)
library(SingleCellExperiment)
library(TENxBrainData)


tenx <- TENxBrainData()
n_cell <- 10000
load("/mnt/federico/BICCN_data/hdf5_example/10x_hvg.Rdata")

set.seed(1234)
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]


res_newWave_genewise_allminibatch <- newFit(dati, K=10,
                                          commondispersion = F,
                                         children = 40,
                                 n_gene_par = 100, n_cell_par = ncol(dati)/10)
