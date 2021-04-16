
cache_dir <- file.path(tempdir(), "ExperimentHub")
dir.create(cache_dir)
ExperimentHub::setExperimentHubOption("CACHE", cache_dir)
system(paste("cp -L /mnt/federico/sassari/*", cache_dir))

library(NewWave)
library(SingleCellExperiment)
library(TENxBrainData)
library(BiocSingular)

tenx <- TENxBrainData()
n_cell <- ncol(tenx)
load("/mnt/federico/BICCN_data/hdf5_example/10x_hvg.Rdata")

set.seed(1234)
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]


res <- BiocSingular::runPCA(counts(dati), rank = 10, BPPARAM = BiocParallel::MulticoreParam(10), BSPARAM=ExactParam())
                                            
