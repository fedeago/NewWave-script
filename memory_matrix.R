
cache_dir <- file.path(tempdir(), "ExperimentHub")
dir.create(cache_dir)
ExperimentHub::setExperimentHubOption("CACHE", cache_dir)
system(paste("cp -L /mnt/federico/sassari/*", cache_dir))
library(NewWave)
library(SingleCellExperiment)
library(TENxBrainData)
#310515


#n_cell <- 10000

tenx <- TENxBrainData()
load("/mnt/federico/BICCN_data/hdf5_example/10x_hvg.Rdata")


n_cell <- 10000
set.seed(1234)
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]

dati <- SingleCellExperiment::SingleCellExperiment(assays=list(counts = as.matrix(counts(dati))))

res_newWave_genewise_allminibatch <- newFit(dati, K=10,
                             commondispersion = F,
                   children = 10,
                   n_gene_par = 100, n_cell_par = ncol(dati)/10)
