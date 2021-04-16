cache_dir <- file.path(tempdir(), "ExperimentHub")
dir.create(cache_dir)
ExperimentHub::setExperimentHubOption("CACHE", cache_dir)
system(paste("cp -L /mnt/federico/sassari/*", cache_dir))


library(cluster)
library(mclust)
library(BiocParallel)
library(TENxBrainData)
library(DelayedMatrixStats)
library(BiocSingular)

tenx <- TENxBrainData()

load("/mnt/federico/BICCN_data/hdf5_example/10x_hvg.Rdata")
n_cell <- 600000
path <- paste0("/mnt/federico/BICCN_data/hdf5_example/res/performance_pca","_",n_cell,".Rdata")

# hvg <- rowVars(counts(tenx))
# names(hvg) <- rowRanges(tenx)@elementMetadata$Symbol
# sort_hvg <- sort(hvg,decreasing=T)
# most_valuable <- sort_hvg[1:1000]
# hvg <- hvg %in% most_valuable
# names(hvg) <- rowRanges(tenx)@elementMetadata$Symbol
# save(hvg, file = "/mnt/federico/BICCN_data/10x_hvg.Rdata")
set.seed(1234)
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]


time_newWave_genemini_all <- system.time(res <- BiocSingular::runPCA(counts(dati), rank = 10, BPPARAM = BiocParallel::MulticoreParam(10), BSPARAM=ExactParam()))



stat_pca <- cbind(time_newWave_genemini_all[3])


save(stat_pca, file=path)