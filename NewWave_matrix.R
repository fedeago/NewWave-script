args = commandArgs(trailingOnly=TRUE)

cache_dir <- file.path(tempdir(), "ExperimentHub")
dir.create(cache_dir)
ExperimentHub::setExperimentHubOption("CACHE", cache_dir)
system(paste("cp -L /mnt/federico/sassari/*", cache_dir))

library(cluster)
library(NewWave)
library(mclust)
library(BiocParallel)
library(microbenchmark)
library(TENxBrainData)
library(DelayedMatrixStats)

tenx <- TENxBrainData()
load("/mnt/federico/BICCN_data/hdf5_example/10x_hvg.Rdata")

n_cell <- ncol(tenx)
path <- paste0("/mnt/federico/BICCN_data/hdf5_example/res/performance_matrix",args,"_",n_cell,".Rdata")
# hvg <- rowVars(counts(tenx))
# names(hvg) <- rowRanges(tenx)@elementMetadata$Symbol
# sort_hvg <- sort(hvg,decreasing=T)
# most_valuable <- sort_hvg[1:1000]
# hvg <- hvg %in% most_valuable
# names(hvg) <- rowRanges(tenx)@elementMetadata$Symbol
# save(hvg, file = "/mnt/federico/BICCN_data/10x_hvg.Rdata")
set.seed(1234)
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]

dati <- SingleCellExperiment::SingleCellExperiment(assays=list(counts = as.matrix(counts(dati))))

time_newWave_genemini_all <- system.time(res_newWave_genewise_allminibatch <- newFit(dati, K=10,
                                                                                      commondispersion = F,
                                                                                     children = 10,
                                                                                     n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_genemini_all <- newAIC(res_newWave_genewise_allminibatch,as.matrix(counts(dati)))
nBIC_genemini_all <- newBIC(res_newWave_genewise_allminibatch,as.matrix(counts(dati)))
# ari_genemini_all <- mean(sapply(1:100,function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_genewise_allminibatch), length(unique(etichette)))$cluster)))
# sil_genemini_all <- mean(silhouette(as.numeric(etichette),daisy(newW(res_newWave_genewise_allminibatch)))[,3])


assign(paste("stat_gene_parmini", args,sep=""), cbind(time_newWave_genemini_all[3],
                                                      nAIC_genemini_all,
                                                      nBIC_genemini_all))

filename <- paste("stat_gene_parmini", args,sep="")
save(list=filename, file=path)
