args = commandArgs(trailingOnly=TRUE)

library(cluster)
library(NewWave)
library(mclust)
library(BiocParallel)
library(microbenchmark)

# load("/mnt/temp1/BICCN_data/1k.Rdata")
# load("/mnt/temp1/BICCN_data/hvg_1e+05.Rdata")
# dati <-dati[hvg,]
#310515
path <- paste0("/mnt/federico/BICCN_data/performance/results/30_thousand/performance_std",args,".Rdata")
n_cell <- 30000
load(paste0("/mnt/federico/BICCN_data/Campionamento/BICCN_",n_cell,".Rdata"))


etichette<-as.factor(dati$label)
####################  NEWWAVE


newwave_commondisp_time <- system.time(res_newWave_commondisp <- newFit(dati,
                                                                           K=10,
                                                                           X = "~batch",  children = 10))
nAIC <- newAIC(res_newWave_commondisp, counts(dati))
nBIC <- newBIC(res_newWave_commondisp, counts(dati))
ari <-  mean(sapply(1:100, function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_commondisp),length(unique(etichette)))$cluster)))
# sil <- mean(silhouette(as.numeric(etichette),daisy(newW(res_newWave_commondisp)))[,3])

assign(paste("stat_std", args,sep=""), cbind(newwave_commondisp_time[3],
                           nAIC,
                           nBIC,
                           # sil,
                           ari))

filename <- paste("stat_std", args,sep="")
save(list=filename, file=path)
