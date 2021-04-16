args = commandArgs(trailingOnly=TRUE)
library(cluster)
library(NewWave)
library(zinbwave)
library(BiocParallel)
library(microbenchmark)
library(mclust)
# load("/mnt/temp1/BICCN_data/1k.Rdata")
# load("/mnt/temp1/BICCN_data/Campionamento/BICCN_data.Rdata")
# load("/mnt/temp1/BICCN_data/hvg_1e+05.Rdata")



path <- paste0("/mnt/temp1/BICCN_data/performance/results/10_thousand/performance_nbwave",args,".Rdata")
n_cell <- 10000
load(paste0("/mnt/temp1/BICCN_data/Campionamento/BICCN_",n_cell,".Rdata"))
etichette<-as.factor(dati$label)

####################  ZINBWAVE


time_zinb <- system.time(res_zinbwave <- zinbFit(dati,
                                                    X = "~batch", K=10, 
                                                    BPPARAM = MulticoreParam(10), zeroinflation=T))

zAIC <- zinbAIC(res_zinbwave,t(counts(dati)))
zBIC <- zinbBIC(res_zinbwave, t(counts(dati)))
ari<- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(res_zinbwave@W,length(unique(etichette)))$cluster)))
sil <- mean(silhouette(as.numeric(etichette),daisy(res_zinbwave@W))[,3])

assign(paste("stat_nbwave", args,sep=""), cbind(time_nb[3],
                                                zAIC,
                                                zBIC,
                                                sil,
                                                ari))
filename <- paste("stat_nbwave", args,sep="")
save(list=filename, file=path)

