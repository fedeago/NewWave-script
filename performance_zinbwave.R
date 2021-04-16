args = commandArgs(trailingOnly=TRUE)
library(cluster)
library(zinbwave)
library(mclust)
library(BiocParallel)
library(microbenchmark)

# load("/mnt/temp1/BICCN_data/1k.Rdata")
# load("/mnt/temp1/BICCN_data/hvg_1e+05.Rdata")
# dati <-dati[hvg,]

#310515
path <- paste0("/mnt/temp1/BICCN_data/performance/results/30_thousand/performance_zinbwave",args,".Rdata")
n_cell <- 30000
load(paste0("/mnt/temp1/BICCN_data/Campionamento/BICCN_",n_cell,".Rdata"))
etichette<-as.factor(dati$label)

####################  NEWWAVE
time_nbwave<- system.time(res_nbwave <- zinbwave(dati,X = "~batch", K=10,BPPARAM = MulticoreParam(10)))


nbAIC <- zinbAIC(res_nbwave,counts(dati))
nbBIC <-zinbBIC(res_nbwave,counts(dati))
ari <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(getW(res_nbwave),length(unique(etichette)))$cluster)))
sil <- mean(silhouette(as.numeric(etichette),daisy(newW(res_newWave_commo_allminibatch)))[,3])


assign(paste("stat_zinbwave", args,sep=""), cbind(time_nbwave[3],
                                                nbAIC,
                                                nbBIC,
                                                sil,
                                                ari))
filename <- paste("stat_zinbwave", args,sep="")
save(list=filename, file=path)