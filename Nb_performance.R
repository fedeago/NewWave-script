library(cluster)
library(NewWave)
library(zinbwave)
library(BiocParallel)
library(microbenchmark)
library(mclust)
# load("/mnt/temp1/BICCN_data/1k.Rdata")
# load("/mnt/temp1/BICCN_data/Campionamento/BICCN_data.Rdata")
# load("/mnt/temp1/BICCN_data/hvg_1e+05.Rdata")



path <- "/mnt/federico/BICCN_data/performance/nb_time_performance_300k.Rdata"
n_cell <- 310515
load(paste0("/mnt/federico/BICCN_data/Campionamento/BICCN_",n_cell,".Rdata"))



etichette<-as.factor(dati$label)

####################  ZINBWAVE

gc()
time_zinb <- microbenchmark(res_zinbwave <- zinbFit(dati,
                                                    X = "~batch", K=10, 
                                                    BPPARAM = MulticoreParam(10), zeroinflation=F), times = 1)

zAIC <- zinbAIC(res_zinbwave,t(counts(dati)))
zBIC <- zinbBIC(res_zinbwave, t(counts(dati)))
# sil <- mean(silhouette(as.numeric(etichette),daisy(res_zinbwave@W))[,3])
ari<- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(res_zinbwave@W,length(unique(etichette)))$cluster)))

stat_zinb <- rbind(zAIC, zBIC,
                   # sil,
                   ari)

save(time_zinb, stat_zinb, file=path)
