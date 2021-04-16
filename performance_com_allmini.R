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
path <- paste0("/mnt/federico/BICCN_data/performance/results/200_thousand/performance_com_allmini",args,".Rdata")
n_cell <- 2e+5
load(paste0("/mnt/federico/BICCN_data/Campionamento/BICCN_",n_cell,".Rdata"))
etichette<-as.factor(dati$label)

####################  NEWWAVE
time_newWave_comm_allmini<- system.time(res_newWave_commo_allminibatch <- newFit(dati, K=10,
                          X = "~batch", 
                      children = 10, n_gene_disp=100,
                n_gene_par = 100, n_cell_par = ncol(dati)/10))
  

nAIC_mini_all <- newAIC(res_newWave_commo_allminibatch,counts(dati))
nBIC_mini_all <-newBIC(res_newWave_commo_allminibatch,counts(dati))
ari_mini_all <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(newW(res_newWave_commo_allminibatch),length(unique(etichette)))$cluster)))
# sil_mini_all <- mean(silhouette(as.numeric(etichette),daisy(newW(res_newWave_commo_allminibatch)))[,3])

 
assign(paste("stat_com_allmini", args,sep=""), cbind(time_newWave_comm_allmini[3],
                                                     nAIC_mini_all,
                                                     nBIC_mini_all,
                                                     # sil_mini_all,
                                                     ari_mini_all))
filename <- paste("stat_com_allmini", args,sep="")
save(list=filename, file=path)


