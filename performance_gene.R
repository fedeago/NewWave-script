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
path <- paste0("/mnt/temp1/BICCN_data/performance/results/30_thousand/performance_gene",args,".Rdata")
n_cell <- 30000
load(paste0("/mnt/temp1/BICCN_data/Campionamento/BICCN_",n_cell,".Rdata"))


etichette<-as.factor(dati$label)


####################  NEWWAVE

time_newWave_gene <- system.time(res_newWave_genewise <- newFit(dati, K=10,
                                                                            X = "~batch", commondispersion = F, 
                                                                            children = 10))


nAIC_gene <- newAIC(res_newWave_genewise,counts(dati))
nBIC_gene <-newBIC(res_newWave_genewise,counts(dati))
ari_gene <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(newW(res_newWave_genewise),length(unique(etichette)))$cluster)))
sil_gene <- mean(silhouette(as.numeric(etichette),daisy(newW(res_newWave_genewise)))[,3])


assign(paste("stat_gene", args,sep=""), cbind(time_newWave_gene[3],
                           nAIC_gene,
                           nBIC_gene,
                           sil_gene,
                           ari_gene))

filename <- paste("stat_gene", args,sep="")
save(list=filename, file=path)
