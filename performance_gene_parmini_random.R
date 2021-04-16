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

n_cell <- as.numeric(args[[2]])
args <- args[[1]]


path <- paste0("/mnt/federico/BICCN_data/performance/results/", n_cell/1000 ,"_thousand/performance_gene_parmini_random",args,".Rdata")

load(paste0("/mnt/federico/BICCN_data/Campionamento/BICCN_",n_cell,".Rdata"))

etichette<-as.factor(dati$label)

####################  NEWWAVE

time_newWave_genemini_all <- system.time(res_newWave_genewise_allminibatch <- newFit(dati, K=10,
                                X = "~batch", commondispersion = F, 
                       children = 10, random_start = T,
                     n_gene_par = 100, n_cell_par = ncol(dati)/10))
  

nAIC_genemini_all <- newAIC(res_newWave_genewise_allminibatch,counts(dati))
nBIC_genemini_all <- newBIC(res_newWave_genewise_allminibatch,counts(dati))
ari_genemini_all <- mean(sapply(1:100,function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_genewise_allminibatch), length(unique(etichette)))$cluster)))
# sil_genemini_all <- mean(silhouette(as.numeric(etichette),daisy(newW(res_newWave_genewise_allminibatch)))[,3])


assign(paste("stat_gene_parmini_random", args,sep=""), cbind(time_newWave_genemini_all[3],
                           nAIC_genemini_all,
                           nBIC_genemini_all,
                           # sil_genemini_all,
                           ari_genemini_all))
filename <- paste("stat_gene_parmini_random", args,sep="")
save(list=filename, file=path)
