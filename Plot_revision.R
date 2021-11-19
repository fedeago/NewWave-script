library(cluster)
library(NewWave)
library(mclust)
library(BiocParallel)
library(microbenchmark)
library(umap)


load(paste0("/mnt/federico/BICCN_data/Campionamento/BICCN_310515.Rdata"))
etichette<-as.factor(dati$label)

####################  NEWWAVE
res_newWave_genewise_allminibatch <- newFit(dati, K=10,
                                            X = "~batch", commondispersion = F, 
                                            children = 10,
                                            n_gene_par = 100, n_cell_par = ncol(dati)/10)
data.umap = umap(res_newWave_genewise_allminibatch@W)

new_res = data.frame(data.umap$layout,etichette,batch=dati$batch)
save(new_res, file="/mnt/federico/BICCN_data/modelloPerPlot.RData")
