library(cluster)
library(NewWave)
library(mclust)
library(BiocParallel)
library(microbenchmark)
library(scater)
library(BiocParallel)
library(BiocSingular)
library(umap)
load(paste0("/mnt/federico/BICCN_data/Campionamento/BICCN_310515.Rdata"))
etichette<-as.factor(dati$label)

####################  NEWWAVE
dati <- logNormCounts(dati)
res_pca <- runPCA(dati,ncomponents=10,BSPARAM = RandomParam(),BPPARAM=MulticoreParam(10))

data.umap = umap(reducedDim(res_pca, "PCA"))

pca_res = data.frame(data.umap$layout,etichette,batch=dati$batch)
save(pca_res, file="/mnt/federico/BICCN_data/modelloPerPlotPCA.RData")
