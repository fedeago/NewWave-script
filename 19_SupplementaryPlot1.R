library(cluster)
library(NewWave)
library(mclust)
library(BiocParallel)
library(microbenchmark)
library(scater)
library(BiocSingular)
library(umap)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(here)

load(here("BICCN_hvg.Rdata"))
dati = all_data
etichette<-as.factor(dati$label)

####################  NEWWAVE
res_newWave_genewise_allminibatch <- newFit(dati, K=10,
                                            X = "~batch", commondispersion = F, 
                                            children = 10,
                                            n_gene_par = 100, n_cell_par = ncol(dati)/10)
data.umap = umap(res_newWave_genewise_allminibatch@W)

new_res = data.frame(data.umap$layout,etichette,batch=dati$batch)

dati <- logNormCounts(dati)
res_pca <- runPCA(dati,ncomponents=10,BSPARAM = RandomParam(),BPPARAM=MulticoreParam(10))

data.umap = umap(reducedDim(res_pca, "PCA"))

pca_res = data.frame(data.umap$layout,etichette,batch=dati$batch)


getPalette = colorRampPalette(brewer.pal(8, "Set2"))
new_res$Cell = as.factor(new_res$etichette)
pca_res$Cell = as.factor(pca_res$etichette)

pNW_b = ggplot(new_res,aes(x=X1,y=X2,col=Batch))+geom_point()+theme(legend.title = element_text("Batch"))
leg <- get_legend(pNW_b)
pNW_b = ggplot(new_res,aes(x=X1,y=X2,col=Batch))+geom_point()+theme(legend.position = "none",legend.key.size = unit(1.1, "cm"),axis.title = element_blank())
pNW_e = ggplot(new_res,aes(x=X1,y=X2,col=Cell))+geom_point()+theme(legend.title = element_text("Cell type"))+scale_color_manual(values=getPalette(20))
leg1 <- get_legend(pNW_e)
pNW_e = ggplot(new_res,aes(x=X1,y=X2,col=Cell))+geom_point()+theme(legend.position = "none",legend.key.size = unit(1.1, "cm"),axis.title = element_blank())+scale_color_manual(values=getPalette(20))
pPCA_b = ggplot(pca_res,aes(x=X1,y=X2,col=Batch))+geom_point()+theme(legend.position = "none",axis.title = element_blank())
pPCA_e = ggplot(pca_res,aes(x=X1,y=X2,col=Cell))+geom_point()+theme(legend.position = "none",axis.title = element_blank())+scale_color_manual(values=getPalette(20))


SuppPlot1 <- ggarrange(pNW_b,pPCA_b, leg,pNW_e,pPCA_e,leg1,  labels = c("A", "B","", "C","D"),   
                 ncol = 3, # Second row with 2 plots in 2 different columns
                 nrow = 2, widths = c(2.5,2.5,1), font.label= list(size=8))
SuppPlot1