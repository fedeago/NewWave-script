library(cluster)
library(NewWave)
library(mclust)
library(BiocParallel)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(microbenchmark)
library(umap)
library(scater)
library(scran)
library(BiocParallel)
library(BiocSingular)
library(Seurat)

getPalette = colorRampPalette(brewer.pal(8, "Set2"))
load("/home/federico/Scrivania/mRNAmix_qc.RData")
sce2_qc$batch = "CEL-seq2"
sce8_qc$batch = "Sort-seq"

sce2_qc1 <- logNormCounts(sce2_qc)
sce8_qc1 <- logNormCounts(sce8_qc)


hvg1<- getTopHVGs(sce2_qc1, n=4000)
hvg2<- getTopHVGs(sce8_qc1, n=4000)

tot_hvg = Reduce(intersect, list(hvg1,hvg2))
tot_hvg = sample(tot_hvg,1000)
sce2_qc <- sce2_qc1[tot_hvg,]
sce8_qc <- sce8_qc1[tot_hvg,]

sce2_qc$label = paste(sce2_qc$H2228_prop,sce2_qc$H1975_prop,sce2_qc$HCC827_prop,sep="_")
sce8_qc$label = paste(sce8_qc$H2228_prop,sce8_qc$H1975_prop,sce8_qc$HCC827_prop,sep="_")


sce_merge = SingleCellExperiment(assays=list(counts=cbind(counts(sce2_qc),
                                                          counts(sce8_qc))))
sce_merge$Batch <- c(sce2_qc$batch,sce8_qc$batch)
sce_merge$cell_type = c(sce2_qc$label,sce8_qc$label)

save(sce_merge,file = "/mnt/federico/BICCN_data/datiMix.Rdata")

sce_merge$Batch<-as.factor(sce_merge$Batch)
sce_merge$cell_type<-as.factor(sce_merge$cell_type)
etichette<-as.factor(sce_merge$cell_type)


finNew <- newFit(sce_merge, K=10, X = "~Batch", commondispersion = F,
                 children = 10,
                 n_gene_par = 100, n_cell_par = ncol(sce_merge)/10)
data.umap = umap(finNew@W)

new_res = data.frame(data.umap$layout,etichette,batch=sce_merge$Batch)
# save(new_res, file="/mnt/federico/BICCN_data/lineePerPlot.RData")

sce_merge <- logNormCounts(sce_merge)
finPCA<- runPCA(sce_merge,ncomponents=10,BSPARAM = RandomParam(),BPPARAM=MulticoreParam(10))

data.umap = umap(reducedDim(finPCA, "PCA"))

pca_res = data.frame(data.umap$layout,etichette,batch=sce_merge$Batch)
# save(pca_res, file="/mnt/federico/BICCN_data/lineePerPlotPCA.RData")


seurat = CreateSeuratObject(counts = counts(sce_merge))
pca = reducedDim(finPCA, "PCA")
new = finNew@W
rownames(new) = colnames(seurat)
seurat[["NewWave"]] <- CreateDimReducObject(embeddings = new, key = "NewWave", assay = DefaultAssay(seurat))
attr(pca,"dimnames") =NULL
attr(pca,"varExplained") =NULL
attr(pca,"percentVar") =NULL
attr(pca,"rotation") =NULL
attr(pca,"dimnames") =NULL

seurat <- FindNeighbors(seurat, dims = 1:10,reduction = "NewWave")
seurat <- FindClusters(seurat, resolution = 1)
rownames(pca) = colnames(seurat)
seurat[["PCA_"]] <- CreateDimReducObject(embeddings = pca, key = "PCA_", assay = DefaultAssay(seurat))

seurat2 <- FindNeighbors(seurat, dims = 1:10,reduction = "PCA_")
seurat2 <- FindClusters(seurat2, resolution = 1)

ari_NewWave <- adjustedRandIndex(as.numeric(etichette),Idents(seurat))
ari_PCA <- adjustedRandIndex(as.numeric(etichette),Idents(seurat2))






ari_NewWave <- mean(sapply(1,function(y)adjustedRandIndex(as.numeric(etichette),kmeans(finNew@W, 7)$cluster)))
ari_PCA <- mean(sapply(1,function(y)adjustedRandIndex(as.numeric(etichette),kmeans(pca, 7)$cluster)))

new_res$Batch = new_res$batch
new_res$Proportion = new_res$etichette
pca_res$Batch = pca_res$batch
pca_res$Proportion = pca_res$etichette
pNW_b = ggplot(new_res,aes(x=X1,y=X2,col=Batch))+geom_point()+theme(legend.title = element_text("Batch"))
leg <- get_legend(pNW_b)
pNW_b = ggplot(new_res,aes(x=X1,y=X2,col=Batch))+geom_point()+theme(legend.position = "none",legend.key.size = unit(1.1, "cm"),axis.title = element_blank())
ggsave(filename="Scrivania/umapNW.png",ggplot(new_res,aes(x=X1,y=X2,col=batch))+geom_point() + scale_fill_discrete(name = "Dose"))
pPCA_b = ggplot(pca_res,aes(x=X1,y=X2,col=Batch))+geom_point()+theme(legend.position = "none",axis.title = element_blank())
ggsave(filename="Scrivania/umapPCA.png",ggplot(pca_res,aes(x=X1,y=X2,col=batch))+geom_point() + scale_fill_discrete(name = "Dose"))


pNW_e = ggplot(new_res,aes(x=X1,y=X2,col=Proportion))+geom_point()+theme(legend.title = element_text("Cell type"))+scale_color_manual(values=getPalette(20))
leg1 <- get_legend(pNW_e)
pNW_e = ggplot(new_res,aes(x=X1,y=X2,col=Proportion))+geom_point()+theme(legend.position = "none",legend.key.size = unit(1.1, "cm"),axis.title = element_blank())+scale_color_manual(values=getPalette(20))
ggsave(filename="Scrivania/umapNWeti.png",ggplot(new_res,aes(x=X1,y=X2,col=etichette))+geom_point() + scale_fill_discrete(name = "Dose"))
pPCA_e = ggplot(pca_res,aes(x=X1,y=X2,col=Proportion))+geom_point()+theme(legend.position = "none",axis.title = element_blank())+scale_color_manual(values=getPalette(20))
ggsave(filename="Scrivania/umapPCAeti.png",ggplot(pca_res,aes(x=X1,y=X2,col=etichette))+geom_point() + scale_fill_discrete(name = "Dose"))


fin <- ggarrange(pNW_b,pPCA_b, leg,pNW_e,pPCA_e,leg1,  labels = c("A", "B","", "C","D"),   
                 ncol = 3, # Second row with 2 plots in 2 different columns
                 nrow = 2, widths = c(2.5,2.5,1), font.label= list(size=8)) 
ggsave(filename="Scrivania/allRevisionMixture.png",fin,width = 9,  height = 10)


