library(NewWave)
library(SingleCellExperiment)
library(TENxBrainData)
library(DelayedMatrixStats)
library(BiocSingular)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(grid)
library(ggpubr)
library(here)

tenx <- TENxBrainData()
n_cell <- 10000
load(here("10x_hvg.Rdata"))

set.seed(1234)
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]
datim <- SingleCellExperiment::SingleCellExperiment(assays=list(counts = as.matrix(counts(dati))))

timeNEWHDF510_10 <- system.time( newFit(dati, K=10,
                                                                                     commondispersion = F,
                                                                                     children = 10,
                                                                                     n_gene_par = 100, n_cell_par = ncol(dati)/10))
timeNEWHDF540_10 <- system.time( newFit(dati, K=10,
                                                                                     commondispersion = F,
                                                                                     children = 40,
                                                                                     n_gene_par = 100, n_cell_par = ncol(dati)/10))

timePCAHDF510_10 <- system.time( BiocSingular::runPCA(counts(dati), rank = 10, BPPARAM = BiocParallel::MulticoreParam(10), BSPARAM=ExactParam()))
timeNEWMATRIX10_10 <- system.time( newFit(datim, K=10,
                                                                                     commondispersion = F,
                                                                                     children = 10,
                                                                                     n_gene_par = 100, n_cell_par = ncol(dati)/10))
stat10 <- rbind(cbind(timeNEWHDF510_10[3],"10000","HDF5_10cores"),cbind(timeNEWHDF540_10[3],"10000","HDF5_40cores"),cbind(timePCAHDF510_10[3],"10000","PCA_10cores"),cbind(timeNEWMATRIX10_10[3],"10000","Matrix_10cores"))

n_cell <- 100000
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]
datim <- SingleCellExperiment::SingleCellExperiment(assays=list(counts = as.matrix(counts(dati))))
timeNEWHDF510_100 <- system.time( newFit(dati, K=10,
                                                                                     commondispersion = F,
                                                                                     children = 10,
                                                                                     n_gene_par = 100, n_cell_par = ncol(dati)/10))
timeNEWHDF540_100 <- system.time( newFit(dati, K=10,
                                                                            commondispersion = F,
                                                                            children = 40,
                                                                            n_gene_par = 100, n_cell_par = ncol(dati)/10))

timePCAHDF510_100 <- system.time( BiocSingular::runPCA(counts(dati), rank = 10, BPPARAM = BiocParallel::MulticoreParam(10), BSPARAM=ExactParam()))
timeNEWMATRIX10_100 <- system.time( newFit(datim, K=10,
                                                                              commondispersion = F,
                                                                              children = 10,
                                                                              n_gene_par = 100, n_cell_par = ncol(dati)/10))
stat100 <- rbind(cbind(timeNEWHDF510_100[3],"100000","HDF5_10cores"),cbind(timeNEWHDF540_100[3],"100000","HDF5_40cores"),cbind(timePCAHDF510_100[3],"100000","PCA_10cores"),cbind(timeNEWMATRIX10_100[3],"100000","Matrix_10cores"))
n_cell <- 200000
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]
datim <- SingleCellExperiment::SingleCellExperiment(assays=list(counts = as.matrix(counts(dati))))
timeNEWHDF510_200 <- system.time( newFit(dati, K=10,
                                                                                     commondispersion = F,
                                                                                     children = 10,
                                                                                     n_gene_par = 100, n_cell_par = ncol(dati)/10))
timeNEWHDF540_200 <- system.time( newFit(dati, K=10,
                                                                            commondispersion = F,
                                                                            children = 40,
                                                                            n_gene_par = 100, n_cell_par = ncol(dati)/10))

timePCAHDF510_200 <- system.time( BiocSingular::runPCA(counts(dati), rank = 10, BPPARAM = BiocParallel::MulticoreParam(10), BSPARAM=ExactParam()))
timeNEWMATRIX10_200 <- system.time( newFit(datim, K=10,
                                                                              commondispersion = F,
                                                                              children = 10,
                                                                              n_gene_par = 100, n_cell_par = ncol(dati)/10))
stat200 <- rbind(cbind(timeNEWHDF510_200[3],"200000","HDF5_10cores"),cbind(timeNEWHDF540_200[3],"200000","HDF5_40cores"),cbind(timePCAHDF510_200[3],"200000","PCA_10cores"),cbind(timeNEWMATRIX10_200[3],"200000","Matrix_10cores"))
n_cell <- 300000
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]
datim <- SingleCellExperiment::SingleCellExperiment(assays=list(counts = as.matrix(counts(dati))))
timeNEWHDF510_300 <- system.time( newFit(dati, K=10,
                                                                                     commondispersion = F,
                                                                                     children = 10,
                                                                                     n_gene_par = 100, n_cell_par = ncol(dati)/10))
timeNEWHDF540_300 <- system.time( newFit(dati, K=10,
                                                                            commondispersion = F,
                                                                            children = 40,
                                                                            n_gene_par = 100, n_cell_par = ncol(dati)/10))

timePCAHDF510_300 <- system.time( BiocSingular::runPCA(counts(dati), rank = 10, BPPARAM = BiocParallel::MulticoreParam(10), BSPARAM=ExactParam()))
timeNEWMATRIX10_300 <- system.time( newFit(datim, K=10,
                                                                              commondispersion = F,
                                                                              children = 10,
                                                                              n_gene_par = 100, n_cell_par = ncol(dati)/10))
stat300 <- rbind(cbind(timeNEWHDF510_300[3],"300000","HDF5_10cores"),cbind(timeNEWHDF540_300[3],"300000","HDF5_40cores"),cbind(timePCAHDF510_300[3],"300000","PCA_10cores"),cbind(timeNEWMATRIX10_300[3],"300000","Matrix_10cores"))
n_cell <- 600000
dati <- tenx[hvg,sample(ncol(counts(tenx)),n_cell)]
datim <- SingleCellExperiment::SingleCellExperiment(assays=list(counts = as.matrix(counts(dati))))
timeNEWHDF510_600 <- system.time( newFit(dati, K=10,
                                                                                     commondispersion = F,
                                                                                     children = 10,
                                                                                     n_gene_par = 100, n_cell_par = ncol(dati)/10))
timeNEWHDF540_600 <- system.time( newFit(dati, K=10,
                                                                            commondispersion = F,
                                                                            children = 40,
                                                                            n_gene_par = 100, n_cell_par = ncol(dati)/10))

timePCAHDF510_600 <- system.time( BiocSingular::runPCA(counts(dati), rank = 10, BPPARAM = BiocParallel::MulticoreParam(10), BSPARAM=ExactParam()))
timeNEWMATRIX10_600 <- system.time( newFit(datim, K=10,
                                                                              commondispersion = F,
                                                                              children = 10,
                                                                              n_gene_par = 100, n_cell_par = ncol(dati)/10))
stat600 <- rbind(cbind(timeNEWHDF510_600[3],"600000","HDF5_10cores"),cbind(timeNEWHDF540_600[3],"600000","HDF5_40cores"),cbind(timePCAHDF510_600[3],"600000","PCA_10cores"),cbind(timeNEWMATRIX10_600[3],"600000","Matrix_10cores"))
dati <- tenx[hvg,]
datim <- SingleCellExperiment::SingleCellExperiment(assays=list(counts = as.matrix(counts(dati))))
timeNEWHDF510_all <- system.time( newFit(dati, K=10,
                                                                                     commondispersion = F,
                                                                                     children = 10,
                                                                                     n_gene_par = 100, n_cell_par = ncol(dati)/10))
timeNEWHDF540_all <- system.time( newFit(dati, K=10,
                                                                            commondispersion = F,
                                                                            children = 40,
                                                                            n_gene_par = 100, n_cell_par = ncol(dati)/10))

timePCAHDF510_all <- system.time( BiocSingular::runPCA(counts(dati), rank = 10, BPPARAM = BiocParallel::MulticoreParam(10), BSPARAM=ExactParam()))
timeNEWMATRIX10_all <- system.time( newFit(datim, K=10,
                                                                              commondispersion = F,
                                                                              children = 10,
                                                                              n_gene_par = 100, n_cell_par = ncol(dati)/10))
statall <- rbind(cbind(timeNEWHDF510_all[3],"1300000","HDF5_10cores"),cbind(timeNEWHDF540_all[3],"1300000","HDF5_40cores"),cbind(timePCAHDF510_all[3],"1300000","PCA_10cores"),cbind(timeNEWMATRIX10_all[3],"1300000","Matrix_10cores"))

stat = rbind(stat10,stat100,stat200,stat300,stat600,stat1000)

stat = as.data.frame(stat)
stat$elapsed = as.numeric(stat$V1)
stat$model = as.factor(stat$V3)

time <- ggplot(stat, aes(y=elapsed/60,x=V2,group=model, col=model))+geom_point()+geom_line()+ylab("Time(in minute)")+xlab("Number of cell(in thousands)")+ labs(fill = "Methods")+theme(legend.position = "right", legend.key.size = unit(1.1, "cm"))+ scale_color_manual(values=c("#000000", "#ffcc00", "#ff6633","#99ff00"))
time

mem_h_10<- read.table("mem_hdf5_10_40.txt", header=TRUE, quote="\"")
h_10 <- mem_h_10[,c(1,9)]
h_10$RSS <- gsub("G", "000000", h_10$RSS)
h_10$RSS <- gsub("M","000", h_10$RSS)
h_10$RSS <- gsub("K","", h_10$RSS)
h_10$RSS <- as.numeric(h_10$RSS)
h_10$index <- as.numeric(as.factor(h_10$Time))
h_10 = h_10 %>% group_by(index) %>% summarise(x = sum(RSS))
h_10$model <-"hdf5_40cores"
h_10$cell <-10

mem_h_100<- read.table("mem_hdf5_100_40.txt", header=TRUE, quote="\"")
h_100 <- mem_h_100[,c(1,9)]
h_100$RSS <- gsub("G", "000000", h_100$RSS)
h_100$RSS <- gsub("M","000", h_100$RSS)
h_100$RSS <- gsub("K","", h_100$RSS)
h_100$RSS <- as.numeric(h_100$RSS)
h_100$index <- as.numeric(as.factor(h_100$Time))
h_100 = h_100 %>% group_by(index) %>% summarise(x = sum(RSS))
h_100$model <-"hdf5_40cores"
h_100$cell <-100

mem_h_200<- read.table("mem_hdf5_200_40.txt", header=TRUE, quote="\"")
h_200 <- mem_h_200[,c(1,9)]
h_200$RSS <- gsub("G", "000000", h_200$RSS)
h_200$RSS <- gsub("M","000", h_200$RSS)
h_200$RSS <- gsub("K","", h_200$RSS)
h_200$RSS <- as.numeric(h_200$RSS)
h_200$index <- as.numeric(as.factor(h_200$Time))
h_200 = h_200 %>% group_by(index) %>% summarise(x = sum(RSS))
h_200$model <-"hdf5_40cores"
h_200$cell <-200

mem_h_300<- read.table("mem_hdf5_300_40.txt", header=TRUE, quote="\"")
h_300 <- mem_h_300[,c(1,9)]
h_300$RSS <- gsub("G", "000000", h_300$RSS)
h_300$RSS <- gsub("M","000", h_300$RSS)
h_300$RSS <- gsub("K","", h_300$RSS)
h_300$RSS <- as.numeric(h_300$RSS)
h_300$index <- as.numeric(as.factor(h_300$Time))
h_300 = h_300 %>% group_by(index) %>% summarise(x = sum(RSS))
h_300$model <-"hdf5_40cores"
h_300$cell <-300

mem_h_600<- read.table("mem_hdf5_600_40.txt", header=TRUE, quote="\"")
h_600 <- mem_h_600[,c(1,9)]
h_600$RSS <- gsub("G", "000000", h_600$RSS)
h_600$RSS <- gsub("M","000", h_600$RSS)
h_600$RSS <- gsub("K","", h_600$RSS)
h_600$RSS <- as.numeric(h_600$RSS)
h_600$index <- as.numeric(as.factor(h_600$Time))
h_600 = h_600 %>% group_by(index) %>% summarise(x = sum(RSS))
h_600$model <-"hdf5_40cores"
h_600$cell <-600

mem_h_all<- read.table("mem_hdf5_1300_40.txt", header=TRUE, quote="\"")
h_all <- mem_h_all[,c(1,9)]
h_all$RSS <- gsub("G", "000000", h_all$RSS)
h_all$RSS <- gsub("M","000", h_all$RSS)
h_all$RSS <- gsub("K","", h_all$RSS)
h_all$RSS <- as.numeric(h_all$RSS)
h_all$index <- as.numeric(as.factor(h_all$Time))
h_all = h_all %>% group_by(index) %>% summarise(x = sum(RSS))
h_all$model <-"hdf5_40cores"
h_all$cell <-1300

hdf5_40 <- rbind(h_10,h_100,h_200,h_300,h_600,h_all) %>% group_by(cell,model) %>% summarise(x = max(x))


### hdf5 10 cores


mem_h_10<- read.table("mem_hdf5_10_10.txt", header=TRUE, quote="\"")
h_10 <- mem_h_10[,c(1,9)]
h_10$RSS <- gsub("G", "000000", h_10$RSS)
h_10$RSS <- gsub("M","000", h_10$RSS)
h_10$RSS <- gsub("K","", h_10$RSS)
h_10$RSS <- as.numeric(h_10$RSS)
h_10$index <- as.numeric(as.factor(h_10$Time))
h_10 = h_10 %>% group_by(index) %>% summarise(x = sum(RSS))
h_10$model <-"hdf5_10cores"
h_10$cell <-10

mem_h_100<- read.table("mem_hdf5_100_10.txt", header=TRUE, quote="\"")
h_100 <- mem_h_100[,c(1,9)]
h_100$RSS <- gsub("G", "000000", h_100$RSS)
h_100$RSS <- gsub("M","000", h_100$RSS)
h_100$RSS <- gsub("K","", h_100$RSS)
h_100$RSS <- as.numeric(h_100$RSS)
h_100$index <- as.numeric(as.factor(h_100$Time))
h_100 = h_100 %>% group_by(index) %>% summarise(x = sum(RSS))
h_100$model <-"hdf5_10cores"
h_100$cell <-100

mem_h_200<- read.table("mem_hdf5_200_10.txt", header=TRUE, quote="\"")
h_200 <- mem_h_200[,c(1,9)]
h_200$RSS <- gsub("G", "000000", h_200$RSS)
h_200$RSS <- gsub("M","000", h_200$RSS)
h_200$RSS <- gsub("K","", h_200$RSS)
h_200$RSS <- as.numeric(h_200$RSS)
h_200$index <- as.numeric(as.factor(h_200$Time))
h_200 = h_200 %>% group_by(index) %>% summarise(x = sum(RSS))
h_200$model <-"hdf5_10cores"
h_200$cell <-200

mem_h_300<- read.table("mem_hdf5_300_10.txt", header=TRUE, quote="\"")
h_300 <- mem_h_300[,c(1,9)]
h_300$RSS <- gsub("G", "000000", h_300$RSS)
h_300$RSS <- gsub("M","000", h_300$RSS)
h_300$RSS <- gsub("K","", h_300$RSS)
h_300$RSS <- as.numeric(h_300$RSS)
h_300$index <- as.numeric(as.factor(h_300$Time))
h_300 = h_300 %>% group_by(index) %>% summarise(x = sum(RSS))
h_300$model <-"hdf5_10cores"
h_300$cell <-300

mem_h_600<- read.table("mem_hdf5_600_10.txt", header=TRUE, quote="\"")
h_600 <- mem_h_600[,c(1,9)]
h_600$RSS <- gsub("G", "000000", h_600$RSS)
h_600$RSS <- gsub("M","000", h_600$RSS)
h_600$RSS <- gsub("K","", h_600$RSS)
h_600$RSS <- as.numeric(h_600$RSS)
h_600$index <- as.numeric(as.factor(h_600$Time))
h_600 = h_600 %>% group_by(index) %>% summarise(x = sum(RSS))
h_600$model <-"hdf5_10cores"
h_600$cell <-600

mem_h_all<- read.table("mem_hdf5_1300_10.txt", header=TRUE, quote="\"")
h_all <- mem_h_all[,c(1,9)]
h_all$RSS <- gsub("G", "000000", h_all$RSS)
h_all$RSS <- gsub("M","000", h_all$RSS)
h_all$RSS <- gsub("K","", h_all$RSS)
h_all$RSS <- as.numeric(h_all$RSS)
h_all$index <- as.numeric(as.factor(h_all$Time))
h_all = h_all %>% group_by(index) %>% summarise(x = sum(RSS))
h_all$model <-"hdf5_10cores"
h_all$cell <-1300

hdf5_10 <- rbind(h_10,h_100,h_200,h_300,h_600,h_all) %>% group_by(cell,model) %>% summarise(x = max(x))

### matrix

mem_h_10<- read.table("mem_matrix_10.txt", header=TRUE, quote="\"")
h_10 <- mem_h_10[,c(1,9)]
h_10$RSS <- gsub("G", "000000", h_10$RSS)
h_10$RSS <- gsub("M","000", h_10$RSS)
h_10$RSS <- gsub("K","", h_10$RSS)
h_10$RSS <- as.numeric(h_10$RSS)
h_10$index <- as.numeric(as.factor(h_10$Time))
h_10 = h_10 %>% group_by(index) %>% summarise(x = sum(RSS))
h_10$model <-"matrix_10cores"
h_10$cell <-10

mem_h_100<- read.table("mem_matrix_100.txt", header=TRUE, quote="\"")
h_100 <- mem_h_100[,c(1,9)]
h_100$RSS <- gsub("G", "000000", h_100$RSS)
h_100$RSS <- gsub("M","000", h_100$RSS)
h_100$RSS <- gsub("K","", h_100$RSS)
h_100$RSS <- as.numeric(h_100$RSS)
h_100$index <- as.numeric(as.factor(h_100$Time))
h_100 = h_100 %>% group_by(index) %>% summarise(x = sum(RSS))
h_100$model <-"matrix_10cores"
h_100$cell <-100

mem_h_200<- read.table("mem_matrix_200.txt", header=TRUE, quote="\"")
h_200 <- mem_h_200[,c(1,9)]
h_200$RSS <- gsub("G", "000000", h_200$RSS)
h_200$RSS <- gsub("M","000", h_200$RSS)
h_200$RSS <- gsub("K","", h_200$RSS)
h_200$RSS <- as.numeric(h_200$RSS)
h_200$index <- as.numeric(as.factor(h_200$Time))
h_200 = h_200 %>% group_by(index) %>% summarise(x = sum(RSS))
h_200$model <-"matrix_10cores"
h_200$cell <-200

mem_h_300<- read.table("mem_matrix_300.txt", header=TRUE, quote="\"")
h_300 <- mem_h_300[,c(1,9)]
h_300$RSS <- gsub("G", "000000", h_300$RSS)
h_300$RSS <- gsub("M","000", h_300$RSS)
h_300$RSS <- gsub("K","", h_300$RSS)
h_300$RSS <- as.numeric(h_300$RSS)
h_300$index <- as.numeric(as.factor(h_300$Time))
h_300 = h_300 %>% group_by(index) %>% summarise(x = sum(RSS))
h_300$model <-"matrix_10cores"
h_300$cell <-300

mem_h_600<- read.table("mem_matrix_600.txt", header=TRUE, quote="\"")
h_600 <- mem_h_600[,c(1,9)]
h_600$RSS <- gsub("G", "000000", h_600$RSS)
h_600$RSS <- gsub("M","000", h_600$RSS)
h_600$RSS <- gsub("K","", h_600$RSS)
h_600$RSS <- as.numeric(h_600$RSS)
h_600$index <- as.numeric(as.factor(h_600$Time))
h_600 = h_600 %>% group_by(index) %>% summarise(x = sum(RSS))
h_600$model <-"matrix_10cores"
h_600$cell <-600

mem_h_all<- read.table("mem_matrix_all.txt", header=TRUE, quote="\"")
h_all <- mem_h_all[,c(1,9)]
h_all$RSS <- gsub("G", "000000", h_all$RSS)
h_all$RSS <- gsub("M","000", h_all$RSS)
h_all$RSS <- gsub("K","", h_all$RSS)
h_all$RSS <- as.numeric(h_all$RSS)
h_all$index <- as.numeric(as.factor(h_all$Time))
h_all = h_all %>% group_by(index) %>% summarise(x = sum(RSS))
h_all$model <-"matrix_10cores"
h_all$cell <-1300

matrix <- rbind(h_10,h_100,h_200,h_300,h_600,h_all) %>% group_by(cell,model) %>% summarise(x = max(x))


### PCA

mem_h_10<- read.table("mem_pca_10.txt", header=TRUE, quote="\"")
h_10 <- mem_h_10[,c(1,9)]
h_10$RSS <- gsub("G", "000000", h_10$RSS)
h_10$RSS <- gsub("M","000", h_10$RSS)
h_10$RSS <- gsub("K","", h_10$RSS)
h_10$RSS <- as.numeric(h_10$RSS)
h_10$index <- as.numeric(as.factor(h_10$Time))
h_10 = h_10 %>% group_by(index) %>% summarise(x = sum(RSS))
h_10$model <-"pca_10cores"
h_10$cell <-10

mem_h_100<- read.table("mem_pca_100.txt", header=TRUE, quote="\"")
h_100 <- mem_h_100[,c(1,9)]
h_100$RSS <- gsub("G", "000000", h_100$RSS)
h_100$RSS <- gsub("M","000", h_100$RSS)
h_100$RSS <- gsub("K","", h_100$RSS)
h_100$RSS <- as.numeric(h_100$RSS)
h_100$index <- as.numeric(as.factor(h_100$Time))
h_100 = h_100 %>% group_by(index) %>% summarise(x = sum(RSS))
h_100$model <-"pca_10cores"
h_100$cell <-100

mem_h_200<- read.table("mem_pca_200.txt", header=TRUE, quote="\"")
h_200 <- mem_h_200[,c(1,9)]
h_200$RSS <- gsub("G", "000000", h_200$RSS)
h_200$RSS <- gsub("M","000", h_200$RSS)
h_200$RSS <- gsub("K","", h_200$RSS)
h_200$RSS <- as.numeric(h_200$RSS)
h_200$index <- as.numeric(as.factor(h_200$Time))
h_200 = h_200 %>% group_by(index) %>% summarise(x = sum(RSS))
h_200$model <-"pca_10cores"
h_200$cell <-200

mem_h_300<- read.table("mem_pca_300.txt", header=TRUE, quote="\"")
h_300 <- mem_h_300[,c(1,9)]
h_300$RSS <- gsub("G", "000000", h_300$RSS)
h_300$RSS <- gsub("M","000", h_300$RSS)
h_300$RSS <- gsub("K","", h_300$RSS)
h_300$RSS <- as.numeric(h_300$RSS)
h_300$index <- as.numeric(as.factor(h_300$Time))
h_300 = h_300 %>% group_by(index) %>% summarise(x = sum(RSS))
h_300$model <-"pca_10cores"
h_300$cell <-300

mem_h_600<- read.table("mem_pca_600.txt", header=TRUE, quote="\"")
h_600 <- mem_h_600[,c(1,9)]
h_600$RSS <- gsub("G", "000000", h_600$RSS)
h_600$RSS <- gsub("M","000", h_600$RSS)
h_600$RSS <- gsub("K","", h_600$RSS)
h_600$RSS <- as.numeric(h_600$RSS)
h_600$index <- as.numeric(as.factor(h_600$Time))
h_600 = h_600 %>% group_by(index) %>% summarise(x = sum(RSS))
h_600$model <-"pca_10cores"
h_600$cell <-600

mem_h_all<- read.table("mem_pca_all.txt", header=TRUE, quote="\"")
h_all <- mem_h_all[,c(1,9)]
h_all$RSS <- gsub("G", "000000", h_all$RSS)
h_all$RSS <- gsub("M","000", h_all$RSS)
h_all$RSS <- gsub("K","", h_all$RSS)
h_all$RSS <- as.numeric(h_all$RSS)
h_all$index <- as.numeric(as.factor(h_all$Time))
h_all = h_all %>% group_by(index) %>% summarise(x = sum(RSS))
h_all$model <-"pca_10cores"
h_all$cell <-1300
pca <- rbind(h_10,h_100,h_200,h_300,h_600,h_all) %>% group_by(cell,model) %>% summarise(x = max(x))


mem = rbind(hdf5_40,hdf5_10,matrix,pca)

memory = ggplot(mem, aes(y =x/1e6, x = cell,  col=model))+geom_point()+geom_line()+ylab("Memory(in GB)")+theme(axis.text.x=element_blank(),axis.title.x=element_blank()
)+ scale_color_manual(values=c("#000000", "#ffcc00", "#ff6633","#99ff00"))
memory

Plot2 <-ggarrange(time,memory, heights =c(1,1.2),
                  common.legend = TRUE, nrow=2, font.label= list(size=8))



