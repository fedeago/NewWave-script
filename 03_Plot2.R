library(cluster)
library(NewWave)
library(mclust)
library(BiocParallel)
library(microbenchmark)
library(zinbwave)

load("/path/to/file/BICCN_hvg.Rdata")


# 10000 cell

## NewWave Commondispersion
set.seed(1234)
dati = all_data[,sample(1:ncol(all_data),10000)]
etichette = dati$label
time_newWave_comm_allmini<- system.time(res_newWave_commo_allminibatch <- newFit(dati, K=10,
                                                                                 X = "~batch", 
                                                                                 children = 10, n_gene_disp=100,
                                                                                 n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_mini_all <- newAIC(res_newWave_commo_allminibatch,counts(dati))
nBIC_mini_all <-newBIC(res_newWave_commo_allminibatch,counts(dati))
ari_mini_all <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(newW(res_newWave_commo_allminibatch),length(unique(etichette)))$cluster)))

stat_com_allmini = cbind(time_newWave_comm_allmini[3],
                                                     nAIC_mini_all,
                                                     nBIC_mini_all,
                                                     ari_mini_all)

## NewWave GenewiseDispersion

time_newWave_genemini <- system.time(res_newWave_genewise_parmini <- newFit(dati, K=10,
                                                                                     X = "~batch", commondispersion = F, 
                                                                                     children = 10,
                                                                                     n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_genemini <- newAIC(res_newWave_genewise_parmini,counts(dati))
nBIC_genemini <- newBIC(res_newWave_genewise_parmini,counts(dati))
ari_genemini <- mean(sapply(1:100,function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_genewise_parmini), length(unique(etichette)))$cluster)))

stat_gene_parmini = cbind(time_newWave_genemini[3],
                              nAIC_genemini,
                              nBIC_genemini,
                              ari_genemini)

## NewWave Standard

newwave_commondisp_time <- system.time(res_newWave_commondisp <- newFit(dati,
                                                                        K=10,
                                                                        X = "~batch",  children = 10))
nAIC <- newAIC(res_newWave_commondisp, counts(dati))
nBIC <- newBIC(res_newWave_commondisp, counts(dati))
ari <-  mean(sapply(1:100, function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_commondisp),length(unique(etichette)))$cluster)))

stat_std =cbind(newwave_commondisp_time[3],
                                             nAIC,
                                             nBIC,
                                             ari)
## ZinbWave Standard
time_zinb <- system.time(res_zinbwave <- zinbFit(dati,
                                                 X = "~batch", K=10, 
                                                 BPPARAM = MulticoreParam(10), zeroinflation=T))

zAIC <- zinbAIC(res_zinbwave,t(counts(dati)))
zBIC <- zinbBIC(res_zinbwave, t(counts(dati)))
ari<- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(res_zinbwave@W,length(unique(etichette)))$cluster)))

stat_zinbwave = cbind(time_zinb[3],
                      zAIC,
                      zBIC,
                                                  ari)

## Zinbwave without zero inflation

time_nbwave<- system.time(res_nbwave <- zinbFit(dati,X = "~batch", K=10,BPPARAM = MulticoreParam(10), zeroinflation=F))


nbAIC <- zinbAIC(res_nbwave,t(counts(dati)))
nbBIC <-zinbBIC(res_nbwave,t(counts(dati)))
ari <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(getW(res_nbwave),length(unique(etichette)))$cluster)))

stat_nb = cbind(time_nbwave[3],
                nbAIC,
                nbBIC,
                ari)
                      
stat_10 <- as.data.frame(rbind(cbind(stat_nb,"ZinbWave\n (Negative Binomial)"),
                            cbind(stat_zinbwave,"ZinbWave"),
                            cbind(stat_std, "NewWave"),
                            cbind(stat_com_allmini,"NewWave\n (common dispersion\n + minibatch)"),
                            cbind(stat_gene_parmini, "NewWave\n (genewise dispersion\n + minibatch)")))

stat_10$cells = "10"

# 30000


# NewWave Commondispersion
set.seed(1234)
dati = all_data[,sample(1:ncol(all_data),30000)]
etichette = dati$label
time_newWave_comm_allmini<- system.time(res_newWave_commo_allminibatch <- newFit(dati, K=10,
                                                                                 X = "~batch", 
                                                                                 children = 10, n_gene_disp=100,
                                                                                 n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_mini_all <- newAIC(res_newWave_commo_allminibatch,counts(dati))
nBIC_mini_all <-newBIC(res_newWave_commo_allminibatch,counts(dati))
ari_mini_all <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(newW(res_newWave_commo_allminibatch),length(unique(etichette)))$cluster)))

stat_com_allmini = cbind(time_newWave_comm_allmini[3],
                         nAIC_mini_all,
                         nBIC_mini_all,
                         ari_mini_all)

## NewWave GenewiseDispersion

time_newWave_genemini <- system.time(res_newWave_genewise_parmini <- newFit(dati, K=10,
                                                                            X = "~batch", commondispersion = F, 
                                                                            children = 10,
                                                                            n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_genemini <- newAIC(res_newWave_genewise_parmini,counts(dati))
nBIC_genemini <- newBIC(res_newWave_genewise_parmini,counts(dati))
ari_genemini <- mean(sapply(1:100,function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_genewise_parmini), length(unique(etichette)))$cluster)))

stat_gene_parmini = cbind(time_newWave_genemini[3],
                          nAIC_genemini,
                          nBIC_genemini,
                          ari_genemini)

## NewWave Standard

newwave_commondisp_time <- system.time(res_newWave_commondisp <- newFit(dati,
                                                                        K=10,
                                                                        X = "~batch",  children = 10))
nAIC <- newAIC(res_newWave_commondisp, counts(dati))
nBIC <- newBIC(res_newWave_commondisp, counts(dati))
ari <-  mean(sapply(1:100, function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_commondisp),length(unique(etichette)))$cluster)))

stat_std =cbind(newwave_commondisp_time[3],
                nAIC,
                nBIC,
                ari)
## ZinbWave Standard
time_zinb <- system.time(res_zinbwave <- zinbFit(dati,
                                                 X = "~batch", K=10, 
                                                 BPPARAM = MulticoreParam(10), zeroinflation=T))

zAIC <- zinbAIC(res_zinbwave,t(counts(dati)))
zBIC <- zinbBIC(res_zinbwave, t(counts(dati)))
ari<- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(res_zinbwave@W,length(unique(etichette)))$cluster)))

stat_zinbwave = cbind(time_zinb[3],
                      zAIC,
                      zBIC,
                      ari)

## Zinbwave without zero inflation

time_nbwave<- system.time(res_nbwave <- zinbFit(dati,X = "~batch", K=10,BPPARAM = MulticoreParam(10), zeroinflation=F))


nbAIC <- zinbAIC(res_nbwave,t(counts(dati)))
nbBIC <-zinbBIC(res_nbwave,t(counts(dati)))
ari <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(getW(res_nbwave),length(unique(etichette)))$cluster)))

stat_nb = cbind(time_nbwave[3],
                nbAIC,
                nbBIC,
                ari)

stat_30 <- as.data.frame(rbind(cbind(stat_nb,"ZinbWave\n (Negative Binomial)"),
                               cbind(stat_zinbwave,"ZinbWave"),
                               cbind(stat_std, "NewWave"),
                               cbind(stat_com_allmini,"NewWave\n (common dispersion\n + minibatch)"),
                               cbind(stat_gene_parmini, "NewWave\n (genewise dispersion\n + minibatch)")))

stat_30$cells = "30"


# 50000 cell

# NewWave Commondispersion
set.seed(1234)
dati = all_data[,sample(1:ncol(all_data),50000)]
etichette = dati$label
time_newWave_comm_allmini<- system.time(res_newWave_commo_allminibatch <- newFit(dati, K=10,
                                                                                 X = "~batch",
                                                                                 children = 10, n_gene_disp=100,
                                                                                 n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_mini_all <- newAIC(res_newWave_commo_allminibatch,counts(dati))
nBIC_mini_all <-newBIC(res_newWave_commo_allminibatch,counts(dati))
ari_mini_all <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(newW(res_newWave_commo_allminibatch),length(unique(etichette)))$cluster)))

stat_com_allmini = cbind(time_newWave_comm_allmini[3],
                         nAIC_mini_all,
                         nBIC_mini_all,
                         ari_mini_all)

## NewWave GenewiseDispersion

time_newWave_genemini <- system.time(res_newWave_genewise_parmini <- newFit(dati, K=10,
                                                                            X = "~batch", commondispersion = F,
                                                                            children = 10,
                                                                            n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_genemini <- newAIC(res_newWave_genewise_parmini,counts(dati))
nBIC_genemini <- newBIC(res_newWave_genewise_parmini,counts(dati))
ari_genemini <- mean(sapply(1:100,function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_genewise_parmini), length(unique(etichette)))$cluster)))

stat_gene_parmini = cbind(time_newWave_genemini[3],
                          nAIC_genemini,
                          nBIC_genemini,
                          ari_genemini)

## NewWave Standard

newwave_commondisp_time <- system.time(res_newWave_commondisp <- newFit(dati,
                                                                        K=10,
                                                                        X = "~batch",  children = 10))
nAIC <- newAIC(res_newWave_commondisp, counts(dati))
nBIC <- newBIC(res_newWave_commondisp, counts(dati))
ari <-  mean(sapply(1:100, function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_commondisp),length(unique(etichette)))$cluster)))

stat_std =cbind(newwave_commondisp_time[3],
                nAIC,
                nBIC,
                ari)
## ZinbWave Standard
time_zinb <- system.time(res_zinbwave <- zinbFit(dati,
                                                 X = "~batch", K=10,
                                                 BPPARAM = MulticoreParam(10), zeroinflation=T))

zAIC <- zinbAIC(res_zinbwave,t(counts(dati)))
zBIC <- zinbBIC(res_zinbwave, t(counts(dati)))
ari<- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(res_zinbwave@W,length(unique(etichette)))$cluster)))

stat_zinbwave = cbind(time_zinb[3],
                      zAIC,
                      zBIC,
                      ari)

## Zinbwave without zero inflation

time_nbwave<- system.time(res_nbwave <- zinbFit(dati,X = "~batch", K=10,BPPARAM = MulticoreParam(10), zeroinflation=F))


nbAIC <- zinbAIC(res_nbwave,t(counts(dati)))
nbBIC <-zinbBIC(res_nbwave,t(counts(dati)))
ari <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(getW(res_nbwave),length(unique(etichette)))$cluster)))

stat_nb = cbind(time_nbwave[3],
                nbAIC,
                nbBIC,
                ari)

stat_50 <- as.data.frame(rbind(cbind(stat_nb,"ZinbWave\n (Negative Binomial)"),
                               cbind(stat_zinbwave,"ZinbWave"),
                               cbind(stat_std, "NewWave"),
                               cbind(stat_com_allmini,"NewWave\n (common dispersion\n + minibatch)"),
                               cbind(stat_gene_parmini, "NewWave\n (genewise dispersion\n + minibatch)")))

stat_50$cells = "50"

# 100000 cell

# NewWave Commondispersion
set.seed(1234)
dati = all_data[,sample(1:ncol(all_data),100000)]
etichette = dati$label
time_newWave_comm_allmini<- system.time(res_newWave_commo_allminibatch <- newFit(dati, K=10,
                                                                                 X = "~batch",
                                                                                 children = 10, n_gene_disp=100,
                                                                                 n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_mini_all <- newAIC(res_newWave_commo_allminibatch,counts(dati))
nBIC_mini_all <-newBIC(res_newWave_commo_allminibatch,counts(dati))
ari_mini_all <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(newW(res_newWave_commo_allminibatch),length(unique(etichette)))$cluster)))

stat_com_allmini = cbind(time_newWave_comm_allmini[3],
                         nAIC_mini_all,
                         nBIC_mini_all,
                         ari_mini_all)

## NewWave GenewiseDispersion

time_newWave_genemini <- system.time(res_newWave_genewise_parmini <- newFit(dati, K=10,
                                                                            X = "~batch", commondispersion = F,
                                                                            children = 10,
                                                                            n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_genemini <- newAIC(res_newWave_genewise_parmini,counts(dati))
nBIC_genemini <- newBIC(res_newWave_genewise_parmini,counts(dati))
ari_genemini <- mean(sapply(1:100,function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_genewise_parmini), length(unique(etichette)))$cluster)))

stat_gene_parmini = cbind(time_newWave_genemini[3],
                          nAIC_genemini,
                          nBIC_genemini,
                          ari_genemini)

## NewWave Standard

newwave_commondisp_time <- system.time(res_newWave_commondisp <- newFit(dati,
                                                                        K=10,
                                                                        X = "~batch",  children = 10))
nAIC <- newAIC(res_newWave_commondisp, counts(dati))
nBIC <- newBIC(res_newWave_commondisp, counts(dati))
ari <-  mean(sapply(1:100, function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_commondisp),length(unique(etichette)))$cluster)))

stat_std =cbind(newwave_commondisp_time[3],
                nAIC,
                nBIC,
                ari)
## ZinbWave Standard
time_zinb <- system.time(res_zinbwave <- zinbFit(dati,
                                                 X = "~batch", K=10,
                                                 BPPARAM = MulticoreParam(10), zeroinflation=T))

zAIC <- zinbAIC(res_zinbwave,t(counts(dati)))
zBIC <- zinbBIC(res_zinbwave, t(counts(dati)))
ari<- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(res_zinbwave@W,length(unique(etichette)))$cluster)))

stat_zinbwave = cbind(time_zinb[3],
                      zAIC,
                      zBIC,
                      ari)

## Zinbwave without zero inflation

time_nbwave<- system.time(res_nbwave <- zinbFit(dati,X = "~batch", K=10,BPPARAM = MulticoreParam(10), zeroinflation=F))


nbAIC <- zinbAIC(res_nbwave,t(counts(dati)))
nbBIC <-zinbBIC(res_nbwave,t(counts(dati)))
ari <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(getW(res_nbwave),length(unique(etichette)))$cluster)))

stat_nb = cbind(time_nbwave[3],
                nbAIC,
                nbBIC,
                ari)

stat_100 <- as.data.frame(rbind(cbind(stat_nb,"ZinbWave\n (Negative Binomial)"),
                               cbind(stat_zinbwave,"ZinbWave"),
                               cbind(stat_std, "NewWave"),
                               cbind(stat_com_allmini,"NewWave\n (common dispersion\n + minibatch)"),
                               cbind(stat_gene_parmini, "NewWave\n (genewise dispersion\n + minibatch)")))

stat_100$cells = "100"

# 200000 cell

# NewWave Commondispersion
set.seed(1234)
dati = all_data[,sample(1:ncol(all_data),200000)]
etichette = dati$label
time_newWave_comm_allmini<- system.time(res_newWave_commo_allminibatch <- newFit(dati, K=10,
                                                                                 X = "~batch",
                                                                                 children = 10, n_gene_disp=100,
                                                                                 n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_mini_all <- newAIC(res_newWave_commo_allminibatch,counts(dati))
nBIC_mini_all <-newBIC(res_newWave_commo_allminibatch,counts(dati))
ari_mini_all <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(newW(res_newWave_commo_allminibatch),length(unique(etichette)))$cluster)))

stat_com_allmini = cbind(time_newWave_comm_allmini[3],
                         nAIC_mini_all,
                         nBIC_mini_all,
                         ari_mini_all)

## NewWave GenewiseDispersion

time_newWave_genemini <- system.time(res_newWave_genewise_parmini <- newFit(dati, K=10,
                                                                            X = "~batch", commondispersion = F,
                                                                            children = 10,
                                                                            n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_genemini <- newAIC(res_newWave_genewise_parmini,counts(dati))
nBIC_genemini <- newBIC(res_newWave_genewise_parmini,counts(dati))
ari_genemini <- mean(sapply(1:100,function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_genewise_parmini), length(unique(etichette)))$cluster)))

stat_gene_parmini = cbind(time_newWave_genemini[3],
                          nAIC_genemini,
                          nBIC_genemini,
                          ari_genemini)

## NewWave Standard

newwave_commondisp_time <- system.time(res_newWave_commondisp <- newFit(dati,
                                                                        K=10,
                                                                        X = "~batch",  children = 10))
nAIC <- newAIC(res_newWave_commondisp, counts(dati))
nBIC <- newBIC(res_newWave_commondisp, counts(dati))
ari <-  mean(sapply(1:100, function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_commondisp),length(unique(etichette)))$cluster)))

stat_std =cbind(newwave_commondisp_time[3],
                nAIC,
                nBIC,
                ari)
## ZinbWave Standard
time_zinb <- system.time(res_zinbwave <- zinbFit(dati,
                                                 X = "~batch", K=10,
                                                 BPPARAM = MulticoreParam(10), zeroinflation=T))

zAIC <- zinbAIC(res_zinbwave,t(counts(dati)))
zBIC <- zinbBIC(res_zinbwave, t(counts(dati)))
ari<- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(res_zinbwave@W,length(unique(etichette)))$cluster)))

stat_zinbwave = cbind(time_zinb[3],
                      zAIC,
                      zBIC,
                      ari)

## Zinbwave without zero inflation

time_nbwave<- system.time(res_nbwave <- zinbFit(dati,X = "~batch", K=10,BPPARAM = MulticoreParam(10), zeroinflation=F))


nbAIC <- zinbAIC(res_nbwave,t(counts(dati)))
nbBIC <-zinbBIC(res_nbwave,t(counts(dati)))
ari <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(getW(res_nbwave),length(unique(etichette)))$cluster)))

stat_nb = cbind(time_nbwave[3],
                nbAIC,
                nbBIC,
                ari)

stat_200 <- as.data.frame(rbind(cbind(stat_nb,"ZinbWave\n (Negative Binomial)"),
                               cbind(stat_zinbwave,"ZinbWave"),
                               cbind(stat_std, "NewWave"),
                               cbind(stat_com_allmini,"NewWave\n (common dispersion\n + minibatch)"),
                               cbind(stat_gene_parmini, "NewWave\n (genewise dispersion\n + minibatch)")))

stat_200$cells = "200"

# 3120000 cell

## NewWave Commondispersion

set.seed(1234)
dati = all_data
etichette = dati$label
time_newWave_comm_allmini<- system.time(res_newWave_commo_allminibatch <- newFit(dati, K=10,
                                                                                 X = "~batch",
                                                                                 children = 10, n_gene_disp=100,
                                                                                 n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_mini_all <- newAIC(res_newWave_commo_allminibatch,counts(dati))
nBIC_mini_all <-newBIC(res_newWave_commo_allminibatch,counts(dati))
ari_mini_all <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(newW(res_newWave_commo_allminibatch),length(unique(etichette)))$cluster)))

stat_com_allmini = cbind(time_newWave_comm_allmini[3],
                         nAIC_mini_all,
                         nBIC_mini_all,
                         ari_mini_all)

## NewWave GenewiseDispersion

time_newWave_genemini <- system.time(res_newWave_genewise_parmini <- newFit(dati, K=10,
                                                                            X = "~batch", commondispersion = F,
                                                                            children = 10,
                                                                            n_gene_par = 100, n_cell_par = ncol(dati)/10))


nAIC_genemini <- newAIC(res_newWave_genewise_parmini,counts(dati))
nBIC_genemini <- newBIC(res_newWave_genewise_parmini,counts(dati))
ari_genemini <- mean(sapply(1:100,function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_genewise_parmini), length(unique(etichette)))$cluster)))

stat_gene_parmini = cbind(time_newWave_genemini[3],
                          nAIC_genemini,
                          nBIC_genemini,
                          ari_genemini)

## NewWave Standard

newwave_commondisp_time <- system.time(res_newWave_commondisp <- newFit(dati,
                                                                        K=10,
                                                                        X = "~batch",  children = 10))
nAIC <- newAIC(res_newWave_commondisp, counts(dati))
nBIC <- newBIC(res_newWave_commondisp, counts(dati))
ari <-  mean(sapply(1:100, function(y)adjustedRandIndex(as.numeric(etichette),kmeans(newW(res_newWave_commondisp),length(unique(etichette)))$cluster)))

stat_std =cbind(newwave_commondisp_time[3],
                nAIC,
                nBIC,
                ari)

## ZinbWave Standard

time_zinb <- system.time(res_zinbwave <- zinbFit(dati,
                                                 X = "~batch", K=10,
                                                 BPPARAM = MulticoreParam(10), zeroinflation=T))

zAIC <- zinbAIC(res_zinbwave,t(counts(dati)))
zBIC <- zinbBIC(res_zinbwave, t(counts(dati)))
ari<- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(res_zinbwave@W,length(unique(etichette)))$cluster)))

stat_zinbwave = cbind(time_zinb[3],
                      zAIC,
                      zBIC,
                      ari)

## Zinbwave without zero inflation

time_nbwave<- system.time(res_nbwave <- zinbFit(dati,X = "~batch", K=10,BPPARAM = MulticoreParam(10), zeroinflation=F))


nbAIC <- zinbAIC(res_nbwave,t(counts(dati)))
nbBIC <-zinbBIC(res_nbwave,t(counts(dati)))
ari <- mean(sapply(1:100, function(y)adjustedRandIndex(etichette,kmeans(getW(res_nbwave),length(unique(etichette)))$cluster)))

stat_nb = cbind(time_nbwave[3],
                nbAIC,
                nbBIC,
                ari)

stat_300 <- as.data.frame(rbind(cbind(stat_nb,"ZinbWave\n (Negative Binomial)"),
                               cbind(stat_zinbwave,"ZinbWave"),
                               cbind(stat_std, "NewWave"),
                               cbind(stat_com_allmini,"NewWave\n (common dispersion\n + minibatch)"),
                               cbind(stat_gene_parmini, "NewWave\n (genewise dispersion\n + minibatch)")))

stat_300$cells = "300"

time_pattern <- rbind(stat_10,stat_30,stat_50,stat_100,stat_200, stat_300)


time_pattern <- data.frame(time_pattern)
time_pattern$cells <- as.factor(as.numeric(time_pattern$cells))
time_pattern$model <- as.factor(time_pattern$V5)
time_pattern$V1 = as.numeric(time_pattern$V1)
time_pattern$ari = as.numeric(time_pattern$ari)

comparison <- ggplot(time_pattern, aes(x = as.numeric(as.character(cells)), y =V1/60, group = model, col=model))+geom_point()+geom_line()+ylab("Time(in hours)")+theme(axis.text.x=element_blank(),axis.title.x=element_blank())+ scale_color_manual(values=c("#336600", "#0000ff","#00ffff","#ff0033","#9900ff"))
ari <- ggplot(time_pattern, aes(x = as.numeric(as.character(cells)), y =ari, group = model, col=model))+geom_point()+geom_line()+ylab("ARI")+xlab("number of cells (thousands)")+ scale_color_manual(values=c("#336600", "#0000ff","#00ffff","#ff0033","#9900ff"))

Plot2 <-ggarrange(comparison,ari, heights =c(1,1.2),
              common.legend = TRUE, nrow=2, font.label= list(size=8))




