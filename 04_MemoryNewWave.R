library(NewWave)
library(SingleCellExperiment)
library(here)

#310515
load(here("BICCN_hvg.Rdata"))

n_cell <- 100000
set.seed(1234)
dati <-all_data[,sample(1:ncol(all_data),n_cell)]


res_zinbwave <- newWave(dati, n_gene_par = 100, n_cell_par = 10000,
                         X = "~batch", K=10, children= 10, commondispersion = F)

