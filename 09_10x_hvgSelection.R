library(scran)
library(scater)
library(batchelor)
library(TENxBrainData)

dati <- TENxBrainData()

dati <- computeLibraryFactors(dati)
dati <- logNormCounts(dati)

hvg <- getTopHVGs(dati, n=1000)

name <- paste0("10x_hvg.Rdata")

save(hvg, file = name)