library(scran)
library(scater)
library(batchelor)
library(TENxBrainData)

dati <- TENxBrainData()

dati <- computeLibraryFactors(dati)
dati <- logNormCounts(dati)
dec <- modelGeneVar(dati)
hvg <- getTopHVGs(dec, n=1000)

names = dati@rowRanges@elementMetadata$Symbol
hvg_names = names[hvg]
hvg = names %in% hvg_names
names(hvg) = names
name <- paste0("10x_hvg.Rdata")

save(hvg, file = name)