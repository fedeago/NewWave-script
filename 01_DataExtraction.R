library(R.utils)
library(dplyr)
library(SingleCellExperiment)
library(Matrix)
library(here)

dir.create(here("CellV2"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v2/mouse/processed/analysis/10X_cells_v2_AIBS/QC.csv", destfile = here("CellV2", "qc.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v2/mouse/processed/analysis/10X_cells_v2_AIBS/barcode.tsv", destfile = here("CellV2/barcode.tsv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v2/mouse/processed/analysis/10X_cells_v2_AIBS/cluster.annotation.csv", destfile = here("CellV2/cluster.annotation.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v2/mouse/processed/analysis/10X_cells_v2_AIBS/cluster.membership.csv", destfile = here("CellV2/cluster.membership.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v2/mouse/processed/analysis/10X_cells_v2_AIBS/features.tsv.gz", destfile = here("CellV2/features.tsv.gz"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v2/mouse/processed/analysis/10X_cells_v2_AIBS/matrix.mtx.gz", destfile = here("CellV2/matrix.mtx.gz"))

cluster_annotation<-read.csv(here("CellV2/cluster.annotation.csv"), row.names = "cluster_id")
cluster_membership <- read.csv(here("CellV2/cluster.membership.csv"), row.names = "X")
qc <- read.csv(here("CellV2/qc.csv"), row.names = "X")
barcode <- read.table(here("CellV2/barcode.tsv"), header = TRUE)
features <- read.table(here("CellV2/features.tsv.gz"))
gunzip(here("CellV2/matrix.mtx.gz"), remove=FALSE)

matrix = readMM(here("CellV2/matrix.mtx"))

matrix = as.matrix(matrix)
barcode$X..x. = gsub('\"', '',barcode$X..x.)
barcode$X..x. = gsub(',', '',barcode$X..x.)

colnames(matrix) = barcode$X..x.
rownames(matrix) = features$V2

label = data.frame(cell = rownames(cluster_membership), label= cluster_annotation$subclass_label[cluster_membership$x])

label = label[which(label$cell %in% qc$x),]

matrix = matrix[,qc$x]

cellV2 = SingleCellExperiment(assay = list(counts = matrix))
cellV2$label = label$label
cellV2$batch = "cell_V2"

# CELLV3
dir.create(here("CellV3"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v3/mouse/processed/analysis/10X_cells_v3_AIBS/QC.csv", destfile = here("CellV3/QC.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v3/mouse/processed/analysis/10X_cells_v3_AIBS/barcode.tsv", destfile = here("CellV3/barcode.tsv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v3/mouse/processed/analysis/10X_cells_v3_AIBS/cluster.annotation.csv", destfile = here("CellV3/cluster.annotation.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v3/mouse/processed/analysis/10X_cells_v3_AIBS/cluster.membership.csv", destfile = here("CellV3/cluster.membership.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v3/mouse/processed/analysis/10X_cells_v3_AIBS/features.tsv.gz", destfile = here("CellV3/features.tsv.gz"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/scell/10x_v3/mouse/processed/analysis/10X_cells_v3_AIBS/matrix.mtx.gz", destfile = here("CellV3/matrix.mtx.gz"))

cluster_annotation<-read.csv(here("CellV3/cluster.annotation.csv"), row.names = "cluster_id")
cluster_membership <- read.csv(here("CellV3/cluster.membership.csv"), row.names = "X")
qc<- read.csv(here("CellV3/QC.csv"), row.names = "X")
barcode = read.table(here("CellV3/barcode.tsv"), header = TRUE)
features = read.table(here("CellV3/features.tsv.gz"))
gunzip(here("CellV3/matrix.mtx.gz"), remove=FALSE)

matrix = readMM(here("CellV3/matrix.mtx"))

matrix = as.matrix(matrix)
barcode$X..x. = gsub('\"', '',barcode$X..x.)
barcode$X..x. = gsub(',', '',barcode$X..x.)

colnames(matrix) = barcode$X..x.
rownames(matrix) = features$V2

label = data.frame(cell = rownames(cluster_membership), label= cluster_annotation$subclass_label[cluster_membership$x])

label = label[which(label$cell %in% qc$x),]

matrix = matrix[,qc$x]

cellV3 = SingleCellExperiment(assay = list(counts = matrix))
cellV3$label = label$label
cellV3$batch = "cell_V3"

# NUCLEUSV2
dir.create(here("NucleusV2"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v2/mouse/processed/analysis/10X_nuclei_v2_AIBS/QC.csv", destfile = here("NucleusV2/QC.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v2/mouse/processed/analysis/10X_nuclei_v2_AIBS/barcode.tsv", destfile = here("NucleusV2/barcode.tsv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v2/mouse/processed/analysis/10X_nuclei_v2_AIBS/cluster.annotation.csv", destfile = here("NucleusV2/cluster.annotation.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v2/mouse/processed/analysis/10X_nuclei_v2_AIBS/cluster.membership.csv", destfile = here("NucleusV2/cluster.membership.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v2/mouse/processed/analysis/10X_nuclei_v2_AIBS/features.tsv.gz", destfile = here("NucleusV2/features.tsv.gz"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v2/mouse/processed/analysis/10X_nuclei_v2_AIBS/matrix.mtx.gz", destfile = here("NucleusV2/matrix.mtx.gz"))


cluster_annotation <-read.csv(here("NucleusV2/cluster.annotation.csv"), row.names = "cluster_id")
cluster_membership <- read.csv(here("NucleusV2/cluster.membership.csv"), row.names = "X")
qc <- read.csv(here("NucleusV2/QC.csv"), row.names = "X")
barcode = read.table(here("NucleusV2/barcode.tsv"), header = TRUE)
features = read.table(here("NucleusV2/features.tsv.gz"))
gunzip(here("NucleusV2/matrix.mtx.gz"), remove=FALSE)

matrix = readMM(here("NucleusV2/matrix.mtx"))

matrix = as.matrix(matrix)
barcode$X..x. = gsub('\"', '',barcode$X..x.)
barcode$X..x. = gsub(',', '',barcode$X..x.)

colnames(matrix) = barcode$X..x.
rownames(matrix) = features$V2

label = data.frame(cell = rownames(cluster_membership), label= cluster_annotation$subclass_label[cluster_membership$x])

label = label[which(label$cell %in% qc$x),]
matrix = matrix[,qc$x]

nucleusV2 = SingleCellExperiment(assay = list(counts = matrix))
nucleusV2$label = label$label
nucleusV2$batch = "nucleus_V2"

## NUCLEUSV3
dir.create(here("NucleusV3"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/QC.csv", destfile = here("NucleusV3/QC.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/barcode.tsv", destfile = here("NucleusV3/barcode.tsv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/cluster.annotation.csv", destfile = here("NucleusV3/cluster.annotation.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/cluster.membership.csv", destfile = here("NucleusV3/cluster.membership.csv"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/features.tsv.gz", destfile = here("NucleusV3/features.tsv.gz"))
download.file("http://data.nemoarchive.org/biccn/lab/zeng/transcriptome/sncell/10x_v3/mouse/processed/analysis/10X_nuclei_v3_AIBS/matrix.mtx.gz", destfile = here("NucleusV3/matrix.mtx.gz"))

cluster_annotation<-read.csv(here("NucleusV3/cluster.annotation.csv"), row.names = "cluster_id")
cluster_membership <- read.csv(here("NucleusV3/cluster.membership.csv"), row.names = "X")
qc <- read.csv(here("NucleusV3/QC.csv"), row.names = "X")
barcode = read.table(here("NucleusV3/barcode.tsv"), header = TRUE)
features = read.table(here("NucleusV3/features.tsv.gz"))
gunzip(here("NucleusV3/matrix.mtx.gz"), remove=FALSE)

matrix = readMM(here("NucleusV3/matrix.mtx"))


matrix = as.matrix(matrix)
barcode$X..x. = gsub('\"', '',barcode$X..x.)
barcode$X..x. = gsub(',', '',barcode$X..x.)

colnames(matrix) = barcode$X..x.
rownames(matrix) = features$V2

label = data.frame(cell = rownames(cluster_membership), label= cluster_annotation$subclass_label[cluster_membership$x])

label = label[which(label$cell %in% qc$x),]
matrix = matrix[,qc$x]
nucleusV3 = SingleCellExperiment(assay = list(counts = matrix))
nucleusV3$label = label$label
nucleusV3$batch = "nucleus_V3"

all_data = cbind(cellV2,cellV3,nucleusV2,nucleusV3)

all_data = save(all_data, file = here("BICCN_data.Rdata"))
