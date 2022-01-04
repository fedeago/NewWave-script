library(cluster)
library(NewWave)
library(mclust)
library(BiocParallel)
library(microbenchmark)
library(zinbwave)


load("/path/to/file/BICCN_hvg.Rdata")

# Variazione dei cores

set.seed(1234)
dati = all_data[,sample(1:ncol(all_data),100000)]

## Tempi computazionali

### 10 cores

time_newWave_genemini_10<- system.time(res_newWave_genewise_allminibatch <- newFit(dati, K=10,
                                                                                   X = "~batch", commondispersion = F, 
                                                                                   children = 10,
                                                                                   n_gene_par = 100, n_cell_par = ncol(dati)/10))
time_zinb10 <- system.time(res_zinbwave <- zinbFit(dati,
                                                   X = "~batch", K=10, 
                                                   BPPARAM = MulticoreParam(10), zeroinflation=T))

### 20 cores

time_newWave_genemini_20 <- system.time(res_newWave_genewise_allminibatch <- newFit(dati, K=10,
                                                                                    X = "~batch", commondispersion = F, 
                                                                                    children = 20,
                                                                                    n_gene_par = 100, n_cell_par = ncol(dati)/10))
time_zinb20 <- system.time(res_zinbwave <- zinbFit(dati,
                                                   X = "~batch", K=10, 
                                                   BPPARAM = MulticoreParam(20), zeroinflation=T))

### 30 cores

time_newWave_genemini_30 <- system.time(res_newWave_genewise_allminibatch <- newFit(dati, K=10,
                                                                                    X = "~batch", commondispersion = F, 
                                                                                    children = 30,
                                                                                    n_gene_par = 100, n_cell_par = ncol(dati)/10))
time_zinb30 <- system.time(res_zinbwave <- zinbFit(dati,
                                                   X = "~batch", K=10, 
                                                   BPPARAM = MulticoreParam(30), zeroinflation=T))

### 40 cores

time_newWave_genemini_40 <- system.time(res_newWave_genewise_allminibatch <- newFit(dati, K=10,
                                                                                    X = "~batch", commondispersion = F, 
                                                                                    children = 40,
                                                                                    n_gene_par = 100, n_cell_par = ncol(dati)/10))
time_zinb40 <- system.time(res_zinbwave <- zinbFit(dati,
                                                   X = "~batch", K=10, 
                                                   BPPARAM = MulticoreParam(40), zeroinflation=T))


time = as.data.frame(c(time_newWave_genemini_10[3],time_zinb10[3],time_newWave_genemini_20[3],time_zinb20[3],time_newWave_genemini_30[3],time_zinb30[3],time_newWave_genemini_40[3],time_zinb40[3]))

time$cores = c("10","10","20","20","30","30","40","40")
time$model = c("NewWave","ZiNB-WAVE")


colnames(time)[1]="Time"
time$model = as.factor(time$model)
time$Time = as.numeric(time$Time)
time$cores = as.factor(time$cores)

time_comp <- ggplot(data = time, aes(y=Time,x=cores,group=model, col=model))+
  geom_point()+
  geom_line()+ylab("Time(in minute)")+
  xlab("Number of cores")+labs(col='Model')+ scale_color_manual(values=c("#ff0033","#00ffff"))
time_comp


mem_h_10<- read.table("/path/to/file/NewWave_mem_10.txt", header=TRUE, quote="\"")
h_10 <- mem_h_10[,c(1,9)]
h_10$RSS <- gsub("G", "000000", h_10$RSS)
h_10$RSS <- gsub("M","000", h_10$RSS)
h_10$RSS <- gsub("K","", h_10$RSS)
h_10$RSS <- as.numeric(h_10$RSS)
h_10$index <- as.numeric(as.factor(h_10$Time))
h_10 = h_10 %>% group_by(index) %>% summarise(x = sum(RSS))
h_10$model <-"NewWave"
h_10$cores <-10

mem_h_100<- read.table("/path/to/file/NewWave_mem_20.txt", header=TRUE, quote="\"")
h_100 <- mem_h_100[,c(1,9)]
h_100$RSS <- gsub("G", "000000", h_100$RSS)
h_100$RSS <- gsub("M","000", h_100$RSS)
h_100$RSS <- gsub("K","", h_100$RSS)
h_100$RSS <- as.numeric(h_100$RSS)
h_100$index <- as.numeric(as.factor(h_100$Time))
h_100 = h_100 %>% group_by(index) %>% summarise(x = sum(RSS))
h_100$model <-"NewWave"
h_100$cores <-20

mem_h_200<- read.table("/path/to/file/NewWave_mem_30.txt", header=TRUE, quote="\"")
h_200 <- mem_h_200[,c(1,9)]
h_200$RSS <- gsub("G", "000000", h_200$RSS)
h_200$RSS <- gsub("M","000", h_200$RSS)
h_200$RSS <- gsub("K","", h_200$RSS)
h_200$RSS <- as.numeric(h_200$RSS)
h_200$index <- as.numeric(as.factor(h_200$Time))
h_200 = h_200 %>% group_by(index) %>% summarise(x = sum(RSS))
h_200$model <-"NewWave"
h_200$cores <-30

mem_h_300<- read.table("/path/to/file/NewWave_mem_10.txt", header=TRUE, quote="\"")
h_300 <- mem_h_300[,c(1,9)]
h_300$RSS <- gsub("G", "000000", h_300$RSS)
h_300$RSS <- gsub("M","000", h_300$RSS)
h_300$RSS <- gsub("K","", h_300$RSS)
h_300$RSS <- as.numeric(h_300$RSS)
h_300$index <- as.numeric(as.factor(h_300$Time))
h_300 = h_300 %>% group_by(index) %>% summarise(x = sum(RSS))
h_300$model <-"NewWave"
h_300$cores <-40

NewWave <- rbind(h_10,h_100,h_200,h_300) %>% group_by(cores,model) %>% summarise(x = max(x))


mem_h_10<- read.table("/path/to/file/zinbwave_mem_10.txt", header=TRUE, quote="\"")
h_10 <- mem_h_10[,c(1,9)]
h_10$RSS <- gsub("G", "000000", h_10$RSS)
h_10$RSS <- gsub("M","000", h_10$RSS)
h_10$RSS <- gsub("K","", h_10$RSS)
h_10$RSS <- as.numeric(h_10$RSS)
h_10$index <- as.numeric(as.factor(h_10$Time))
h_10 = h_10 %>% group_by(index) %>% summarise(x = sum(RSS))
h_10$model <-"ZinbWave"
h_10$cores <-10

mem_h_100<- read.table("/path/to/file/zinbwave_mem_20.txt", header=TRUE, quote="\"")
h_100 <- mem_h_100[,c(1,9)]
h_100$RSS <- gsub("G", "000000", h_100$RSS)
h_100$RSS <- gsub("M","000", h_100$RSS)
h_100$RSS <- gsub("K","", h_100$RSS)
h_100$RSS <- as.numeric(h_100$RSS)
h_100$index <- as.numeric(as.factor(h_100$Time))
h_100 = h_100 %>% group_by(index) %>% summarise(x = sum(RSS))
h_100$model <-"ZinbWave"
h_100$cores <-20

mem_h_200<- read.table("/path/to/file/zinbwave_mem_30.txt", header=TRUE, quote="\"")
h_200 <- mem_h_200[,c(1,9)]
h_200$RSS <- gsub("G", "000000", h_200$RSS)
h_200$RSS <- gsub("M","000", h_200$RSS)
h_200$RSS <- gsub("K","", h_200$RSS)
h_200$RSS <- as.numeric(h_200$RSS)
h_200$index <- as.numeric(as.factor(h_200$Time))
h_200 = h_200 %>% group_by(index) %>% summarise(x = sum(RSS))
h_200$model <-"ZinbWave"
h_200$cores <-30

mem_h_300<- read.table("/path/to/file/zinbwave_mem_10.txt", header=TRUE, quote="\"")
h_300 <- mem_h_300[,c(1,9)]
h_300$RSS <- gsub("G", "000000", h_300$RSS)
h_300$RSS <- gsub("M","000", h_300$RSS)
h_300$RSS <- gsub("K","", h_300$RSS)
h_300$RSS <- as.numeric(h_300$RSS)
h_300$index <- as.numeric(as.factor(h_300$Time))
h_300 = h_300 %>% group_by(index) %>% summarise(x = sum(RSS))
h_300$model <-"ZinbWave"
h_300$cores <-40

ZinbWave <- rbind(h_10,h_100,h_200,h_300) %>% group_by(cores,model) %>% summarise(x = max(x))


memory = rbind(NewWave,ZinbWave)
mem_comp <- ggplot(data = memory, aes(y=x/1e6,x=cores,group=model, col=model))+
  geom_point()+
  geom_line()+ylab("Memory(in GB)")+
  xlab("Number of cores")+labs(col='Model')+ scale_color_manual(values=c("#ff0033","#00ffff"))

Plot3 <-ggarrange(time_comp,mem_comp, common.legend = TRUE, ncol=2, font.label= list(size=8))
