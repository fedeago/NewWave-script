library(ggplot2)
library(RColorBrewer)
library(wesanderson)

getPalette = colorRampPalette(brewer.pal(8, "Set2"))
load("/mnt/federico/BICCN_data/modelloPerPlotPCA.RData")
pca_res$Batch = as.factor(gsub("_AIBS","",x = as.character(pca_res$batch)))
load("/mnt/federico/BICCN_data/modelloPerPlot.RData")
new_res$Batch = as.factor(gsub("_AIBS","",x = as.character(new_res$batch)))
new_res$Cell = as.factor(new_res$etichette)
pca_res$Cell = as.factor(pca_res$etichette)

pNW_b = ggplot(new_res,aes(x=X1,y=X2,col=Batch))+geom_point()+theme(legend.title = element_text("Batch"))
leg <- get_legend(pNW_b)
pNW_b = ggplot(new_res,aes(x=X1,y=X2,col=Batch))+geom_point()+theme(legend.position = "none",legend.key.size = unit(1.1, "cm"),axis.title = element_blank())
ggsave(filename="Scrivania/umapNW.png",ggplot(new_res,aes(x=X1,y=X2,col=batch))+geom_point() + scale_fill_discrete(name = "Dose"))

pNW_e = ggplot(new_res,aes(x=X1,y=X2,col=Cell))+geom_point()+theme(legend.title = element_text("Cell type"))+scale_color_manual(values=getPalette(20))
leg1 <- get_legend(pNW_e)
pNW_e = ggplot(new_res,aes(x=X1,y=X2,col=Cell))+geom_point()+theme(legend.position = "none",legend.key.size = unit(1.1, "cm"),axis.title = element_blank())+scale_color_manual(values=getPalette(20))
ggsave(filename="Scrivania/umapNWEtichette.png",ggplot(new_res,aes(x=X1,y=X2,col=etichette))+geom_point()+scale_color_manual(values=getPalette(20)))



pPCA_b = ggplot(pca_res,aes(x=X1,y=X2,col=Batch))+geom_point()+theme(legend.position = "none",axis.title = element_blank())
ggsave(filename="Scrivania/umapPCA.png",ggplot(pca_res,aes(x=X1,y=X2,col=batch))+geom_point())
pPCA_e = ggplot(pca_res,aes(x=X1,y=X2,col=Cell))+geom_point()+theme(legend.position = "none",axis.title = element_blank())+scale_color_manual(values=getPalette(20))
ggsave(filename="Scrivania/umapPCAEtichette.png",ggplot(pca_res,aes(x=X1,y=X2,col=etichette))+geom_point()+scale_color_manual(values=getPalette(20)))



fin <- ggarrange(pNW_b,pPCA_b, leg,pNW_e,pPCA_e,leg1,  labels = c("A", "B","", "C","D"),   
                 ncol = 3, # Second row with 2 plots in 2 different columns
                 nrow = 2, widths = c(2.5,2.5,1), font.label= list(size=8)) 
ggsave(filename="Scrivania/allRevision.png",fin,width = 9,  height = 10)


