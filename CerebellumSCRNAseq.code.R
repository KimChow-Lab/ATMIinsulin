library(Seurat)
library(RColorBrewer)
library(presto)
library(dplyr)
library(tibble)
library(fgsea)

setwd("Project/Cerebellum")
CB4AdultMouseOldVersion=readRDS("rawData/cb_annotated_object.RDS") #24409 genes across 611034 samples
CB4AdultMouse <- UpdateSeuratObject(object = CB4AdultMouseOldVersion) #update seurat object from v2 to v3
CB4AdultMouseMyUMAP=RunUMAP(CB4AdultMouse,dims=c(1:50))

colorPalette=colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(unique(CB4AdultMouseMyUMAP$cluster))) 
pdf("CB4AdultMouseCellTypeMyUMAP.pdf",width=8)
DimPlot(CB4AdultMouseMyUMAP,group.by="cluster",reduction="umap")&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())
dev.off()

#extract the purkinje cells
PCs=subset(CB4AdultMouseMyUMAP,cluster %in% "Purkinje") 

#define the region group
PCs$RegionGroup=ifelse(PCs$region%in% c("I","II","III","CUL","VI","VII","VIII","IX","X"),"Vermis","Hemispere")
PCs<- RunTSNE(PCs, reduction = "pca", dims = 1:50)
pdf("RegionGroup4PCs8TSNE.pdf",width=8)
DimPlot(PCs,group.by="RegionGroup",reduction="tsne",cols=c("LightGrey","Orange"),raster=TRUE,pt.size=1)&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())
dev.off()


Vermis=subset(PCs,RegionGroup %in% "Vermis")
Vermis$regionGroup=ifelse(Vermis$region %in% c("I","II","III","CUL"),"Anterior","Posterior")
Vermis<- RunTSNE(Vermis, reduction = "pca", dims = 1:50)
pdf("RegionGroup4Vermis8TSNE.pdf",width=8)
DimPlot(Vermis,group.by="regionGroup",reduction="tsne",cols=c("Khaki","MediumTurquoise"),raster=TRUE,pt.size=1)&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())
dev.off()



Metabolism=read.table("MetabolismPathway4MouseGene.txt",header=T,sep="\t")
fgsea_sets<- split(Metabolism$Symbol,Metabolism$PathwayName)
VermisRegionDEG <- wilcoxauc(Vermis, 'regionGroup')
table(VermisRegionDEG$group)
regionGroup="Posterior"
clusterCell<- VermisRegionDEG %>% dplyr::filter(group == regionGroup) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0,minSize = 10)
fwrite(fgseaRes, file="VermisRegionKEGG4Graph8GSEA.txt", sep="\t", sep2=c("", " ", ""))


data=read.table("D:/Cerebellum/RegionSpecific8Kozareva/Purkinje/GSEA/Metabolism_PostvsAnter_PCsVermis.txt",header=T,sep="\t")
data=data[order(data$pval),]
data=data[1:10,]
PathwayList=factor(data$pathway,levels=rev(unique(data$pathway)))
ggplot(data, aes(x=-log10(pval), y=PathwayList,color=NES,size=size)) +
  geom_point() + scale_color_gradient2(low="navy",mid="white",high = "red")+
  theme_bw()
