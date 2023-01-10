library(Seurat)
library(ggplot2)
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

#expression of the aldoc between anterior and posterior
pdf("AldocExprVlnPlot.pdf",width=3,height=3)
VlnPlot(object = Vermis, features = c("Aldoc"),group.by="regionGroup")&theme(axis.title.x = element_blank(),axis.title.y = element_blank())
dev.off()

#aldoc distribution between anterior and posterior
Vermis$AdlocGroup=ifelse(Vermis$subcluster %in% c("Purkinje_Anti_Aldoc_1","Purkinje_Anti_Aldoc_2"),"AldocNeg","AldocPos")
tmp=paste0(Vermis$regionGroup,"_",Vermis$AdlocGroup,sep="")
tmp=data.frame(table(tmp))
colnames(tmp)=c("Group","GroupCount")
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(tmp$Group),'_')))
resultTmp=data.frame(Region=TmpInfo[,1],Aldoc=TmpInfo[,2],"Number"=tmp$GroupCount)
RegionCount=data.frame(table(Vermis$regionGroup))
colnames(RegionCount)=c("Region","RegionCount")
Result=merge(resultTmp,RegionCount,by="Region")
g=ggplot(Result, aes(Region, Number, fill=Aldoc)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(expand=c(0,0))+theme_bw()+
  scale_fill_manual(values=c("SlateBlue","Salmon"))+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("AldocDisInRegion.pdf",width=6)
print(g)
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

data=read.table("Metabolism_PostvsAnter_PCsVermis.txt",header=T,sep="\t") 
data=data[order(data$pval),]
data=data[1:10,]
PathwayList=factor(data$pathway,levels=rev(unique(data$pathway)))
pdf("VermisRegionKEGG4Graph8GSEATop10Terms.pdf",height=3,width=5)
ggplot(data, aes(x=-log10(pval), y=PathwayList,color=NES,size=size)) +
  geom_point() + scale_color_gradient2(low="navy",mid="white",high = "red")+
  theme_bw()
dev.off()
