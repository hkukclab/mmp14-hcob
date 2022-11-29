library(Seurat)
library(dplyr)
library(patchwork)

MMPhet.data <-Read10X(data.dir = "count_nelson10x_MMP14-Het/outs/filtered_feature_bc_matrix/")
MMPhet<-CreateSeuratObject(counts = MMPhet.data, project = "MMPhet",min.cells=1,min.features=1)
MMPhet<-NormalizeData(MMPhet)
MMPhet<-ScaleData(MMPhet,features=rownames(MMPhet))

MMPhetTDT.data <-Read10X(data.dir = "count_nelson10x_WPRE4_CRE_GFP_TDT_ERCC_LACZ_IRES_ref310_MMP14-Het/outs/raw_feature_bc_matrix/")
MMPhetTDT<-CreateSeuratObject(counts = MMPhetTDT.data, project = "MMPhet",min.cells=0,min.features=0)


indexo<-match(colnames(MMPhet),colnames(MMPhetTDT))

for(exo in rownames(MMPhetTDT)){
	thisExo<-as.integer(GetAssayData(MMPhetTDT)[exo ,indexo])
	cat(exo, table(thisExo>0),"\r\n")
	if(sum(thisExo>0)>0)
		MMPhet[[exo]]<-thisExo
}

MMPhet[["TDT_status"]]<-MMPhet@meta.data$TDT>0
MMPhet[["WPRE_status"]]<-MMPhet@meta.data$WPRE>0

############################################################
############################################################
MMPhet<-FindVariableFeatures(MMPhet)

MMPhet<-RunPCA(MMPhet)
MMPhet<-FindNeighbors(MMPhet, dims = 1:10)
MMPhet<-FindClusters(MMPhet, resolution = 0.25)
MMPhet<-RunUMAP(MMPhet, reduction= "pca", dims= 1:40)
MMPhet<-RunTSNE(MMPhet, reduction= "pca", dims.use = 1:40, do.fast = T)

TSNEPlot(MMPhet, pt.size = 2,label=T,label.size=8)+
	UMAPPlot(MMPhet, pt.size = 2,label=T,label.size=8)+
	UMAPPlot(MMPhet, group.by="TDT_status", cols=c("grey", "red"), pt.size = 1,label=F,label.size=12)

LSIG.MMPhet<-GETSIG(MMPhet)



FeaturePlot(MMPhet,label=T, features=c("TDT","WPRE", "Col1a1", "Runx2", "Sox9", "Col2a1", "Col10a1", "Ptprc", "Hba-x"),
		cols=c("lightgrey","brown"),reduction="tsne",ncol=3,pt.size=1)
MTREADS<-colSums(as.matrix(GetAssayData(MMPhet, slot="counts")[rownames(MMPhet)[grep("^Mt", rownames(MMPhom))], ]))
boxplot((MTREADS/MMPhet@meta.data$nCount_RNA)~MMPhet@meta.data$seurat_clusters, outline=F)
VlnPlot(MMPhet, features=c("nFeature_RNA", "nCount_RNA"), pt.size=0)

FeaturePlot(MMPhet,label=T, features=c("Pecam1"),
		cols=c("lightgrey","brown"),reduction="tsne",ncol=1,pt.size=1)

############################################################
############################################################
MMPhet_2<-subset(MMPhet,idents=c(9, 12, 13))
MMPhet_2<-NormalizeData(MMPhet_2)
MMPhet_2<-ScaleData(MMPhet_2,features=rownames(MMPhet))
MMPhet_2[["First"]]<-paste0("C_",MMPhet_2@meta.data$seurat_clusters)

MMPhet_2<-FindVariableFeatures(MMPhet_2)
MMPhet_2<-RunPCA(MMPhet_2)
MMPhet_2<-FindNeighbors(MMPhet_2, dims = 1:10)
MMPhet_2<-FindClusters(MMPhet_2, resolution = 0.25)
MMPhet_2<-RunUMAP(MMPhet_2, reduction= "pca", dims= 1:40)
MMPhet_2<-RunTSNE(MMPhet_2, reduction= "pca", dims.use = 1:40, do.fast = T)

LSIG.MMPhet_2<-GETSIG(MMPhet_2)


TSNEPlot(MMPhet_2, pt.size = 2,label=T,label.size=12)+
	UMAPPlot(MMPhet_2, pt.size = 2,label=T,label.size=12)

	UMAPPlot(MMPhet_2, group.by="TDT_status", cols=c("grey", "red"),  pt.size = 2,label=F,label.size=12)

FeaturePlot(MMPhet_2,label=T, features=c("Pecam1","WPRE", "Lepr", 
	"Col1a1", "Runx2", "Sox9", "Col2a1", "Col10a1", "Ihh", "Pthlh", 
	"Bglap", "Ifitm5", "Spp1", "Sost", "Phex", "Dmp1", "Prg4", "Cilp", "Aspn"),
		cols=c("lightgrey","brown"),reduction="umap",ncol=5,pt.size=1)


FeaturePlot(MMPhet_2,label=T, features=c("Myh11","Tagln", "Lepr", "Col1a1", "Fgfr3", "Col2a1", "Sox9", "Grem1"),
		cols=c("lightgrey","brown"),reduction="umap",ncol=3,pt.size=1)


FeaturePlot(MMPhet_2,label=T, features=c("TDT","WPRE", "Lepr", 
	"Col1a1", "Runx2", "Sox9", "Col2a1", "Col10a1", "Ihh", "Pthlh", 
	"Bglap", "Ifitm5", "Spp1", "Sost", "Phex", "Dmp1"),
		cols=c("lightgrey","brown"),reduction="umap",ncol=3,pt.size=1)


############################################################
############################################################
library(scales)

MMPhet_3<-subset(MMPhet_2,idents=c(0, 2, 3, 6, 8))
MMPhet_3<-NormalizeData(MMPhet_3)
MMPhet_3<-ScaleData(MMPhet_3,features=rownames(MMPhet))
MMPhet_3[["First"]]<-paste0("C_",MMPhet_3@meta.data$seurat_clusters)

MMPhet_3<-FindVariableFeatures(MMPhet_3)
MMPhet_3<-RunPCA(MMPhet_3)
MMPhet_3<-FindNeighbors(MMPhet_3, dims = 1:10)
MMPhet_3<-FindClusters(MMPhet_3, resolution = 0.25)
MMPhet_3<-RunUMAP(MMPhet_3, reduction= "pca", dims= 1:40)
MMPhet_3<-RunTSNE(MMPhet_3, reduction= "pca", dims.use = 1:40, do.fast = T)

tabTdt3<-t(table(MMPhet_3@meta.data$seurat_clusters, MMPhet_3@meta.data$TDT_status))[c(2,1),]

bpp<-barplot(tabTdt3)
text(bpp, colSums(tabTdt3)+2, percent(tabTdt3[1, ]/colSums(tabTdt3)), srt=0, offset=0, pos=3, cex=1, xpd=T)

LSIG.MMPhet_3<-GETSIG(MMPhet_3)


TSNEPlot(MMPhet_3, pt.size = 2,label=T,label.size=12)+
	UMAPPlot(MMPhet_3, pt.size = 2,label=T,label.size=12)

	UMAPPlot(MMPhet_3, group.by="TDT_status", cols=c("grey", "red"),  pt.size = 2,label=F,label.size=12)

FeaturePlot(MMPhet_3,label=T, features=c("Pecam1","WPRE", "Lepr", 
	"Col1a1", "Runx2", "Sox9", "Col2a1", "Col10a1", "Ihh", "Pthlh", 
	"Bglap", "Ifitm5", "Spp1", "Sost", "Phex", "Dmp1", "Prg4", "Cilp", "Aspn"),
		cols=c("lightgrey","brown"),reduction="umap",ncol=5,pt.size=1)


FeaturePlot(MMPhet_3,label=T, features=c("Myh11","Tagln", "Lepr", "Col1a1", "Fgfr3", "Col2a1", "Sox9", "Grem1"),
		cols=c("lightgrey","brown"),reduction="umap",ncol=3,pt.size=1)


FeaturePlot(MMPhet_3,label=T, features=c( "Sox9", "Col9a1", "Col2a1", "Col10a1", "Ihh", "Fgfr3", "Pthlh", 
	"Cxcl12", "Lpl", "Lepr", "Grem1", "Steap4", 
	"Col1a1", "Runx2",
	"Bglap", "Ifitm5", "Spp1", "Sost", "Phex", "Dmp1"),
		cols=c("lightgrey","brown"),reduction="umap",ncol=5,pt.size=1)

save(MMPhet, MMPhet_2, MMPhet_3, LSIG.MMPhet,
	LSIG.MMPhet_2, LSIG.MMPhet_3, file="MMPhet.RData")







