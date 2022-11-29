library(Seurat)
library(dplyr)
library(patchwork)

HET.data <-Read10X(data.dir = "count_IRX35_P6Irx3-5CKO-het/outs/filtered_feature_bc_matrix")
HET<-CreateSeuratObject(counts = HET.data, project = "IRX35",min.cells=1,min.features=1)



HET<-NormalizeData(HET)
HET<-ScaleData(HET)
hist(HET$nFeature_RNA,breaks=1e2,xlab="Num of genes per cell",
		main="distributions of num-of-genes per cell in P6 Het")
abline(v=1600,col=2,lwd=3)
plot(HET$nFeature_RNA,GetAssayData(object = HET, slot = 'counts')["Ctcf",])
Ctcf<-GetAssayData(object = HET, slot = 'counts')["Ctcf",]
lm(Ctcf~HET$nFeature_RNA)
#############################################

WPRE.data <-Read10X(data.dir = "count_IRX35_WPRE4_CRE_GFP_TDT_ERCC_LACZ_IRES_ref310_P6Irx3-5CKO-het/outs/raw_feature_bc_matrix")
WPRE<-CreateSeuratObject(counts = WPRE.data, project = "IRX35",min.cells=1,min.features=1)
#WPRE<-RenameCells(WPRE, add.cell.id = "wpre")
Ctcf<-GetAssayData(object = WPRE, slot = 'counts')["WPRE",]

WPRE_CELLS<-colnames(WPRE)
table(WPRE_CELLS%in%colnames(HET))
indHW<-match(colnames(HET),colnames(WPRE))
for(EXO in rownames(WPRE)){
	tmpx<-GetAssayData(object = WPRE, slot = 'counts')[EXO,][indHW]
	tmpx[is.na(tmpx)]<-0
	HET[[EXO]]<-as.numeric(tmpx)
}

#############################################
HET2<-subset(HET, subset = nFeature_RNA>1600)
hist(HET2$nFeature_RNA,breaks=1e2,xlab="Num of genes per cell",
		main="distributions of num-of-genes per cell in Het Col1 or Col10+")
Col10a1<-GetAssayData(object = HET2, slot = 'counts')["Col10a1",]
Col1a1<-GetAssayData(object = HET2, slot = 'counts')["Col1a1",]
Mmp13<-GetAssayData(object = HET2, slot = 'counts')["Mmp13",]
Ibsp<-GetAssayData(object = HET2, slot = 'counts')["Ibsp",]

plot(Col1a1+1,Mmp13+1,log="xy")
plot(Col1a1+1,HET2[["WPRE"]][,1]+1)
plot(Ibsp+1,HET2[["WPRE"]][,1]+1,log="xy")
plot(Ibsp+1,Col1a1+1,log="xy",pch=16,col=as.integer(HET2[["TDT"]][,1]>0)+1)


#############################################
HET2<-NormalizeData(HET2)
HET2<-ScaleData(HET2)
HET2<-FindVariableFeatures(HET2)

HET2<-RunPCA(HET2)
HET2<-FindNeighbors(HET2, dims = 1:10)
HET2<-FindClusters(HET2, resolution = 0.25)
HET2<-RunUMAP(HET2, reduction= "pca", dims= 1:40)
HET2<-RunTSNE(HET2, reduction= "pca", dims.use = 1:40, do.fast = T)

TSNEPlot(HET2, pt.size = 2,label=T,label.size=12)
FeaturePlot(HET2,label=T, features=c("Col1a1","Mmp13","Ibsp","Ihh",
		"Col9a1","Col10a1","Ptprc","Gypc","Tek","Cd34","Pecam1",
		"Tagln","Acta2","Twist1","Thy1","Pdgfra","Nt5e","Eng",
		"Ucma","Epyc","Tnmd","Prg4","WPRE","TDT"),
		cols=c("lightgrey","brown"),reduction="tsne",ncol=6,pt.size=1/2)
FeaturePlot(HET2,label=F, features=c("Tek","Cd34","Pecam1",
		"Tagln","Acta2","Twist1","Thy1","Pdgfra"),
		cols=c("lightgrey","brown"),reduction="tsne",pt.size=1/2)
FeaturePlot(HET2,label=T, features=c("Thy1"),
		cols=c("lightgrey","brown"),reduction="tsne",pt.size=2)
#############################################
HET.markers <- FindAllMarkers(HET2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

LSIG<-split(as.character(HET.markers$gene),HET.markers$cluster)

LSIG2<-sapply(levels(HET.markers$cluster),function(i){
	thisDEG<-HET.markers[HET.markers$cluster==i,]
	thisDEG<-thisDEG[(thisDEG[,3]-thisDEG[,4])>0.25& thisDEG[,2]>0.5 ,] #
	thisDEG<-thisDEG[order(thisDEG[,3]-thisDEG[,4],decreasing=T),]
	as.character(thisDEG$gene)
})

TFCD10<-unlist(sapply(LSIG2,function(x)head(intersect(x[!grepl("Gm|Rik|^a$",x)],c(surface.mouse,TF.mouse)),n=10)))
TOP10<-unlist(sapply(LSIG2,function(x)head(x[!grepl("Gm|Rik|^a$",x)],n=10)))
TOPCOL10<-unlist(sapply(LSIG2,function(x)head(x[grepl("^Col[0-9]+",x)&!grepl("Gm|Rik|^a$",x)],n=10)))


pdf("heatmap.top10TFCD.pdf",height=16,width=16)
	DoHeatmap(HET2,features = TOP10)
	DoHeatmap(HET2,features = TFCD10)
	DoHeatmap(HET2,features = TOPCOL10)

dev.off()

#############################################
HET2[["HET2"]]<-HET2[["seurat_clusters"]]
HET3<-subset(HET2,idents=c(0,1,9,11))
HET3<-NormalizeData(HET3)
HET3<-ScaleData(HET3)
HET3<-FindVariableFeatures(HET3)
HET3<-RunPCA(HET3)
HET3<-RunUMAP(HET3, reduction= "pca", dims= 1:40)
HET3<-RunTSNE(HET3, reduction= "pca", dims.use = 1:40, do.fast = T)
HET3<-FindNeighbors(HET3, dims = 1:10)
HET3<-FindClusters(HET3, resolution = 0.25)
table(HET3@meta.data[,"HET2"])

pdf("HET3.tsne.umap.pdf",width=32,height=16)
TSNEPlot(HET3,label=T,label.size=8, pt.size = 2,group.by="HET2")+
	UMAPPlot(HET3,label=T,label.size=8, pt.size = 2,group.by="HET2")+
	TSNEPlot(HET3,label=T,label.size=8, pt.size = 2)+
	UMAPPlot(HET3,label=T,label.size=8, pt.size = 2)
dev.off()
#############################################
HET.markers3 <- FindAllMarkers(HET3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

pdf("P10.TSNE.wo.HSC.LYMPH.pdf")
	DimPlot(object =HET3,reduction="tsne",label=T,label.size=12, pt.size = 2)
	DimPlot(object =HET3,reduction="umap",label=T,label.size=12, pt.size = 2)
dev.off()

#HET.markers3 <- FindAllMarkers(HET3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
LSIG3<-split(as.character(HET.markers3$gene),HET.markers3$cluster)

LSIG32<-sapply(levels(HET.markers3$cluster),function(i){
	thisDEG<-HET.markers3[HET.markers3$cluster==i,]
	thisDEG<-thisDEG[ abs(thisDEG[,"pct.1"]-thisDEG[,"pct.2"])>0.25& thisDEG[,2]>0.5 ,] #
	thisDEG<-thisDEG[order( abs(thisDEG[,"pct.1"]-thisDEG[,"pct.2"]),decreasing=T),]
	#thisDEG<-thisDEG[(thisDEG[,3]-thisDEG[,4])>0.25& thisDEG[,2]>0.5 ,] #
	#thisDEG<-thisDEG[order(thisDEG[,3]-thisDEG[,4],decreasing=T),]
	as.character(thisDEG$gene)
})

TFCD103<-unlist(sapply(LSIG32,function(x)head(intersect(x[!grepl("Gm|Rik|^a$",x)],c(surface.mouse,TF.mouse)),n=10)))
TOP103<-unlist(sapply(LSIG32,function(x)head(x[!grepl("Gm|Rik|^a$",x)],n=10)))
TOPCOL103<-unlist(sapply(LSIG32,function(x)head(x[grepl("^Col[0-9]+",x)&!grepl("Gm|Rik|^a$",x)],n=10)))


pdf("heatmap.C0.1.9.11.pdf",height=10,width=16)
	DoHeatmap(HET3,features = TFCD103)
	DoHeatmap(HET3,features = TOP103)
	DoHeatmap(HET3,features = TOPCOL103)

dev.off()


FeaturePlot(HET3,label=T, features=c("Col1a1","Mmp13","Ibsp","Ihh",
		"Col9a1","Col10a1","Ptprc","Gypc","Tek","Cd34","Pecam1",
		"Tagln","Acta2","Twist1","Thy1","Pdgfra","Nt5e","Eng",
		"Ucma","Epyc","Tnmd","Prg4","WPRE","TDT"),
		cols=c("lightgrey","brown"),reduction="tsne",ncol=6,pt.size=1/2)
FeaturePlot(HET3,label=T, features=c("Plin2","Aoc3","Fabp4",
		"Pparg","Col3a1","Postn","Bglap","Bglap2"),
		cols=c("lightgrey","brown"),reduction="tsne",ncol=3,pt.size=1)
FeaturePlot(HET3,label=T, features=c("Top2a","Mki67","Birc5"),
		cols=c("lightgrey","brown"),reduction="tsne",ncol=3,pt.size=1)
VlnPlot(HET3,features=c("CRE","EGFP","WPRE","TDT"),
		pt.size=0,ncol=2)
##############################################
HET3[["HET3"]]<-HET3[["seurat_clusters"]]
HET4<-subset(HET3,idents=c(0,2,4,6,7,1,3))
HET4<-NormalizeData(HET4)
HET4<-ScaleData(HET4)
HET4<-FindVariableFeatures(HET4)
HET4<-RunPCA(HET4)
HET4<-RunUMAP(HET4, reduction= "pca", dims= 1:40)
HET4<-RunTSNE(HET4, reduction= "pca", dims.use = 1:40, do.fast = T)
HET4<-FindNeighbors(HET4, dims = 1:10)
#HET4<-FindClusters(HET4, resolution = 0.25)
table(HET4[["HET3"]])

TSNEPlot(HET4,label=T,label.size=8, pt.size = 2,group.by="HET3")+
	UMAPPlot(HET4,label=T,label.size=8, pt.size = 2,group.by="HET3")+
	TSNEPlot(HET4,label=T,label.size=8, pt.size = 2)+
	UMAPPlot(HET4,label=T,label.size=8, pt.size = 2)

table(HET4[["HET3"]][,1],HET4[["seurat_clusters"]][,1])

FeaturePlot(HET4,label=T, features=c("Top2a","Mki67","Birc5",
		"WPRE","TDT","CRE"),
		cols=c("lightgrey","brown"),reduction="tsne",ncol=3,pt.size=1,label.size=8)

FeaturePlot(HET4,label=T, features=c("Top2a","Mki67","Birc5",
		"WPRE","TDT","CRE"),
		cols=c("lightgrey","brown"),slot="counts",reduction="tsne",ncol=3,pt.size=1,label.size=8)

FeaturePlot(HET4,label=T, features=c("Col9a1","Col2a1","Sox9","Ucma","Epyc",
		"Cnmd","Col10a1","Ihh"),
		cols=c("lightgrey","brown"),slot="counts",reduction="tsne",ncol=4,pt.size=1,label.size=8)


FeaturePlot(HET4,label=T, features=c("Col1a1","Mmp13","Ibsp","Runx2","Sp7","Spp1","Smpd3","Alpl",
		"WPRE","TDT","CRE","nCount_RNA","nFeature_RNA"),
		cols=c("lightgrey","brown"),reduction="tsne",ncol=5,pt.size=1,label.size=8)


FeaturePlot(HET4,label=T, features=c("Col1a1","Mmp13","Ibsp","Ihh",
		"Col9a1","Col10a1","Ptprc","Gypc","Tek","Cd34","Pecam1",
		"Pthlr","Acta2","Twist1","Thy1","Pdgfra","Nt5e","Eng",
		"Ucma","Epyc","Tnmd","Prg4","WPRE","TDT"),
		cols=c("lightgrey","brown"),reduction="tsne",ncol=6,pt.size=1)
###########################################
tabTDT<-table(HET4[["HET3"]][,1],HET4$TDT>0)
normtabTDT<-apply(tabTDT,1,function(x)x*100/sum(x))
bpp<-barplot(normtabTDT)
text(bpp,5,paste0(signif(normtabTDT[2,],2),"%"),srt=90,offset=0,pos=4,col="white",cex=2)
pie(tabTDT["1",])
pie(tabTDT["3",])

tabWPRE<-table(HET4[["HET3"]][,1],HET4$WPRE>2)
barplot(quantile(HET4$WPRE[HET4[["HET3"]][,1]==7],seq(100)/100),
	xlab="percentile",ylab="WPRE cutoff")
abline(h=quantile(HET4$WPRE[HET4[["HET3"]][,1]==7],0.25),col=2,lwd=2)
text(5,quantile(HET4$WPRE[HET4[["HET3"]][,1]==7],0.25)+3,
	quantile(HET4$WPRE[HET4[["HET3"]][,1]==7],0.25),col=2)

tabWPRE<-table(HET4[["HET3"]][,1],HET4$WPRE>14)
normtabWPRE<-apply(tabWPRE,1,function(x)x*100/sum(x))
bpp<-barplot(normtabWPRE)
text(bpp,5,paste0(signif(normtabWPRE[2,],2),"%"),srt=90,offset=0,pos=4,col="white",cex=2)

boxplot(HET4$WPRE~HET4[["HET3"]][,1],xlab="Clusters",ylab="WPRE")

vec<-HET4$WPRE>14
vec[vec]<-"WPRE>14.25"
vec[vec=="FALSE"]<-"WPRE<=14.25"

HET4[["WPREgt14"]]<-vec

TSNEPlot(HET4,label=T,label.size=8, pt.size = 1,group.by="WPREgt14")+
	TSNEPlot(HET4,label=T,label.size=8, pt.size = 1,group.by="HET3")

save(HET4,HET.markers4,file="HET4.RData")

##############################################
##############################################
