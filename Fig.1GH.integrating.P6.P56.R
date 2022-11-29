library(Seurat)
library(tidyverse)
library(patchwork)
library(scales)

source("helpers.GETSIG.R")

load("P6/HET3.RData")

UMAPPlot(HET3,label=T,label.size=8, pt.size = 2)
HET3[["HET3"]]<-HET3[["seurat_clusters"]]
HET4<-subset(HET3,idents=c(0,1,2,3,5,6,8))
HET4<-NormalizeData(HET4)
HET4<-ScaleData(HET4)
HET4<-FindVariableFeatures(HET4)
HET4<-RunPCA(HET4)
HET4<-RunUMAP(HET4, reduction= "pca", dims= 1:40)
HET4<-RunTSNE(HET4, reduction= "pca", dims.use = 1:40, do.fast = T)
HET4<-FindNeighbors(HET4, dims = 1:10)
HET4<-HET4(HET4, resolution = 0.25)

HET4[["oneBased"]]<-paste0("N",formatC(as.integer(Idents(HET4)), width=2, flag="0"))

###################################
load("mmp14-data/MMPhet.RData")


P6Wk8<-RunCCA(HET4,MMPhet_3)
P6Wk8[["OldClust"]]<-P6Wk8@meta.data$oneBased

P6Wk8<-NormalizeData(P6Wk8)
P6Wk8<-ScaleData(P6Wk8)
P6Wk8<-FindVariableFeatures(P6Wk8)
P6Wk8<-RunPCA(P6Wk8)
P6Wk8<-RunUMAP(P6Wk8, reduction= "cca", dims= 1:10)
P6Wk8<-RunTSNE(P6Wk8, reduction= "cca", dims.use = 1:40, do.fast = T)
P6Wk8<-FindNeighbors(P6Wk8, reduction= "cca", dims = 1:10)
P6Wk8<-FindClusters(P6Wk8, resolution = 0.25)

P6Wk8[["ioneBased"]]<-paste0("i",formatC(as.integer(Idents(P6Wk8)), width=2, flag="0"))
Idents(P6Wk8)<-"ioneBased"

UMAPPlot(HET4,label=T,label.size=8, pt.size = 2)
UMAPPlot(P6Wk8,label=T, label.size=8, pt.size = 2, group.by="ioneBased")

table(P6Wk8@meta.data$ioneBased, P6Wk8@meta.data$orig.ident)


LSIG.HET4<-GETSIG(HET4)

#######################################
LSIG.P6Wk8<-GETSIG(P6Wk8)

table(P6Wk8@meta.data$oneBased, P6Wk8@meta.data$ioneBased)


FeaturePlot(P6Wk8, features=c("Sost"))
FeaturePlot(HET4, features=c("Gypa", "Hba-a1"))
FeaturePlot(P6Wk8, features=c("Hba-a1"))
FeaturePlot(MMPhet_3, features=c("Gypa", "Hba-a1", "Hba-a2", "nFeature_RNA"))

plot(P6Wk8@reductions$umap@cell.embeddings)
abline(h=-6.1)

P6Wk8[["nonRBCs"]]<- P6Wk8@reductions$umap@cell.embeddings[,2] > -6.1 & P6Wk8@meta.data$ioneBased != "i08"

P6Wk8.2<-subset(P6Wk8, subset=nonRBCs)
UMAPPlot(P6Wk8.2,label=T, label.size=8, pt.size = 2, group.by="ioneBased")

P6Wk8.2<-NormalizeData(P6Wk8.2)
P6Wk8.2<-ScaleData(P6Wk8.2)
P6Wk8.2<-FindVariableFeatures(P6Wk8.2)
P6Wk8.2<-RunPCA(P6Wk8.2)
P6Wk8.2<-RunUMAP(P6Wk8.2, reduction= "cca", dims= 1:10)
P6Wk8.2<-RunTSNE(P6Wk8.2, reduction= "cca", dims.use = 1:40, do.fast = T)
P6Wk8.2<-FindNeighbors(P6Wk8.2, reduction= "cca", dims = 1:10)
P6Wk8.2<-FindClusters(P6Wk8.2, resolution = 0.25)

UMAPPlot(P6Wk8.2,label=T, label.size=8, pt.size = 2)
P6Wk8.2[["ioneBased2"]]<-paste0("i",formatC(as.integer(Idents(P6Wk8.2)), width=1, flag="0"))
Idents(P6Wk8.2)<-"ioneBased2"

pdf("umap.p6.annaP56.merge.chondro-osteo.pdf", width=8)
	UMAPPlot(P6Wk8.2,label=T, label.size=8, pt.size = 2, group.by="ioneBased2")
	UMAPPlot(P6Wk8.2,label=T, label.size=8, pt.size = 2, group.by="orig.ident")
dev.off()

table(P6Wk8.2@meta.data$ioneBased2, P6Wk8.2@meta.data$orig.ident)

#####################################################
plot(P6Wk8.3@reductions$umap@cell.embeddings)
abline(v=-3)
P6Wk8.3<-subset(P6Wk8.2, idents=c("i2", "i3", "i5"))
P6Wk8.3[["nonChondro"]]<- P6Wk8.3@reductions$umap@cell.embeddings[,1] < -3
P6Wk8.3<-subset(P6Wk8.3, subset=nonChondro)

UMAPPlot(P6Wk8.3,label=T, label.size=8, pt.size = 2)

P6Wk8.3<-NormalizeData(P6Wk8.3)
P6Wk8.3<-ScaleData(P6Wk8.3)
P6Wk8.3<-FindVariableFeatures(P6Wk8.3)
P6Wk8.3<-RunPCA(P6Wk8.3)
P6Wk8.3<-RunUMAP(P6Wk8.3, reduction= "cca", dims= 1:10)
P6Wk8.3<-RunTSNE(P6Wk8.3, reduction= "cca", dims.use = 1:40, do.fast = T)
P6Wk8.3<-FindNeighbors(P6Wk8.3, reduction= "cca", dims = 1:10)
P6Wk8.3<-FindClusters(P6Wk8.3, resolution = 0.25)

	UMAPPlot(P6Wk8.3,label=T, label.size=8, pt.size = 2)
	TSNEPlot(P6Wk8.3,label=T, label.size=8, pt.size = 2)
	FeaturePlot(P6Wk8.3, features=LSIG.P6Wk83[["i5"]][1:12], reduction="tsne", label=T)

	FeaturePlot(P6Wk8.3, features=c("Sost", "Bglap", "Ifitm5", "Phex", "Dmp1", 
		"Lpl", "Steap4", "Grem1", "Alpl", "Col4a1"))

P6Wk8.3[["ioneBased3"]]<-paste0("i",formatC(as.integer(Idents(P6Wk8.3)), width=1, flag="0"))
Idents(P6Wk8.3)<-"ioneBased3"
LSIG.P6Wk83<-GETSIG(P6Wk8.3)

table(P6Wk8.3@meta.data$ioneBased3, P6Wk8.3@meta.data$orig.ident)

#####################################################
P6Wk8.4<-subset(P6Wk8.3, idents=c("i1", "i2", "i3", "i4"))

P6Wk8.4<-NormalizeData(P6Wk8.4)
P6Wk8.4<-ScaleData(P6Wk8.4)
P6Wk8.4<-FindVariableFeatures(P6Wk8.4)
P6Wk8.4<-RunPCA(P6Wk8.4)
P6Wk8.4<-RunUMAP(P6Wk8.4, reduction= "cca", dims= 1:10)
P6Wk8.4<-RunTSNE(P6Wk8.4, reduction= "cca", dims.use = 1:40, do.fast = T)
P6Wk8.4<-FindNeighbors(P6Wk8.4, reduction= "cca", dims = 1:10)
P6Wk8.4<-FindClusters(P6Wk8.4, resolution = 0.25)

P6Wk8.4[["ioneBased4"]]<-factor(paste0("i",formatC(as.integer(Idents(P6Wk8.4)), width=1, flag="0")))
Idents(P6Wk8.4)<-"ioneBased4"

LSIG.P6Wk84<-GETSIG(P6Wk8.4)
pdf("umap.p6.annaP56.merge.osteo.pdf", width=8.5)
	UMAPPlot(P6Wk8.4,label=T, label.size=8, pt.size = 2)
	TSNEPlot(P6Wk8.4,label=T, label.size=8, pt.size = 2)
dev.off()

tab1<-table(P6Wk8.4@meta.data$ioneBased3, P6Wk8.4@meta.data$orig.ident)
pieP6<-as.integer(tab1[,1])/sum(tab1[,1])
pieP56<-as.integer(tab1[,2])/sum(tab1[,2])
names(pieP6)<-paste0(rownames(tab1), " (", percent(pieP6), ")")
names(pieP56)<-paste0(rownames(tab1), " (", percent(pieP56), ")")

pdf("piechart.osteogenic.sub-populations.per.timepoint.pdf")
	pie(pieP6, col=hue_pal()(4))
	pie(pieP56, col=hue_pal()(4))
dev.off()

	FeaturePlot(P6Wk8.4, features=c("Sost", "Bglap", "Ifitm5", "Phex", "Dmp1", 
		"Lpl", "Steap4", "Grem1", "Alpl", "Col4a1", "Mki67"), reduction="tsne")

P6Wk8.4[["HCnonHC"]]<-"nonHC"
P6Wk8.4[["HCnonHC"]][which(P6Wk8.4@meta.data$TDT_status),1]<-"HC"
P6Wk8.4[["HCnonHC"]][which(P6Wk8.4@meta.data$WPRE>14),1]<-"HC"

pdf("umap.p6.annaP56.merge.osteo.HCnonHC.pdf", width=8.5)
	UMAPPlot(P6Wk8.4,label=T, label.size=8, pt.size = 2, group.by="HCnonHC")
dev.off()


pdf("featureplit.p6.annaP56.merge.osteo.split.by.timepoints.pdf", width=12, height=3)
	FeaturePlot(P6Wk8.4, features=c("Grem1", "Mmp14", "Pth1r", "Col1a1", "Bglap", "Phex", 
		"Dmp1", "Sost"), split.by="orig.ident", by.col = F, cols=c("grey", "brown"))
dev.off()

tab2<-table(P6Wk8.4@meta.data$orig.ident, P6Wk8.4@meta.data$HCnonHC)


tab2<-table(paste0(P6Wk8.4@meta.data$ioneBased3, "_", P6Wk8.4@meta.data$orig.ident), P6Wk8.4@meta.data$HCnonHC)



pdf("HC.derived.pdf")
	barplot(apply(tab2, 1, function(x)x/sum(x)), las=2)

	ggplot(as.data.frame(tab2)%>%group_by(Var1)%>%
		mutate(Num=sum(Freq), pct=100*Freq/Num, LAB=paste0(round(pct,1),"%")),
		aes(x=1,y=pct,label=LAB,fill=Var2)) + geom_bar(stat="identity", width=1) + coord_polar(theta="y") +
		xlim(c(0.0,1.5)) + facet_wrap(vars(Var1)) + geom_text() +
		scale_fill_brewer(palette="Set1") + theme_void()

dev.off()


LSIG.P6Wk84<-GETSIG(P6Wk8.4)

Idents(P6Wk8.4)<-"ioneBased4"

x1<-as.character(Idents(P6Wk8.4))
x1<-factor(x1, levels=c("i3", "i1", "i4", "i2"))
levels(x1)<-c("i3-prog", "i1-immOb", "i4-cyc", "i2-mOb")
Idents(P6Wk8.4)<-x1
DoHeatmap(P6Wk8.4, features=as.vector(sapply(LSIG.P6Wk84[c(3,1,4,2)], head, 10)))

#####################################################
	Idents(P6Wk8.4)<-"HCnonHC"
	#MARKERS<- FindMarkers(P6Wk8.4, ident.1 = "HC", ident.2 = "nonHC", only.pos = F, min.pct = 0, logfc.threshold = 0.01)
	DEGs<-MARKERS[(MARKERS$p_val_adj<0.05 & abs(MARKERS$avg_log2FC)>0.5), ]
	upDEGs<-DEGs[DEGs$avg_log2FC>0, ]
	downDEGs<-DEGs[DEGs$avg_log2FC<0, ]

pdf("Volcano.HC.vs.nonHC.pdf", width=10, height=10)
	plot(MARKERS$avg_log2FC, -log10(MARKERS$p_val), pch=16, col="grey")
	points(upDEGs$avg_log2FC, -log10(upDEGs$p_val), pch=16, col="red")
	points(downDEGs$avg_log2FC, -log10(downDEGs$p_val), pch=16, col="blue")
	text(DEGs$avg_log2FC, -log10(DEGs$p_val), rownames(DEGs))
dev.off()

#####################################################


matMmpPth<-t(as.matrix( GetAssayData(P6Wk8.4, slot="counts")[c("Mmp14", "Pth1r"),]))

tab3<-table(paste0(matMmpPth[,1]>0, "__", matMmpPth[,2]>0), P6Wk8.4@meta.data$HCnonHC)
pieHC<-as.integer(tab3[,1])/sum(tab3[,1])
pieNonHC<-as.integer(tab3[,2])/sum(tab3[,2])
names(pieHC)<-paste0(rownames(tab3), " (", percent(pieHC), ")")
names(pieNonHC)<-paste0(rownames(tab3), " (", percent(pieNonHC), ")")


pdf("piechart.Mmp14.Pth1r.HC.vs.nonHC.pdf")
	pie(pieHC, col=hue_pal()(6)[2:5])
	pie(pieNonHC, col=hue_pal()(6)[2:5])
dev.off()

#####################################################
#save(P6Wk8, P6Wk8.2, P6Wk8.3, P6Wk8.4, 
	LSIG.HET4, LSIG.P6Wk8, LSIG.P6Wk83, LSIG.P6Wk84, 
	file="P6.P56.Integrate.RData")



