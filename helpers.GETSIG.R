GETSIG<-function(OBJ){
	MARKERS<- FindAllMarkers(OBJ, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
	sapply(levels(MARKERS$cluster),function(i){
		thisDEG<-MARKERS[MARKERS$cluster==i&MARKERS$p_val_adj<0.05,]
		thisDEG<-thisDEG[(thisDEG[,3]-thisDEG[,4])>0.25 | thisDEG[,2]>0.5 ,] #
		thisDEG<-thisDEG[order(thisDEG[,3]-thisDEG[,4],decreasing=T),]
		as.character(thisDEG$gene)
	})
}
GETSIGi<-function(OBJ, GRP1){
	MARKERS<- FindMarkers(OBJ,ident.1=GRP1, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)

	thisDEG<-MARKERS[MARKERS$p_val_adj<0.05,]
	thisDEG<-thisDEG[(thisDEG[,3]-thisDEG[,4])>0.25 | thisDEG[,2]>0.5 ,] #
	thisDEG<-thisDEG[order(thisDEG[,3]-thisDEG[,4],decreasing=T),]
	as.character(rownames(thisDEG))
}
GETDEGs<-function(OBJ){
	MARKERS<- FindAllMarkers(OBJ, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
	flag1<-MARKERS$p_val_adj<0.05
	flag2<-(MARKERS[,3]-MARKERS[,4])>0.25
	flag3<-MARKERS[,2]>0.5
	MARKERS[which(flag1 & flag2 & flag3), , drop=F]
}
SIGcomp<-function(LSIG.A, LSIG.B){
	sapply(LSIG.A, function(x){
		sapply(LSIG.B, function(y){
			length(intersect(x,y))/length(union(x,y))
		})
	})
}
