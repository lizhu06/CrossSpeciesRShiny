heatmapPathway <- function(arsPvalue, cluster.assign,scatter.index=NULL){
  
  ## cluster.assign is the result from consensus clustering
  
  C <- length(unique(cluster.assign)) #number of clusters
  if(!is.null(scatter.index)){
    cluster.assign[scatter.index] <- -1
  }
  
  dataOrder <- -log10(arsPvalue)[unlist(sapply(c(1:C,-1),function(x) which(cluster.assign==x))),]
  colnames(dataOrder) <- colnames(arsPvalue)
  #rownames(dataOrder) <- sapply(rownames(dataOrder),function(x) substr(x,1,15))
  row.colors <- sapply(1:cluster.assign,function(x) {
                        ifelse(x==-1,"gray",topo.colors(C)[cluster.assign[x]])
                     })
  
  par(cex.main=1, font.lab=2, font.axis=2)
  
  pdf(paste("Heatmap_Pathway",C,"clusters.pdf",sep="_"))
  hm<-heatmap.2(dataOrder, symm=F,main=NULL,
                cexCol=0.6,cexRow=0.2,adjCol= c(NA,-1),
                rowsep=c(0,cumsum(unlist(sapply(c(1:C,-1),function(x) sum(cluster.assign.cc==x)))) ),
                sepwidth=c(0.1, 0.3),  # width of the borders
                sepcolor=c('white'),scale='none',
                symbreaks=T,key=T, keysize=1,symkey=F, 
                dendrogram=c('none'),density.info="none", 
                trace="none",Rowv=F,Colv=T,
                srtCol=50,RowSideColors=row.colors,
                col=greenred,breaks=seq(0,max(dataOrder),by=0.01),
                key.ytickfun=function(){
                  side = 2
                } )
  dev.off()
  
  return(hm)

}  