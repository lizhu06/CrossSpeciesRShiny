heatmapModel <- function(arsPvaluePath,model.name, cluster.assign,
  pathway.name, folder){
  ## delta.mat is a matrix of pairwise -log(arsPvalue) of M rows and M columns 
  ## with model names 
  ## cluster.assign: results from SA algorithm
  ## pathway.name: pathway name
  
  M <- length(model.name)
  distF <- -log10(arsPvaluePath)
  delta.mat <- matrix(NA,nrow=M,ncol=M)
  rownames(delta.mat) <- colnames(delta.mat) <- model.name 

  for(i in 1:(M-1)){
    for(j in (i+1):M){
      name1 <- rownames(delta.mat)[i]
      name2 <- rownames(delta.mat)[j]
      delta.mat[name1,name2] <- delta.mat[name2,name1] <- distF[paste(name1,name2,sep="")]
    }
  }
  diag(delta.mat) <- max(delta.mat,na.rm=T)
  
  png(paste(folder, "/", pathway.name,'_clustHeatmap.png',sep=""))
  hm <-heatmap.2(delta.mat[order(cluster.assign),order(cluster.assign)], 
                 main=pathway.name,
                 cexCol=1,cexRow=1,
                 colsep=cumsum(table(cluster.assign)),
                 rowsep=cumsum(table(cluster.assign)),
                 sepwidth=c(0.05, 0.05),  # width of the borders
                 sepcolor=c('white'),
                 symbreaks=T,key=T, keysize=1,symkey=F, 
                 dendrogram=c('none'),density.info="none", 
                 trace="none",Rowv=F,Colv=F,
                 srtCol=50, symm=F,
                 col=greenred,breaks=seq(0,round(max(delta.mat)),by=0.01) )
  dev.off()

  return(hm)

}