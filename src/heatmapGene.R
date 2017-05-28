heatmapGene <- function(signPM.list, pathway.list, pathway.name, folder){

  pathway.genes <- pathway.list[[pathway.name]]
  M <- length(signPM.list)
  model.name <- names(signPM.list)
  std.genes <- intersect(names(signPM.list[[1]]),pathway.genes) 
  G <- length(std.genes)
  
  mat <- matrix(0,nrow=G,ncol=M)
  rownames(mat) <- std.genes
  colnames(mat) <- model.name
  
  for (m in 1:M){
    m.genes <- intersect(names(signPM.list[[m]]),std.genes)
    mat[m.genes,m] <- signPM.list[[m]][m.genes]
  }
 
  png(paste(folder, '/Heatmap', pathway.name,'.png',sep=""))
  par(cex.main=1)
  hm<-heatmap.2(mat, symm=F,main=pathway.name,
                cexCol=1,cexRow=0.4,
                colsep=c(1:ncol(mat)),
                rowsep=c(1:nrow(mat)),
                sepwidth=c(0.01, 0.2),  # width of the borders
                sepcolor=c('black'),scale='none',
                symbreaks=T,key=T, keysize=1,symkey=F, 
                dendrogram=c('row'),density.info="none", 
                trace="none",Rowv=T,Colv=F,
                col=bluered,breaks=seq(-1,1,by=0.001),
                srtCol=0,adjCol = c(NA,0.5))
  dev.off()
  
  return(hm)
}    