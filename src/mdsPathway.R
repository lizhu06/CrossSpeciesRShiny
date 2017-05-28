mdsPathway <- function(arsPvalue,cluster.assign,scatter.index=NULL,folder) {

  ## plot MDS for all pathways
  ## arsPvalue is a matrix of K (pathways) rows and choose(M,2) columns 
  ## cluster.assign is the result from consensus clustering
  
  C <- length(unique(cluster.assign)) #number of clusters
  if(C>9){
    warning("Too many clusters, not enough colors")
  }
  
  dist.mat <-  dist(-log10(arsPvalue),method = "euclidean",
                    upper = TRUE, diag = TRUE)
  
  fit <- cmdscale(dist.mat,k=2)
  x <- fit[,1]
  y <- fit[,2]
  xlimit <- ifelse(abs(min(x))>abs(max(x)),abs(min(x)),abs(max(x)))
  ylimit <- ifelse(abs(min(y))>abs(max(y)),abs(min(y)),abs(max(y)))
  xcenter <- tapply(x,as.factor(cluster.assign),mean)
  ycenter <- tapply(y,as.factor(cluster.assign),mean)
  
  if(!is.null(scatter.index)){
    cluster.assign[scatter.index] <- -1
  }

  unique.color <- rainbow(C,s=0.5,v=1,alpha=1)
  unique.shape <- 1:C
  sizes <- shapes <- colors <- cluster.assign
  for(i in 1:(C+1)){
    colors[cluster.assign==i] <- unique.color[i]
    shapes[cluster.assign==i] <- unique.shape[i]
    sizes[cluster.assign==i] <- 2
    if(i== (C+1)){
      colors[cluster.assign== -1] <- "gray50"
      shapes[cluster.assign== -1] <- 20 
      sizes[cluster.assign== -1] <- 1
    }
  }
  
  png(paste(folder, "/MDS_pathway_",C,"_clusters.png",sep=""))
  p <-    ggplot() +
          ggtitle("") +
          xlab("Coordinate 1") + ylab("Coordinate 2") + 
          xlim(c(-xlimit,xlimit)) + ylim(c(-ylimit,ylimit)) + 
          geom_point(aes(x, y), shape=shapes, 
                     color = colors ,size=sizes) +
          geom_point(aes(xcenter,ycenter),
                     shape=unique.shape, color = unique.color, 
                     size =5) + 
          theme_bw() + 
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                plot.title = element_text(size = 15, hjust=0.5,face="bold"),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12))
  print(p)
  dev.off()
  
  return(p)
  
}  
  