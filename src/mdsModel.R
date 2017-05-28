mdsModel <- function(arsPath,model.name,pathway.name) {
  ## for each pathway, plot MDS for all models
  ## arsP is a vector of choose(M,2) elements (named in paste(name1,name2,sep="")) 
  ## model.name is a vector of model names
  M <- length(model.name)
  distF <- ARStransform(arsPath)
  d <- matrix(NA,nrow=M,ncol=M)
  rownames(d) <- colnames(d) <- model.name
  
  for(i in 1:(M-1)){
    for(j in (i+1):M){
      name1 <- rownames(d)[i]
      name2 <- rownames(d)[j]
      d[name1,name2] <- d[name2,name1] <- distF[paste(name1,name2,sep="")]
    }
  }
  
  diag(d) <- 0
  dist <- as.dist(d,upper = TRUE, diag = TRUE)
  fit <- sammon(d=dist, y= jitter(cmdscale(dist, 2)), k=2) # k is the number of dim
  
  x <- fit$points[,1]
  y <- fit$points[,2]
  xlimit <- ifelse(abs(min(x))>abs(max(x)),abs(min(x)),abs(max(x)))
  ylimit <- ifelse(abs(min(y))>abs(max(y)),abs(min(y)),abs(max(y)))
  
  color <- rainbow(M,s=0.5,v=1,alpha=1)
  
    p <-  ggplot() +
          ggtitle(pathway.name) +
          xlab("Coordinate 1") + ylab("Coordinate 2") + 
          xlim(c(-xlimit-0.5,xlimit+0.5)) + ylim(c(-ylimit-0.5,ylimit+0.5)) + 
          geom_point(aes(x, y), color = color  ,size=6) +
          geom_text_repel(aes(x, y, label = rownames(d),fontface="bold"),size=8) + 
          theme(plot.title = element_text(size = 15, hjust=0.5,face="bold"),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12))
   print(p)
  return(p)
}


