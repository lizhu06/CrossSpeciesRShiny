permGlobal <- function(delta1,delta2,B=10){
  G <- nrow(delta1)
  top <- 500
  out <- rep(NA,B)
  for(b in 1:B){
    delta1perm <- delta1[sample(1:nrow(delta1),G,replace = F),]
    delta2perm <- delta2[sample(1:nrow(delta2),G,replace = F),]
    permDE1 <- sample(1:G,top,replace = F)
    permDE2 <- sample(1:G,top,replace = F)
    out[b] <- round(ARSF(delta1perm,delta2perm,
                      permDE1,permDE2),digits=3)
  }  
  return(out)
}



permPathway <- function(delta1,delta2,
                    select.pathways, pathway.size,B=10){
  K <- length(select.pathways)
  G <- nrow(delta1)
  
  out <- matrix(NA,nrow=B,ncol=K)
  colnames(out) <- select.pathways

  for(b in 1:B){
    delta1perm <- delta1[sample(1:nrow(delta1),nrow(delta1),replace = F),]
    delta2perm <- delta2[sample(1:nrow(delta2),nrow(delta2),replace = F),]
    
    out[b,] <- sapply(1:K,function(k){
      size <- pathway.size[k]
      index <- sample(1:G,size,replace=F)
      d1.select <- delta1perm[index,]
      d2.select <- delta2perm[index,]
      if(length(index)<=5){
      	  d1.permDE <- d2.permDE <- 1:length(index)
        } else {
          d1.permDE <- sample(1:length(index),round(0.5*size),replace = F)
          d2.permDE <- sample(1:length(index),round(0.5*size),replace = F)
      }
      kout <- round(ARSF(d1.select,d2.select,
                         d1.permDE,d2.permDE),digits=3)
      return(kout)
    },simplify = T)  ## each x K matrix
  }  

  return(out)
}




