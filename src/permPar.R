permParGlobal <- function(delta1,delta2,B=10,cpu=2){
  each <- B/cpu
  G <- nrow(delta1)
  top <- 500
  
  sfInit(parallel=T,cpus=cpu,type="SOCK")
  
  parFun <- function(xxx) {

    delta1perm <- delta1[sample(1:nrow(delta1),G,replace = F),]
    delta2perm <- delta2[sample(1:nrow(delta2),G,replace = F),]
    permDE1 <- sample(1:G,top,replace = F)
    permDE2 <- sample(1:G,top,replace = F)
    
    out <- round(ARSF(delta1perm,delta2perm,
                      permDE1,permDE2),digits=3)
    return(out)
  }  
  
  sfExport(list=c("delta1","delta2","G","top"))
  
  result<-sfLapply(1:cpu, parFun) 
  sfStop()
  
  out <- unlist(result)
  return(out)
}



permParPathway <- function(delta1,delta2,
                    select.pathways, pathway.size,
                    B=10,cpu=2){
  each <- B/cpu
  K <- length(select.pathways)
  G <- nrow(delta1)
  
  sfInit(parallel=T,cpus=cpu,type="SOCK")
  
  parFun <- function(xxx) {

    delta1perm <- delta1[sample(1:nrow(delta1),nrow(delta1),replace = F),]
    delta2perm <- delta2[sample(1:nrow(delta2),nrow(delta2),replace = F),]
    
    out <- t(replicate(each,sapply(1:K,function(k){
      size <- pathway.size[k]
      index <- sample(1:G,size,replace=F)
      d1.select <- delta1perm[index,]
      d2.select <- delta2perm[index,]
      d1.permDE <- sample(1:length(index),round(0.2*size),replace = F)
      d2.permDE <- sample(1:length(index),round(0.2*size),replace = F)
      
      kout <- round(ARSF(d1.select,d2.select,
                         d1.permDE,d2.permDE),digits=3)
      return(kout)
    },simplify = T)))  ## each x K matrix
    
    colnames(out) <- select.pathways
    return(out)
  }  
  
  sfExport(list=c("delta1","delta2","select.pathways",
                  "pathway.size","G","K","each"))
  
  result<-sfLapply(1:cpu, parFun) 
  sfStop()
  
  out <- do.call(rbind, result)
  return(out)
}




