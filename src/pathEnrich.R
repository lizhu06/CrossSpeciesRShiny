pathEnrich <- function(dat, deIndex, pathway.list){
  ## dat is data matrix (either rawData or pData) from the individual model for pathway analysis, deIndex is the DE index
  geneDat <- rownames(dat)
  K <- length(pathway.list)
  pathwayName <- names(pathway.list)
  pathwaySize <- sapply(pathway.list,length)

  bk <- intersect(unique(unlist(pathway.list)),geneDat)
  de <- intersect(unique(unlist(pathway.list)),geneDat[deIndex])
  
  out <- gsa.fisher(x= de, background = bk ,
                    pathway = pathway.list)
  
  overlapSize <- as.integer(out$DE_in_Set)+as.integer(out$NonDE_in_Set)
  DEnum <- as.integer(out$DE_in_Set)
  OR <- (as.numeric(out$DE_in_Set)*as.numeric(out$NonDE_not_in_Set))/
    (as.numeric(out$DE_not_in_Set)*as.numeric(out$NonDE_in_Set))
  OR <- sapply(OR, function(x) {
    if(is.na(x) || is.infinite(x)) {
      x <- 0 } else{
        x <- x
      }
  }, simplify = T)
  
  logOR <- ifelse(OR==0,1,log(OR))
  pvalue <- as.numeric(out$pvalue)
  
  summary <- data.frame(pathway.size= pathwaySize, overlap.size= overlapSize, 
                        DE.inpathway = DEnum, Odds.ratio = OR, 
                        log.odds.ratio = logOR, p.value=pvalue)
  rownames(summary) <- pathwayName
  return(summary)
}


metaPath <- function(pathSummary, pathway.list){

  M <- length(pathSummary)
  K <- length(pathway.list)
  pathwayName <- names(pathway.list)
  pathwaySize <- sapply(pathway.list,length)
  
  pdat <- Reduce("cbind",lapply(pathSummary,function(x) x$p.value))
  rownames(pdat) <- pathwayName
  p.fisher <- apply(pdat,1,fisher)
  q.fisher <- p.adjust(p.fisher,method="BH")
  
  overlapMinSize <- apply(Reduce("cbind",lapply(pathSummary,
                                                function(x) x$overlap.size)),
                          1,min)
  
  overlapMedDE <-  apply(Reduce("cbind",lapply(pathSummary,
                                               function(x) x$DE.inpathway)),1,
                         function(x) (
                           if(M %% 2 ==0){
                             (sort(x)[M/2]+sort(x)[(M/2+1)])/2
                           } else { sort(x)[(M+1)/2] } ))
  
  meta.summary <- data.frame(pathway.size= pathwaySize, 
                             min.overlap.size= overlapMinSize, 
                             med.DE.inpathway = overlapMedDE, 
                             p.fisher=as.numeric(p.fisher), 
                             q.fisher=as.numeric(q.fisher))
  rownames(meta.summary) <- pathwayName
  return(meta.summary)
  
}


pathSelect <- function(meta.summary,pathwaysize.lower.cut=5,
                       pathwaysize.upper.cut=200,
                       overlapsize.cut=5,med.de.cut=3,
                       qfisher.cut = 0.05){
  select.ind <- which(meta.summary$pathway.size >= pathwaysize.lower.cut & 
                        meta.summary$pathway.size <= pathwaysize.upper.cut &
                        meta.summary$min.overlap.size >=overlapsize.cut &
                        meta.summary$med.DE.inpathway >=med.de.cut &
                        meta.summary$q.fisher < qfisher.cut) 
  select.pathways <- rownames(meta.summary)[select.ind]
  return(select.pathways)
}

