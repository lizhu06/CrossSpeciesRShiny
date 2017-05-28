deSelect <- function(pData, p.cut=NA,topDE.number = 1000){
   pvalue <- pData[,1]
   qvalue <- p.adjust(pvalue,method="BH")
  if (is.na(p.cut)) { 
    DEindex <- which(rank(pvalue,ties.method="first") %in% 1:topDE.number)
    names(DEindex) <- rownames(pData)[DEindex]
  } else {
     DEindex <- which(qvalue < p.cut)
     names(DEindex) <- rownames(pData)[DEindex]
  }
  return(DEindex)
}