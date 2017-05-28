orthMatch <- function(geneset1,geneset2,ortholog,reference){
  # column 1 of ortholog corresponds to geneset1
  # column 2 of ortholog corresponds to geneset2
  if(reference=="1") {
    G <- length(geneset1)
    ref <- geneset1
    other <- geneset2
    ref_orth <- ortholog[,1]
    other_orth <- ortholog[,2]
  } else if(reference=="2") {
    G <- length(geneset2) 
    ref <- geneset2
    other <- geneset1
    ref_orth <- ortholog[,2]
    other_orth <- ortholog[,1]
  }

  match_ref <- match_other <- c()
  
  for(g in 1:G){
    if(g%%100==0) print(g)
    ref_gene <-  ref[g]
    if(sum(ref_orth==ref_gene)==0){
      g <- g + 1 
    } else{
      other_orth_gene <- other_orth[which(ref_orth==ref_gene)]
      if(sum(other_orth_gene %in% other)==0){
        g <- g + 1
      } else{
        other_gene <- other_orth_gene[other_orth_gene %in% other][1]  
        match_ref <- c(match_ref,ref_gene)
        match_other <- c(match_other,other_gene)
      }
    }    
  }
  
  if(reference=="1") {
    match1 <- match_ref
    match2 <- match_other
  } else if(reference=="2") {
    match2 <- match_ref
    match1 <- match_other
  }
  return(list(match1,match2))
}