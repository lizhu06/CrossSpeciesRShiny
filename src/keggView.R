keggView <- function(mat,kegg.name,folder){
  ## human only 
  #set your working dir, automatically save there
  #mat is two column signed PM
  #kegg.name <- gsub("KEGG ","",pathway.name)  
  #database = c("org.Hs.eg.db")
  #also need library("biomaRt") library("KEGG.db")

  #require(database)
  kegg.name1 <- gsub("KEGG ","",kegg.name) 
  xx <- unlist(as.list(KEGGPATHID2NAME))
  pathwayID <- names(xx)[which(xx==kegg.name1)]

  model.name <- colnames(mat)
  std.genes <- rownames(mat)
  
  out <- unlist(mget(x=std.genes,envir=org.Hs.egALIAS2EG))
  IDs <- sapply(std.genes,function(x) out[grep(x, names(out))[1]]  ) 
  geneDat <- mat
  rownames(geneDat) <- IDs
  
  current_dir <- getwd()
  dir.path <- paste(folder, "/topo_", kegg.name, sep="")
  if (!file.exists(dir.path)) dir.create(dir.path)
  setwd(dir.path)
  pv.out <- pathview(gene.data = geneDat, pathway.id = pathwayID, species = "hsa", out.suffix = "", kegg.native = T, key.pos = "bottomright", map.null=T,cex = 0.15)
  setwd(current_dir)
  return(pv.out)
  
}    