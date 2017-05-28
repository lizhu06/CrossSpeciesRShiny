textMine <- function(hashtb,pathways,cluster.assign){
  ##cluster.assign with pathway names (w/o scatterness)
  
  tmk <- TextMine(hashtb=hashtb, pathways= pathways,
                  pathway=names(cluster.assign), result=cluster.assign)
  C <- length(unique(cluster.assign))
  tm_filtered <- list()
  for (i in 1:C){ 
    tm_filtered[[i]] <- tmk[[i]][which((as.numeric(tmk[[i]][,4]) < 0.05)), ]
  }
  
  pathway.summary <- lapply(1:C, function(x) names(which(cluster.assign==x)))

  dir.path <- "textmine"
  if (!file.exists(dir.path)) dir.create(dir.path)
  setwd(dir.path)
  writeTextOut(tm_filtered,C,pathway.summary)
  setwd("../")
}


