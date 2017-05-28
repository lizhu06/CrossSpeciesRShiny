clustPathway <- function(arsPvalue) {
  #require(ConsensusClusterPlus)
  #set your working dir, automatically save there
  
  results = ConsensusClusterPlus(d=t(-log10(arsPvalue)),maxK=10,reps=50,pItem=0.8,pFeature=1,title="Pathway clustering",clusterAlg="hc",innerLinkage="ward.D2",finalLinkage="ward.D2",seed=15213,plot="png")
  
  return(results) ## a list of K elements (each k represents number of clusters)
}