bayesP <- function(pData, seed=12345){
  p <- pData[,1]
  lfc <- pData[,2]
  G <- nrow(pData)
  z <- PtoZ(p,lfc)
  iteration <- 2000
  burnin <- 1000
  S <- 500
  names(z) <- rownames(pData)
  prop <- SelectGamma(p)
  if(prop <= 0.3) {
     gamma <- G*0.3
  } else{
     gamma <- G*prop
  }
  MCMCout <- MCMC(z, iteration, gamma) 
  signdelta <- MCMCout$Y[,-c(1:burnin)]
  
  set.seed(seed)
  signdelta.sub <- signdelta[,sample(1:ncol(signdelta),S)]
  
  rownames(signdelta.sub) <- names(z)
  return(signdelta.sub) #the full sign delta for a dataset (subsample 500)
}  