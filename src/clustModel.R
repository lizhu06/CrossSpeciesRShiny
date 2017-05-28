clustModel <- function(arsPvaluePath,model.name,Tm=10,P=0.5,
                       mu=0.9,epsilon=1e-5,N=1000,seed=12345){
  ## Total possible configurations = choose(n,1) + choose(n,2) + ...
  ## clustering models using SA algorithm
  ## delta.mat is a matrix of pairwise -log10(arsPvalue) of M rows and M columns 
  ## with model names
  M <- length(model.name)
  distF <- -log10(arsPvaluePath)
  delta.mat <- matrix(NA,nrow=M,ncol=M)
  rownames(delta.mat) <- colnames(delta.mat) <- model.name 

  for(i in 1:(M-1)){
    for(j in (i+1):M){
      name1 <- rownames(delta.mat)[i]
      name2 <- rownames(delta.mat)[j]
      delta.mat[name1,name2] <- delta.mat[name2,name1] <- distF[paste(name1,name2,sep="")]
    }
  }
  diag(delta.mat) <- max(delta.mat,na.rm=T)
    
  n <- nrow(delta.mat) # =M
  obs.name <- rownames(delta.mat)
  
  ## initialize
  hc <- hclust(d=dist((max(delta.mat)-delta.mat)^2))
  K <- 3
  a <- cutree(hc, k = K)
  names(a) <- obs.name
  delta.est <- Est_mean(delta.mat,a) 
  names(delta.est) <- c(0,1:K)
  
  count <- 0 #initial
  Jc <- E_tot(delta.mat,a,delta.est)
  pi <- exp(-Jc/Tm) ## Boltzmann dist
  
  while((count < N) && (Tm >= epsilon)) {
    ##New trial
    u <- runif(1)
    if(u>0.5){
      a_new <- Split(a)
    } else{
      a_new <- Relocate(a)
    }
    names(a_new) <- obs.name
    K_new <- length(unique(a_new))
    delta.est_new <- Est_mean(delta.mat,a_new)
    Jn <- E_tot(delta.mat, a_new, delta.est_new)
    pi_new <- exp(-Jn/Tm)
    
    if(Jn < Jc) { 
      ##accept
      Jc <- Jn;
      a <- a_new;
      K <- K_new;
      delta.est <- delta.est_new;
    } else {
      count <- count + 1;
      r <- min(1,pi_new/pi); ## acceptance prob.
      u <- runif(1);
      if(u>r) {
        ##not accept
        Tm <- Tm*mu
      } else {
        Jc <- Jn;
        a <- a_new;
        K <- K_new;
        delta.est <- delta.est_new;
      } 
    }
  }
  return(a) #cluster assigned 
}
