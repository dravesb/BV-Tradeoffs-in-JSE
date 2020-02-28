#--------------------------------------
#
#      Determine Critical Value
#       for Omni Semi Par
#
#--------------------------------------

get_T_threshold <- function(X_known, X.scaled, alpha = 0.05){
  
  no.iters <- 500
  T.stat <- numeric(no.iters)
  
  for(i in 1:no.iters){
    #sample A1, A2 
    A1 <- sampP(tcrossprod(X_known))
    A2 <- sampP(tcrossprod(X_known))
    
    #make Omni and embedd for Lhat and Z
    Lhat <- ase(make_omni(list(A1, A2)), 2)
    
    #align
    L_hat <- procrustes(Lhat, X.scaled)$X.new
    
    #make X1 and X2
    X1 <- L_hat[1:(nrow(X_known)),]
    X2 <- L_hat[(nrow(X_known) + 1):(2 * nrow(X_known)),]
    
    #get T value
    T.stat[i] <- sum(apply(X1- X2, 1, function(x) sum(x^2)))
  }
  
  threshold <- unname(quantile(T.stat, probs = 1 - alpha))
  return(threshold)
}

get_T_threshold_parallel <- function(X_known, X.scaled, alpha = 0.05){
  
  no.iters <- 500
  T.stat <- numeric(no.iters)
  
  #set up parallelization
  cores <- detectCores()
  cl <- makeCluster(cores[1] -1)
  registerDoParallel(cl)
  
  #fetch Monte Carlo estimates
  T.stat <- foreach(i = 1:no.iters, 
                    .combine = 'c',
                    .packages = c('MCMCpack')) %dopar% {
    
    #source outside functions
    source('/Users/benjamindraves/Documents/Work/github/BJSE/power_analysis/code/basic_functions.R')
    
    #sample A1, A2 
    A1 <- sampP(tcrossprod(X_known))
    A2 <- sampP(tcrossprod(X_known))
    
    #make Omni and embedd for Lhat and Z
    Lhat <- ase(make_omni(list(A1, A2)), 2)
    
    #align
    L_hat <- procrustes(Lhat, X.scaled)$X.new
    
    #make X1 and X2
    X1 <- L_hat[1:(nrow(X_known)),]
    X2 <- L_hat[(nrow(X_known) + 1):(2 * nrow(X_known)),]
    
    #get T value
    T.stat.here <- sum(apply(X1- X2, 1, function(x) sum(x^2)))
    
    #return T.stat.here
    T.stat.here
  }
  
  #close cluster
  stopCluster(cl)
  
  threshold <- unname(quantile(T.stat, probs = 1 - alpha))
  return(threshold)
}
