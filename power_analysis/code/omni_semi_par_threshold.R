#--------------------------------------
#
#      Determine Critical Value
#       for Omni Semi Par
#
#--------------------------------------

get_T_threshold <- function(X_known, alpha = 0.05){
  
  no.iters <- 500
  T.stat <- numeric(no.iters)
  
  for(i in 1:no.iters){
    #sample A1, A2 
    A1 <- sampP(tcrossprod(X_known))
    A2 <- sampP(tcrossprod(X_known))
    
    #make Omni and embedd for Lhat and Z
    Lhat <- ase(make_omni(list(A1, A2)), 2)
    
    #make X1 and X2
    X1 <- Lhat[1:(nrow(X_known)),]
    X2 <- Lhat[(nrow(X_known) + 1):(2 * nrow(X_known)),]
    
    #get T value
    T.stat[i] <- sum(apply(X1- X2, 1, norm2))
  }
  
  threshold <- unname(quantile(T.stat, probs = 1 - alpha))
  return(threshold)
}