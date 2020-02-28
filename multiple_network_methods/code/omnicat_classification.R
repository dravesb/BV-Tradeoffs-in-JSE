#--------------------------------------
#
#       Omnicat classiciation
#
#--------------------------------------

#source basic function
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/basic_functions.R")

omnicat_classes <- function(adj_matrices, d, K){
  #create omnibus matrix
  Atil <- make_omni(adj_matrices)
  
  #get embedding
  Xhat <- ase(Atil, d)
  
  #get graph parameters
  m <- length(adj_matrices)
  n <- dim(adj_matrices[[1]])[1]
  
  #get Xcat
  Xcat <- matrix(NA, nrow = n, ncol = m*d)
  start <- 1; stop <- n
  
  for(i in 1:m){
    #fill Xcat
    Xcat[,c( (d*i - 2):(d*i))] <- Xhat[start:stop,]
    
    #update pointers
    start <- stop + 1; stop <- start + n - 1
    
  }
  
  #cluster based on k means
  clusters <- Mclust(Xcat, G = K, modelNames = "VVV")$classification
  
  #return clusters
  return(clusters)
  
}

