#--------------------------------------
#
#       Omnibus embedding
#
#--------------------------------------

#source basic function
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/basic_functions.R")

omni_classes <- function(adj_matrices, d, K){
  #create omnibus matrix
  Atil <- make_omni(adj_matrices)
  
  #get embedding
  Xhat <- ase(Atil, d)
  
  #get graph parameters
  m <- length(adj_matrices)
  n <- dim(adj_matrices[[1]])[1]
  
  #get 1/m * [I I ... I]^T which is nm x n stack of 1/m*I(nxn) matrices
  Im_stack <- kronecker(rep(1, m), diag(n)/m)
  
  #get sum of each entry
  Xbar <- crossprod(Im_stack, Xhat)
  
  #cluster based on k means
  #clusters <- kmeans(Xbar, centers = K)$cluster
  clusters <- Mclust(Xbar, G = K, modelNames = "VVV")$classification
  
  #return clusters
  return(clusters)
  
}

