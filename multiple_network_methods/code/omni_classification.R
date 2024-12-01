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
  clusters <- Mclust(Xbar, G = K, modelNames = "VVV")$classification
  if(is.null(clusters)){
    clusters <- kmeans(Xbar, centers = K)$cluster 
  }
  
  # get distance between clusters
  d12 <- get_mahalanobis(Xbar[clusters == 1, ], Xbar[clusters == 2, ])
  d13 <- get_mahalanobis(Xbar[clusters == 1, ], Xbar[clusters == 3, ])
  d23 <- get_mahalanobis(Xbar[clusters == 2, ], Xbar[clusters == 3, ])
  
  #return clusters
  return(list(clusters = clusters, distances = c(d12, d13, d23)))
  
  
}

