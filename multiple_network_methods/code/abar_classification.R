#--------------------------------------
#
#       Abar Classification
#
#--------------------------------------

#source basic function
source("~//Documents/Work/github/BJSE/multiple_network_methods/code/basic_functions.R")

abar_classes <- function(adj_matrices, d, K){
  
  #get Abar
  Abar <- 1/length(adj_matrices) * Reduce('+', adj_matrices)
  
  #get ase
  Xbar <- ase(Abar, d)
  
  #cluster
  #clusters <- kmeans(V, centers = K)
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

