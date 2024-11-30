#--------------------------------------
#
#       ASE Classification
#
#--------------------------------------

#source basic function
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/basic_functions.R")

ase_classes <- function(adj_mat, d, K){
  
  #get graph level ASE
  Xhat <- ase(adj_mat, d)
  
  #cluster
  #clusters <- kmeans(V, centers = K)
  clusters <- Mclust(Xhat, G = K, modelNames = "VVV")$classification
  if(is.null(clusters)){
    clusters <- kmeans(Xhat, centers = K)$cluster 
  }
  
  
  # get distance between clusters
  d12 <- get_mahalanobis(Xhat[clusters == 1, ], Xhat[clusters == 2, ])
  d13 <- get_mahalanobis(Xhat[clusters == 1, ], Xhat[clusters == 3, ])
  d23 <- get_mahalanobis(Xhat[clusters == 2, ], Xhat[clusters == 3, ])
  
  #return clusters
  return(list(clusters = clusters, distances = c(d12, d13, d23)))
  
}

