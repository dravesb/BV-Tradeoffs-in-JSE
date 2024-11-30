#------------------------------------
#
#       Implement JE
#
#------------------------------------

#source basic function
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/basic_functions.R")

#source JE function
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/joint_embedding.R")

je_classes <-  function(adj_matrices, d, K){
 
  #get H matrix
  H <- multidembed(A = adj_matrices, d, Innitialize = 1, maxiter = 100, large.and.sparse = F)$h
  
  #Cluster based on the rows of H
  clusters <- Mclust(H, G = K, modelNames = "VVV")$classification
  if(is.null(clusters)){
    clusters <- kmeans(H, centers = K)$cluster 
  }
  
  
  # get distance between clusters
  d12 <- get_mahalanobis(H[clusters == 1, ], H[clusters == 2, ])
  d13 <- get_mahalanobis(H[clusters == 1, ], H[clusters == 3, ])
  d23 <- get_mahalanobis(H[clusters == 2, ], H[clusters == 3, ])
  
  #return clusters
  return(list(clusters = clusters, distances = c(d12, d13, d23)))
   
}