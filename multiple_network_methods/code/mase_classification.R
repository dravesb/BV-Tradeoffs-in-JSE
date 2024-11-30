#--------------------------------------
#
#       MASE Classification
#
#--------------------------------------

#source basic function
source("~//Documents/Work/github/BJSE/multiple_network_methods/code/basic_functions.R")

mase_classes <- function(adj_matrices, d, K){

    #get graph level ASE
    U <- do.call("cbind", lapply(adj_matrices, function(x) ase(x,d)))
    
    #get leading left singular vectors
    V <- svd(U, nu = d, nv = 0)$u
    
    #cluster
    clusters <- Mclust(V, G = K, modelNames = "VVV")$classification
    if(is.null(clusters)){
      clusters <- kmeans(V, centers = K)$cluster
    }
    
    
    # get distance between clusters
    d12 <- get_mahalanobis(V[clusters == 1, ], V[clusters == 2, ])
    d13 <- get_mahalanobis(V[clusters == 1, ], V[clusters == 3, ])
    d23 <- get_mahalanobis(V[clusters == 2, ], V[clusters == 3, ])
    
    #return clusters
    return(list(clusters = clusters, distances = c(d12, d13, d23)))
  
}

