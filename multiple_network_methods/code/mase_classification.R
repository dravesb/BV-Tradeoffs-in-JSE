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
    clusters <- kmeans(V, centers = K)
    
    #return clusters
    return(clusters$cluster)
  
}

