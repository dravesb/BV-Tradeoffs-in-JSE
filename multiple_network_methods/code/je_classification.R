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
  classes <- Mclust(H, G = K, modelNames = "VVV")$classification
  
  #return clusters
  return(classes)
   
}