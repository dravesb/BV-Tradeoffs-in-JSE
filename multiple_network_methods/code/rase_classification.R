#--------------------------------------
#
#       RASE Classification
#
#--------------------------------------

#source basic function
source("~//Documents/Work/github/BJSE/multiple_network_methods/code/basic_functions.R")

uase_ise <- function(A_list, d, diag.augment = T, ...){
  # diagonally augment
  if(diag.augment & Reduce('+', lapply(A_list, function(A) sum(abs(diag(A))))) == 0) {
    A_list_aug <- lapply(A_list, function(A){
      deg <- colSums(A)
      n <- ncol(A)
      diag(A) <- deg/(n-1)
      A
    })
  }
  
  # construct A and decompose
  A_hat <- Reduce('cbind', A_list_aug)
  A_hat.svd <- irlba::irlba(A_hat, nv = d, nu = 0)
  
  # RASE 
  #X_rase <- sqrt(sqrt(length(A_list))) * A_hat.svd$v %*% diag(sqrt(A_hat.svd$d[1:d]))
  X_rase <- length(A_list)^(0.25) * A_hat.svd$v %*% diag(sqrt(A_hat.svd$d[1:d]))
  
  #break into graph latent positions
  X_list <- list(); n <- nrow(A_list[[1]])
  for(g in 1:length(A_list)){
    start <- n*(g-1) + 1
    stop <- n*g
    X_list[[g]] <- as.matrix(X_rase[start:stop, ], nrow = n)
  }
  return(list(X_list = X_list))
}

rase_classes <- function(adj_matrices, d, K){
  
  #get Abar
  Xbar <- Reduce('+', uase_ise(adj_matrices, d)$X_list)/length(adj_matrices)
  
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

