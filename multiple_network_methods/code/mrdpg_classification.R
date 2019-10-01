#------------------------------------
#
#       Implement MRDPG 
#
#------------------------------------


#source basic function
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/basic_functions.R")

#some helpful functions for MRDPG 
plus <- function(x) max(c(x,0))
get_A_plus <- function(A){
  #get eigendecomposition of A
  eigen_system <- eigen(A)
  
  #threshold values
  D_plus <- diag(sapply(eigen_system$values, plus))
  
  #return A plus
  eigen_system$vectors %*% tcrossprod(D_plus, eigen_system$vectors)
  
}
mrdpg_U <- function(adj_matrices, d){
  
  #get basic parameters
  n <- dim(adj_matrices[[1]])[1]
  m <- length(adj_matrices)
  
  #initialize Lambda_g = I  for g = 1, 2, ..., m 
  Lambda <- rep(list(diag(d)), m)
  
  #Initialize U_old to be an orthogonal n x d matrix
  U_old <- rbind(diag(d), matrix(0, nrow = n - d, ncol = d))
  
  #Get A^g+
  A_plus <- lapply(adj_matrices, get_A_plus)
  
  #set convergence criterion
  conv.crit <- Inf
  
  #begin iteration
  while(conv.crit > 1e-6){
    
    #want to create the sum of ( A_+^g * U_old * Lambda_g)
    sum_mat <- do.call("cbind", lapply(A_plus, function(x) x %*% U_old)) %*% do.call("rbind", Lambda)
    
    #get SVD and set B and C 
    svd_sum_mat <- svd(sum_mat)
    B <- svd_sum_mat$u
    C <- svd_sum_mat$v
    
    #update U
    U <- tcrossprod(B, C) 
    
    #get Z matrices
    Z <- lapply(A_plus, function(x) crossprod(U, x) %*% U)
    
    #update Lambdas
    Lambda_new <- lapply(Z, function(z) diag(sapply(diag(z), plus)))
    
    #update convergence criterion
    conv.crit <- norm(U_old - U, type = "F")/norm(U_old, type = "F")
    
    #update U
    U_old <- U
    
  }
  
  return(U_old)
}


#get classes function
mrdpg_classes <- function(adj_matrices, d, K){
  
  #get commom eigen space
  U <- mrdpg_U(adj_matrices, d)
  
  #kmeans on the rows of U
  classes <- kmeans(U, centers = K)
  
  #return classes
  return(classes$cluster)
}




