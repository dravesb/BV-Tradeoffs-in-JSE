#--------------------------------------------------------------
#            Theoretical MSE Functions - 2D
#
#
#     - Bias for: ASE, Abar, Omni, Omnibar
#     - Variance for: ASE, Abar, Omni, Omnibar 
#
#
#--------------------------------------------------------------

#------------------------------
#       Bias Functions
#-----------------------------
ase_bias <- function(X, C_list){
  #get n and m
  n <- nrow(X)
  m <- length(C_list)
  d <- nrow(C_list[[1]])
  
  #return bias
  return(matrix(0, nrow = n*m, ncol = d))
  
}
abar_bias <- function(X, C_list){
  #calculate Cbar
  m <- length(C_list)
  Cbar <- Reduce('+', C_list)/m
  
  #create bias matrix
  return(Reduce(rbind, lapply(C_list, function(C)  X %*% (sqrt(Cbar) - sqrt(C)))))
}
omni_bias <- function(X, C_list, S_list){
  #fetch m 
  m <- length(S_list)
  
  #create bias matrix
  return(Reduce(rbind, lapply(1:m, function(g)  X %*% (S_list[[g]] - sqrt(C_list[[g]])))))
}
omnibar_bias <- function(X, C_list, S_list){
  #fetch m and d
  m <- length(S_list)
  
  #Set up Sbar matrix 
  Sbar <- Reduce('+', S_list)/m
  
  #create bias matrix
  return(Reduce(rbind, lapply(1:m, function(g)  X %*% (Sbar - sqrt(C_list[[g]])))))
}

#------------------------------
#     Trace Variance Functions
#-----------------------------
#helper variance functions
Sigma_tilde <- function(z, C, L.here = L){
  #fetch latent positions
  l1 <- L.here[1,]; l2 <- L.here[2,]; d <- ncol(L)
  
  #set up constants
  c1 <- tcrossprod(crossprod(z, C), l1) - (tcrossprod(crossprod(z, C), l1))^2
  c2 <- tcrossprod(crossprod(z, C), l2) - (tcrossprod(crossprod(z, C), l2))^2
  
  #set up matrices
  mat1 <- matrix(c(l1[1]^2, l1[1]*l1[2], l1[1]*l1[2], l1[2]^2),
                 ncol = d, nrow = d)
  mat2 <- matrix(c((c1+c2)/2, (c1-c2)/2, (c1-c2)/2, (c1+c2)/2),
                 ncol = d, nrow = d)
  
  #return sigma_tilde
  return(mat1*mat2)
  
}
Sigma_sum_omni <- function(S_list, C_list, graph, comm, L.here = L){
  #fetch m
  m <- length(S_list)
  
  #define mSbar
  S <- Reduce('+', S_list)
  
  #calculate variance
  V1 <- (S_list[[graph]] + S) %*% Sigma_tilde(L.here[comm,], C_list[[graph]]) %*% (S_list[[graph]] + S)
  V2 <- 0
  for(g in (1:m)[-graph]){
    V2 <- V2 + S_list[[g]] %*% Sigma_tilde(L.here[comm,], C_list[[g]]) %*% S_list[[g]] 
  }
  
  #return variance
  return(V1 + V2)
  
}
Sigma_sum_omnibar <- function(S_list, C_list, comm, L.here = L){
  
  #set Sbar
  m <- length(S_list)
  Sbar <- Reduce('+', S_list)/m
  
  #calculate variance
  mat <- 0
  for(g in 1:m){
    mat <- mat + (Sbar + S_list[[g]]) %*% Sigma_tilde(L.here[comm,], C_list[[g]]) %*% (Sbar + S_list[[g]])
  }
  
  #return variance
  return(mat)
  
}

#define trace variance functions by method
ase_var <- function(X, C_list, L.here = L){
  #fetch n, m, and d 
  n <- nrow(X)
  d <- ncol(X)
  m <- length(C_list)
  
  #set Delta and community labels
  Delta <- diag(L.here[1,]^2)
  comm_labels <- ifelse(sign(X[,2]) > 0, 1, 2)
  
  #get trace of community variance
  var.here <- matrix(NA, ncol = 1, nrow = n*m)
  for(graph in 1:m){
    for(vertex in 1:n){
      var.here[n*(graph - 1) + vertex, 1] <- sum(diag(solve(Delta %*% sqrt(C_list[[graph]])) 
                                                      %*% Sigma_tilde(L.here[comm_labels[vertex],], C_list[[graph]], L.here) 
                                                      %*% solve(Delta %*% sqrt(C_list[[graph]]))))
    }
  }
  
  #return variance
  return(var.here)
}
abar_var <- function(X, C_list, L.here = L){
  #fetch n, m, and d
  m <- length(C_list)
  n <- nrow(X)
  d <- ncol(X)
  
  #set Delta, community labels, and Cbar
  Delta <- diag(L.here[1,]^2)
  comm_labels <- ifelse(sign(X[,2]) > 0, 1, 2)
  Cbar <- Reduce('+', C_list)/m
  
  #get trace of community variance
  var.here <- matrix(NA, ncol = 1, nrow = n*m)
  for(graph in 1:m){
    for(vertex in 1:n){
      var.here[n*(graph - 1) + vertex] <- sum(diag(0.25 * solve(sqrt(Cbar) %*% Delta) 
                                                   %*% (Sigma_tilde(L.here[comm_labels[vertex],], C_list[[1]], L.here) + 
                                                          Sigma_tilde(L.here[comm_labels[vertex],], C_list[[2]], L.here)) 
                                                   %*% solve(sqrt(Cbar) %*% Delta)
                                                   ))
      
    }
  }
  
                     
  
  #return variance
  return(var.here)
  
}
omni_var <- function(X, C_list, S_list, L.here = L){
  #fetch n,m, and d
  n <- nrow(X)
  m <- length(C_list)
  d <- ncol(X)
  
  #set Delta, community labels, and S2
  Delta <- crossprod(L.here)/d
  comm_labels <- ifelse(sign(X[,2]) > 0, 1, 2)
  S2 <- Reduce('+',(lapply(S_list, function(x) x^2)))
  
  #get trace of community variance
  var.here <- matrix(NA, ncol = 1, nrow = n*m)
  for(graph in 1:m){
    for(vertex in 1:n){
      var.here[n*(graph - 1) + vertex] <- 0.25 * sum(diag(solve(Delta %*% S2) 
                                                   %*% Sigma_sum_omni(S_list, C_list, graph, comm_labels[vertex], L.here) 
                                                   %*% solve(Delta %*% S2)
                                                   ))
    }
  }
  
  #return variance
  return(var.here)
}
omnibar_var <- function(X, C_list, S_list, L.here = L){
  #fetch n, m, and d
  n <- nrow(X)
  m <- length(C_list)
  d <- ncol(X)
  
  #set Delta, community labels, and S2
  Delta <- crossprod(L.here)/2
  comm_labels <- ifelse(sign(X[,2]) > 0, 1, 2)
  S2 <- Reduce('+', lapply(S_list, function(x) x^2))
  
  #get trace of community variance
  var.here <- matrix(NA, ncol = 1, nrow = n*m)
  for(graph in 1:m){
    for(vertex in 1:n){
      var.here[n*(graph - 1) + vertex] <- sum(diag(0.25 * solve(Delta %*% S2) 
                                                   %*% Sigma_sum_omnibar(S_list, C_list, comm_labels[vertex]) 
                                                   %*% solve(Delta %*% S2)
                                                   ))
    }
  }
  
  #return variance
  return(var.here)
  
}

#------------------------------
#       MSE Functions
#-----------------------------

#define MSE functions by method
ase_mse <- function(n, X, C_list){
  #fetch m
  m <- length(C_list)
  d <- ncol(X)
  n <- nrow(X)
  
  #get bias matrix
  bias2.here <- apply(ase_bias(X, C_list), 1, function(x) sum(x^2))
  
  #get variance matrix
  var.here <- apply(ase_var(X, C_list), 1, function(x) x/n)
  
  #return MSE matrix 
  return(bias2.here + var.here)
  
}
abar_mse <- function(n, X, C_list){
  #fetch n,m,d
  m <- length(C_list)
  d <- ncol(X)
  n <- nrow(X)
  
  #get bias matrix
  bias2.here <- apply(abar_bias(X, C_list), 1, function(x) sum(x^2))
  
  #get variance matrix
  var.here <- apply(abar_var(X, C_list), 1, function(x) x/n)
  
  #return MSE matrix 
  return(bias2.here + var.here)
  
}
omni_mse <- function(n, X, C_list, S_list){
  #fetch n,m,d
  m <- length(C_list)
  d <- ncol(X)
  n <- nrow(X)
  
  #get bias matrix
  bias2.here <- apply(omni_bias(X, C_list, S_list), 1, function(x) sum(x^2))
  
  #get variance matrix
  var.here <- apply(omni_var(X, C_list, S_list), 1, function(x) x/n)
  
  #return MSE matrix 
  return(bias2.here + var.here)
}
omnibar_mse <- function(n, X, C_list, S_list){
  #fetch n,m,d
  m <- length(C_list)
  d <- ncol(X)
  n <- nrow(X)
  
  #get bias matrix
  bias2.here <- apply(omnibar_bias(X, C_list, S_list), 1, function(x) sum(x^2))
  
  #get variance matrix
  var.here <- apply(omnibar_var(X, C_list, S_list), 1, function(x) x/n)
  
  #return MSE matrix 
  return(bias2.here + var.here)
  
}