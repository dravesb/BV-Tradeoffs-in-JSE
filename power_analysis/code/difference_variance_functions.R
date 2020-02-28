#--------------------------------------
#
#      Variance Functions for 
#       difference (Cor 4)
#
#--------------------------------------

#helper functions
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
Sigma_sum_hat <- function(Xbar, vertex){
  #fetch n, d
  n <- nrow(Xbar)
  d <- ncol(Xbar)
  
  #get inner sum estimate
  total <- diag(0, nrow = d, ncol = 2)
  for(j in (1:n)[-vertex]){
    total <- total + as.numeric(crossprod(Xbar[vertex,], Xbar[j,]) - crossprod(Xbar[vertex,], Xbar[j,])^2) * tcrossprod(Xbar[j,])
  }
  
  #return sum 
  return(1/(n-1) * total)
  
}

#known variance function
Sigma_D <- function(X, S_list, C_list, L.here = L){
  #fetch n,m, and d
  n <- nrow(X)
  m <- length(C_list)
  d <- ncol(X)
  
  #set Delta, community labels, and S2
  Delta <- crossprod(L.here)/d
  comm_labels <- ifelse(sign(X[,2]) > 0, 1, 2)
  S2 <- Reduce('+',(lapply(S_list, function(x) x^2)))
  S <- Reduce('+', S_list)
  
  #set up variance list 
  var.here <- list()
  for(i in 1:n){
    var.here[[i]] <- 0.25 * solve(S2 %*% Delta) %*% S %*% 
      (Sigma_tilde(L.here[comm_labels[i],], C_list[[1]]) + Sigma_tilde(L.here[comm_labels[i],], C_list[[2]])) %*% 
      S %*% solve(S2 %*% Delta)
  }
  
  #return variances for each row of D
  return(var.here)
  
  
}

#estimated variane function
Sigma_D_hat <- function(L_hat){
  #fetch n,m,d
  d <- ncol(L_hat)
  m <- 2
  n <- nrow(L_hat)/m
  
  #set S2_Delta_hat and Omnibar
  S2_Delta_hat <- crossprod(L_hat)/(n*m)
  Omnibar <- 0.5*(L_hat[1:n,] + L_hat[(n+1):(2*n),])
  
  #set up variance list 
  var.here <- list()
  for(i in 1:n){
    var.here[[i]] <- 0.5 * solve(S2_Delta_hat) %*% Sigma_sum_hat(Omnibar,i) %*% solve(S2_Delta_hat)
  }
  
  #return estimated variances for each row of D
  return(var.here)
  
  
}
