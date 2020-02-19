#--------------------------------------------------------------
#            Theoretical MSE Functions
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
#       Variance Functions
#-----------------------------

#helper variance functions
Sigma_tilde <- function(c, p){
  p^2*c - p^3*c^2
}
Sigma_sum_omni <- function(S_list, c, graph){
  #set p 
  p <- 0.5
  m <- length(S_list)
  
  #define mSbar
  S <- sum(unlist(S_list))
  
  #calculate variance
  V1 <- (S_list[[graph]] + S) * Sigma_tilde(c, p) * (S_list[[graph]] + S)
  V2 <- 0
  for(g in (1:m)[-graph]){
    V2 <- V2 + S_list[[g]] * Sigma_tilde(c,p) * S_list[[g]] 
  }
  
  #return variance
  return(V1 + V2)
  
}
Sigma_sum_omnibar <- function(S_list, C_list, Sbar){
  #set p 
  p <- 0.5
  #m <- length(S_list)
  
  #calculate variance
  mat <- 0
  for(g in 1:m){
    mat <- mat + (Sbar + S_list[[g]]) * Sigma_tilde(C_list[[g]], p) * (Sbar + S_list[[g]])
  }
  
  #return variance
  return(mat)
  
}

#define variance functions by method
ase_var <- function(X, C_list){
  #fetch n, m, and d 
  n <- nrow(X)
  d <- ncol(X)
  m <- length(C_list)
  
  #set p 
  p <- 0.5
  
  #calculate variance
  Delta <- p
  var.here <- matrix(rep(sapply(1:m, 
                                function(graph) (Delta * sqrt(C_list[[graph]]))^(-1) * Sigma_tilde(C_list[[graph]], p) * (Delta * sqrt(C_list[[graph]]))^(-1)), 
                                #function(graph) 1 - C_list[[graph]]^2*p), 
                         each = n),
                     ncol = d, nrow = n*m)
  #return variance
  return(var.here)
}
abar_var <- function(X, C_list){
  #fetch n, m, and d
  m <- length(C_list)
  n <- nrow(X)
  d <- ncol(X)
  
  #set p and m
  p <- 0.5
  
  #calcualte Delta and C bar
  Delta <- p
  Cbar <- Reduce('+', C_list)/m
  
  #calculate variance
  var.here <- matrix(rep(sapply(1:m, 
                                function(graph) 0.25 * (sqrt(Cbar)*Delta)^(-1) * (Sigma_tilde(C_list[[1]], p) +  Sigma_tilde(C_list[[2]], p)) * (sqrt(Cbar)*Delta)^(-1)), 
                                #function(graph) 1 - Cbar2^2*p), 
                         each = n),
                     nrow = n*m, ncol = d
  )
                     
  
  #return variance
  return(var.here)
  
}
omni_var <- function(X, C_list, S_list){
  #fetch n,m, and d
  n <- nrow(X)
  m <- length(C_list)
  d <- ncol(X)
  
  #set p 
  p <- 0.5

  #calcualte S^2 and Delta
  Delta <- p
  S2 <- sum(unlist(lapply(S_list, function(x) x^2)))
  
  #calculate variance
  var.here <- matrix(rep(sapply(1:m , function(graph) 0.25 * (Delta * S2)^(-1) * Sigma_sum_omni(S_list, C_list[[graph]], graph) * (Delta * S2)^(-1)), 
                         each = n),
                     nrow = n*m, ncol = d)
  
  #return variance
  return(var.here)
}
omnibar_var <- function(X, C_list, S_list){
  #fetch n, m, and d
  n <- nrow(X)
  m <- length(C_list)
  d <- ncol(X)
  
  #set p
  p <- 0.5
  
  #calcualte S^2 and Delta
  Delta <- p
  S2 <- Reduce('+', lapply(S_list, function(x) x^2))
  Sbar <- Reduce('+', S_list)/m
  
  #calculate variance
  var.here <- matrix(rep(sapply(1:m , 
                                function(graph) 0.25 * (Delta * S2)^(-1) * Sigma_sum_omnibar(S_list, C_list, Sbar) * (Delta * S2)^(-1)), 
                         each = n),
                     nrow = n*m, ncol = d)
  
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
  bias2.here <- apply(ase_bias(X, C_list), 1, function(x) x^2)
  
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
  bias2.here <- apply(abar_bias(X, C_list), 1, function(x) x^2)
  
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
  bias2.here <- apply(omni_bias(X, C_list, S_list), 1, function(x) x^2)
  
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
  bias2.here <- apply(omnibar_bias(X, C_list, S_list), 1, function(x) x^2)
  
  #get variance matrix
  var.here <- apply(omnibar_var(X, C_list, S_list), 1, function(x) x/n)
  
  #return MSE matrix 
  return(bias2.here + var.here)
  
}