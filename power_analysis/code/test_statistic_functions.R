#--------------------------------------
#
#      Test Statistics Functions 
#
#--------------------------------------

#Semi Par Test statistic
get_T_stat <- function(L_hat){
  #fetch n, m ,d
  d <- ncol(L_hat)
  m <- 2
  n <- nrow(L_hat)/m
  
  #get D_hat
  D_hat <- L_hat[1:n, ] - L_hat[(n+1):(2*n), ]
  
  #get test statistic
  T_stat <- sum(apply(D_hat, 1, function(x) sum(x^2)))
  
  #return statistic
  return(T_stat)
}

#Known Variance Test statistic
get_W_stat <- function(L_hat, X, S_list, C_list){
  #fetch n, m ,d
  d <- ncol(L_hat)
  m <- 2
  n <- nrow(L_hat)/m
  
  #get D_hat
  D_hat <- L_hat[1:n, ] - L_hat[(n+1):(2*n), ]
  
  #fetch vertex level variance
  var.here <- Sigma_D(X, S_list, C_list)
  
  #unstable inversion
  var.here.invert <- lapply(var.here, function(x) solve(x))
  
  #get vertex level test stat
  Wi_list <- lapply(1:n, function(i) n*crossprod(D_hat[i,], var.here.invert[[i]]) %*% D_hat[i,])
  
  #get and return W 
  W <- as.numeric(Reduce('+', Wi_list))
  return(W)
  
}

#Unknown Variance Test statistic
get_W_hat_stat <- function(L_hat){
  #fetch n, m ,d
  d <- ncol(L_hat)
  m <- 2
  n <- nrow(L_hat)/m
  
  #get D_hat
  D_hat <- L_hat[1:n, ] - L_hat[(n+1):(2*n), ]
  
  #fetch vertex level variance
  var.here <- Sigma_D_hat(L_hat)
  
  #project onto psd cone and then psuedo inverse
  var.here.invert <- lapply(var.here, function(X) ginv(psd_proj(X)))
  
  #get vertex level test stat
  Wi_hat_list <- lapply(1:n, function(i) n*crossprod(D_hat[i,], var.here.invert[[i]]) %*% D_hat[i,])
  
  #get and return W 
  W_hat <- as.numeric(Reduce('+', Wi_hat_list))
  return(W_hat)
  
}