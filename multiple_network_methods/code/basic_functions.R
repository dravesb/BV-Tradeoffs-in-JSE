#--------------------------------------
#
#       Basic Function
#
#--------------------------------------

#sampling functions
getP <- function(X) tcrossprod(X)
sampBern <- function(p) rbinom(1,1,p) 
sampP <- function(P){
  A <- apply(P, c(1,2), sampBern)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) <- 0
  A
} 

#estimating functions
H <- function(g,m){
  ones <- rep(1, m)
  e <- diag(m)[,g]
  .5 * (tcrossprod(ones,e) + tcrossprod(e,ones))
}
H1 <- function(x){
  ones <- rep(1, length(x))
  .5 * (tcrossprod(ones,x) + tcrossprod(x,ones))
}
make_omni <- function(mats){
  m <- length(mats)
  Reduce("+", lapply(1:m, function(x) kronecker(H(x,m), mats[[x]])))
}
norm2 <- function(u){
  sqrt(sum(u^2))
} 
normalize.cols<- function(A){
  norm.vec <- function(u) u/norm2(u) #define vector normalization func.
  if(ncol(A) == 1) return(norm.vec(A[,1]))
  apply(A, 2, norm.vec) # vectorize 
} 
ase <- function(A,d){
  E <- eigen(A)
  U <- normalize.cols(as.matrix(E$vectors[,1:d], ncol = d))
  S <- diag(x = sign(E$values[1:d]), ncol = d, nrow = d)*diag(sqrt(abs(E$values[1:d])), nrow = d, ncol = d)
  U %*% S
}

#clustering functions
get_mc <- function(x, true){
  
  #get some helpful parameters
  n <- length(x)
  
  #unchanged 123
  mc1 <- sum((x - true) != 0)
  
  #swap 213
  swap <- x
  swap[x == 1] <- 2; swap[x == 2] <- 1; swap[x == 3] <- 3
  mc2 <- sum((swap - true) != 0)
  
  #swap 321
  swap <- x
  swap[x == 1] <- 3; swap[x == 2] <- 2; swap[x == 3] <- 1
  mc3 <- sum((swap - true) != 0)
  
  #swap 132
  swap <- x
  swap[x == 1] <- 1; swap[x == 2] <- 3; swap[x == 3] <- 2
  mc4 <- sum((swap - true) != 0)
  
  #swap 231
  swap <- x
  swap[x == 1] <- 2; swap[x == 2] <- 3; swap[x == 3] <- 1
  mc5 <- sum((swap - true) != 0)
  
  #swap 312
  swap <- x
  swap[x == 1] <- 3; swap[x == 2] <- 1; swap[x == 3] <- 2
  mc6 <- sum((swap - true) != 0)
  
  return(min(c(mc1, mc2, mc3, mc4, mc5, mc6)/n))
}

#mahalanobis distances
get_mahalanobis <- function(X, Y){
  #get sample sizes
  nx <- nrow(X)
  ny <- nrow(Y)
  
  #get covariances
  CX <- cov(X)
  CY <- cov(Y)
  
  #pool covariances
  Cpool <- ((nx - 1)*CX + (ny - 1)*CY)/(nx + ny -2)
  
  #get column means
  xbar <- colMeans(X)
  ybar <- colMeans(Y)
  
  #calcualte R^{-T}(x - y)
  Rc <- chol(Cpool)
  fact <- forwardsolve(t(Rc), xbar - ybar)
  
  #return distance
  norm2(fact)
  
  #return average distance
  #return(norm2(xbar - ybar))
  
  #return average Mah distance
  #return(0.5 *(mahalanobis(xbar - ybar, center = FALSE, CX, inverted = TRUE) + mahalanobis(xbar - ybar, center = FALSE, CY, inverted = TRUE)))
  
  #take cholesky decomposition
  #Rx <- chol(CX)
  #Ry <- chol(CY)
  
  
  #standarize variances
  #xstan <- backsolve(Rx, xbar)
  #ystan <- backsolve(Ry, ybar)
  
  #get return distance
  #norm2(xstan - ystan)
  
}

