#------------------------------------------------------------------------------------------------
#             Basic Functions
#
# sampP: samples adjacency matrix from P 
# make_omni: takes list of adjacency matrix and makes Omnibus matrix
# ase: Adjacency spectral embedding
# get_S_list: Takes a list of C matrices and maps them to S matrices
#------------------------------------------------------------------------------------------------

#sampling functions
sampP <- function(P) {
  #set up matrix
  n = ncol(P)
  A = matrix(0, n, n)
  
  #sample from A 
  A[upper.tri(A)] <- 1*(runif(sum(upper.tri(P))) < P[upper.tri(P)])
  A = A + t(A)
  return(A)
}

#Make Omni function
make_omni <- function(mats){
  #H(x) = (1x^T + x1^T)/2
  H <- function(g,m){
    ones <- rep(1, m)
    e <- diag(m)[,g]
    .5 * (tcrossprod(ones,e) + tcrossprod(e,ones))
  }
  
  #sum up each kronecker 
  m <- length(mats)
  Reduce("+", lapply(1:m, function(x) kronecker(H(x,m), mats[[x]])))
}

#ASE Embedding
ase <- function(A,d){
  #function to normalize columns of A
  norm2 <- function(u){
    sqrt(sum(u^2))
  } 
  normalize.cols<- function(A){
    norm.vec <- function(u) u/norm2(u) #define vector normalization func.
    if(ncol(A) == 1) return(norm.vec(A[,1]))
    apply(A, 2, norm.vec) # vectorize 
  } 
  
  #construct ASE 
  E <- eigen(A)
  U <- normalize.cols(as.matrix(E$vectors[,1:d], ncol = d))
  S_sqrt <- diag(x = sign(E$values[1:d]), ncol = d, nrow = d)*diag(sqrt(abs(E$values[1:d])), nrow = d, ncol = d)
  U %*% S_sqrt
}

#vector mse functions (assumed x centered)
get_mse <- function(x){
  sum(x^2)
}

#Maps C matrices to S matrices
get_Slist <- function(C_list){
  #fetch m
  m <- length(C_list)
  d <- nrow(C_list[[1]])
  
  #V (d x m)  matrix with v1, v2, ..., vd in the __rows__ of V
  V <- Reduce('cbind', lapply(C_list, diag))
  
  #Alpha (m x d)  matrix with alpha^1, alpha^2, ..., alpha^d in the __columns__ of Alpha
  Alpha <- apply(V, 1, function(v){
    H.here <- 0.5* (tcrossprod(rep(1, m), v) + tcrossprod(v, rep(1, m)))
    eigen.system <- eigen(H.here)
    u <- eigen.system$vectors[, 1]
    s <- eigen.system$values[1]
    s_half <- sign(s) * sqrt(abs(s))
    alpha <- abs(u * s_half) #chosen to have positive values 
    return(alpha)
  })
  
  #Set up S matrices by taking the diagnol of the __rows__ of Alpha
  S_list <- list()
  for(g in 1:m) S_list[[g]] <- diag(Alpha[g, ], nrow = d, ncol = d)
  
  return(S_list)
  

}


psd_proj <- function(M){
  #decomp 
  eigen_decomp <- eigen(M)
  U <- eigen_decomp$vectors
  S <- eigen_decomp$values
  
  #threshold eigenvalues
  S.thres <- diag(sapply(S, function(x) max(c(0,x))))
  
  #return
  return(tcrossprod(U %*% S.thres, U))
  
}


