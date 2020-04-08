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

gs.dim.select <- function(X, k=NULL, edge.attr=NULL, n = 3, threshold = FALSE, plot = FALSE) {
  
  if (class(X) != "numeric" || class(X) != "matrix" || class(X) != "array" || class(X) != "igraph") {
    stop("You have passed neither a 2-D matrix/array, an igraph object, nor a 1-D numeric vector.")
  }
  if (class(X) != 'numeric') {
    if (class(X) == "graph") {
      X <- gs.as_adj(X, edge.attr=edge.attr)
    }
    if (length(dim.X) > 2) {
      stop("You have input an array with more than 2 dimensions.")
    }
    if (is.null(k)) {
      k <- floor(log2(min(dim(X))))
    }
    if (dim.X[1] > 1 && dim.X[2] > 1) {
      d <- gs.embed(X, k)$D
    } else {
      d <- X
    }
  }
  tryCatch({
    d <- as.numeric(d)
  }, warning=function(w) {stop("The singular values have invalid entries and cannot be cast to numeric.")})
  
  
  if (sum("numeric" == apply(as.array(d),1,class)) != length(d)) { stop("Input object 'sdev' contains non-numeric values.") }
  d <- sort(d, decreasing=TRUE)
  if (is.numeric(threshold)) {
    d <- d[d > threshold]
    p <- length(d)
    if (p == 0) {
      stop(sprintf("The singular values do not have any elements larger than the threshold value, %f. The maximum singular value is %.3f.",
                   threshold, max(d)))
    }
  } else {
    p <- length(d)
  }
  
  if (class(n) != "numeric") {
    stop("Input object 'n' is not a numeric value.")
  } else if (length(n) != 1) {
    stop("Input object 'n' is not a numeric value.")
  } else {
    if (n <= 0 || n >= length(d)) {
      stop("Input object 'n' must be in the appropriate interval. n must be between 0 and %d, but n is %d.", length(d), n)
    }
  }
  
  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
      (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  
  q <- which.max(lq)
  if (n > 1 && q < (p-1)) {
    q <- c(q, q + gs.dim.select(d[(q+1):p], n=n-1, plot=FALSE)$elbow)
  }
  
  if (!is.logical(plot)) { stop("Input object 'plot' is not a logical.") }
  
  out <- list(value=d[q], elbow=q)
  if (plot) {
    sv.obj <- data.frame(Dimension=1:p, Value=d)
    elbow.obj <- data.frame(Dimension=q, Value=d[q])
    plot <- ggplot(sv.obj, aes(x=Dimension, y=Value)) +
      geom_line(color='blue') +
      geom_point(data=elbow.obj, aes(x=Dimension, y=Value), color='red') +
      xlab("Dimension") +
      ylab("Singular Value") +
      theme_bw()
    out$plot <- plot
  }
  
  return(out)
}

