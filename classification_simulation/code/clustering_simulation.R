#--------------------------------------
#
#       Classification Simulations
#
#--------------------------------------

#load packages
library(mclust)

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
make_omni <- function(A,B,C){
  kronecker(H(1,2), A) + kronecker(H(2,2), B)
}
make_omni2 <- function(mats){
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
  
  n <- length(x)
  mc1 <- sum(abs(x - samp))
  mc2 <- sum(abs(ifelse(x == 1, 2, 1) - samp))
  
  return(min(c(mc1/n, mc2/n)))
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


#-----------------------------------------
#
#    Set up of Base Model
#
#-----------------------------------------

#set up blocks
B <- matrix(c(.25, .05, .05, .25), byrow = T, nrow = 2)
b_ase <- ase(B, 2)

#get latent positions
x1 <- b_ase[1,]; x2 <- b_ase[2,]

#set prior probabilities
pi <- .5

#get rotation (eigenvectors of Delta)
Delta <- pi * tcrossprod(x1) + (1 - pi)*tcrossprod(x2) 
R <- eigen(Delta)$vectors

#Apply rotation to x1 and x2
x1_til <- t(R) %*% x1 
x2_til <- t(R) %*% x2

#set base latent positions
L <- rbind(t(x1_til), t(x2_til))

#-----------------------------------------
#
#     Set up C values
#
#-----------------------------------------

#converges to ER with p = .3 (1,1) --> (2,0)
C <- function(t){
  diag(c(t + 1, -t + 1))
}

#-----------------------------------------
#
#    Set up Bias matrices 
#
#-----------------------------------------

S <- function(C){
  #bias matrices
  v1 <- c(1, C[1,1])
  v2 <- c(1, C[2,2])
  
  #get embeddings
  a1 <- ase(H1(v1), 1)[,1]
  a2 <- ase(H1(v2), 1)[,1]  
  
  #define scaling matrices
  S1 <- diag(c(a1[1], a2[1]))
  S2 <- diag(c(a1[2], a2[2]))
  
  #return results
  return(list(S1, S2))
}

#-----------------------------------------
#
#     Set up Community 
#     Detection Simulation
#
#-----------------------------------------

#set up base parameters 
#net_size <- c(25, 50, 75, 100, 150, 200) #network sizes 
net_size <- round(10^(seq(log(25, base = 10), log(250, base = 10), length.out = 10))) #network sizes 
net_size <- sort(c(net_size, 100))

t <- seq(0, 1, length.out = 11)[-11]
mc_runs <- 200 #number of iterations

#set up storage
here <- 1
df_list <- list()

#set seed 
set.seed(1985)

for(i in 1:length(net_size)){
  for(j in 1:length(t)){
    
    #iterate over mc_runs
    for(k in 1:mc_runs){
    
    #sample rows of X
    samp <- sample(1:2, net_size[i], replace = TRUE) 
    X <- L[samp, ]
    
    #Get bias scaling matrices 
    tmp <- S(C(t[j]))
    S1 <- tmp[[1]]
    S2 <- tmp[[2]]
    
    #set up P matrices
    P1 <- tcrossprod(X)
    P2 <- tcrossprod(X %*% C(t[j]), X)
    
      #sample A1 and A2 
      A1 <- sampP(P1)
      A2 <- sampP(P2)
      
      #----------------------
      #   ASE Clustering
      #----------------------
      
      #embedd individually 
      X_hat <- ase(A1, 2)
      Y_hat <- ase(A2, 2)
      
      #run GMM on rows
      groups1 <- Mclust(X_hat, G = 2, modelNames = "VVV")$classification
      groups2 <- Mclust(Y_hat, G = 2, modelNames = "VVV")$classification
      
      #get misclassification rates
      ase_mc1 <- get_mc(groups1, samp)
      ase_mc2 <- get_mc(groups2, samp)
      
      #get mahalanobis distances
      ase_graph1_dist <- get_mahalanobis(X_hat[samp == 1,], X_hat[samp == 2,])
      ase_graph2_dist <- get_mahalanobis(Y_hat[samp == 1,], Y_hat[samp == 2,])
      #----------------------
      #   Abar Clustering
      #----------------------
      
      #embedd A bar
      X_hat <- ase((A1 + A2)/2, 2)
      
      #run GMM on rows
      groups1 <- Mclust(X_hat, G = 2, modelNames = "VVV")$classification
      
      #get misclassification rates
      abar_mc1 <- abar_mc2 <- get_mc(groups1, samp)
      
      #get mahalanobis distances
      abar_graph1_dist <- abar_graph2_dist <- get_mahalanobis(X_hat[samp == 1,], X_hat[samp == 2,])
      
      #----------------------
      #   Omni Clustering
      #----------------------
      
      #construct Omni + embed
      Atil <- make_omni(A1, A2)
      Lhat <- ase(Atil, 2)
      
      #run GMM on rows
      groups1 <- Mclust(Lhat[1:net_size[i],], G = 2, modelNames = "VVV")$classification
      groups2 <- Mclust(Lhat[(net_size[i]+1):(2*net_size[i]), ], G = 2, modelNames = "VVV")$classification
      
      #get misclassification rates
      omni_mc1 <- get_mc(groups1, samp)
      omni_mc2 <- get_mc(groups2, samp)
      
      #get mahalanobis distances
      omni_graph1_dist <- get_mahalanobis(Lhat[1:net_size[i],][samp == 1,], Lhat[1:net_size[i],][samp == 2,])
      omni_graph2_dist <- get_mahalanobis(Lhat[(net_size[i]+1):(2*net_size[i]), ][samp == 1,], 
                                          Lhat[(net_size[i]+1):(2*net_size[i]), ][samp == 2,])
      
      #----------------------
      #   Omni bar Clustering
      #----------------------
     
      #get average
      Lmean <- .5*(Lhat[1:net_size[i],] + Lhat[(net_size[i]+1):(2*net_size[i]), ])
      
      #run GMM  on rows
      groups <- Mclust(Lmean, G = 2, modelNames = "VVV")$classification
      
      #get misclassification rates
      omnibar_mc1 <- omnibar_mc2 <- get_mc(groups, samp)
      
      #get mahalanobis distances
      omnibar_graph1_dist <- omnibar_graph2_dist <- get_mahalanobis(Lmean[samp == 1,], Lmean[samp == 2,])
      
      #----------------------
      #   Store result
      #----------------------
      df_list[[here]] <- matrix(c(k, net_size[i], t[j], "ASE", "Graph 1", ase_graph1_dist, ase_mc1,
                                  k, net_size[i], t[j], "ASE", "Graph 2", ase_graph2_dist, ase_mc2,
                                      
                                  k, net_size[i], t[j], "Omni", "Graph 1", omni_graph1_dist, omni_mc1,
                                  k, net_size[i], t[j], "Omni", "Graph 2", omni_graph2_dist, omni_mc2,
                                      
                                  k, net_size[i], t[j], "Abar", "Graph 1", abar_graph1_dist, abar_mc1,
                                  k, net_size[i], t[j], "Abar", "Graph 2", abar_graph2_dist, abar_mc2,
                                      
                                  k, net_size[i], t[j], "Omnibar", "Graph 1", omnibar_graph1_dist, omnibar_mc1,
                                  k, net_size[i], t[j], "Omnibar", "Graph 2", omnibar_graph2_dist, omnibar_mc2), 
                                              nrow = 8, ncol = 7, byrow = TRUE)
    
      #update counter 
      here <- here + 1
    }      
  }
  print(i)
}

#----------------------------------------
#
#       Data Cleaning + Storage
#
#----------------------------------------
library(dplyr)
df <- as.data.frame(do.call("rbind", df_list))
#recast variables
df[,1] <- as.numeric(as.character(df[,1]))
df[,2] <- as.numeric(as.character(df[,2]))
df[,3] <- as.numeric(as.character(df[,3]))
df[,6] <- as.numeric(as.character(df[,6]))
df[,7] <- as.numeric(as.character(df[,7]))

#name columns
colnames(df) <- c("sim_number","net_size", "t","Method", "Graph","dist","mc_rate")

#collapse over number of simulations
plotdf <- df %>% 
  group_by(net_size, t, Method, Graph) %>%
  summarize(MC_Rate = mean(mc_rate, na.rm = TRUE),
            mc_se = sd(mc_rate, na.rm = TRUE)/sqrt(n()),
            Dist = mean(dist, na.rm = TRUE),
            dist_se = sqrt(var(dist, na.rm = TRUE)/n()), 
            count = n())

#------------------------------
#  save data frame
#-----------------------------
setwd("~/Documents/Work/github/BJSE/classification_simulation/data")
write.csv(plotdf, "plotdf_gmm.csv")
