#-----------------------------------------
#
#    Simulation Functions
#
#-----------------------------------------

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
  kronecker(H(1,3), A) + kronecker(H(2,3), B) + kronecker(H(3,3), C)
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
pm <- function(A,d){
  E <- eigen(A)
  U <- normalize.cols(as.matrix(E$vectors[,1:d], ncol = d))
  S.inv <- diag(x = sign(E$values[1:d]), ncol = d, nrow = d)*diag(1/(sqrt(abs(E$values[1:d]))), nrow = d, ncol = d)
  U %*% S.inv
}

#expansion functions
library(expm) 
expansion <- function(C, theta, s){
  
  J <- matrix(c(0,1,-1,0), nrow = 2)
  tot <- matrix(0, nrow = 2, ncol = 2)
  
  for(l in 0:s){ # loop over sum variable
    for(k in 0:l){ #loop over other variable
      tot <- tot + (theta^l)/(factorial(k) * factorial(l-k)) * crossprod(J%^%k, C)%*%(J%^%(l-k))
    }
  }
  return(tot)
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

library(MCMCpack)
#Define three different C matrices
C1 <- diag(c(.75, .5))
C2 <- diag(c(.5, .75))
C3 <- diag(c(1, 0))

#-----------------------------------------
#
#     Set up S matrices 
#
#-----------------------------------------

#v vectors
v1 <- c(C1[1,1], C2[1,1], C3[1,1])
v2 <- c(C1[2,2], C2[2,2], C3[2,2])

#ase's
alpha1 <- ase(H1(v1), 1)[,1]
alpha2 <- ase(H1(v2), 1)[,1]

#S matrices
S1 <- diag(c(alpha1[1], alpha2[1]))
S2 <- diag(c(alpha1[2], alpha2[2]))
S3 <- diag(c(alpha1[3], alpha2[3]))

#Bias 
B1 <- L %*% (S1 - sqrt(C1))
B2 <- L %*% (S2 - sqrt(C2))
B3 <- L %*% (S3 - sqrt(C3))

#Make Bias Dataframes
bias_dat <- cbind(rbind(B1, B2, B3),
                  rep(c("Graph 1", "Graph 2", "Graph 3"), each = 2), 
                  rep(c("Community 1", "Community 2"), 3))
bias_dat <- data.frame(x = as.numeric(bias_dat[,1]),
                       y = as.numeric(bias_dat[,2]),
                       Graph = bias_dat[,3],
                       Community = bias_dat[,4])


#-----------------------------------------
#
#     Do simulation
#
#-----------------------------------------
library(dplyr)
net_size <- c(250,500,750,1000)
mc_runs <- 25
expan.size <- 0


#storage
df <- matrix(NA, ncol = 14, nrow = 3*mc_runs*sum(net_size))
start <- stop <-  0
colnames(df) <-c("Network.Size","Graph","Community","Balance",
                 "Theta", "Det",
                 "Xhat","Yhat",
                 "XC1","XC2",
                 "XS1","XS2",
                 "RX", "RY")

set.seed(1985)


for(i in 1:length(net_size)){
  #iterate over group sizes
  #com.balance <- round(seq(net_size[i]*(pi - 1/sqrt(net_size[i])), 
  #                         net_size[i]*(pi + 1/sqrt(net_size[i])),
  #                         length.out = mc_runs))

  #even group sizes
  com.balance <- net_size[i]/2
    
  for(j in 1:mc_runs){
    
    #-----------------------------------------------------------
    #                 Set up X matrices
    #-----------------------------------------------------------
    
    #sample with replacement from L
    ind <- c(rep(1, com.balance), rep(2, net_size[i] - com.balance))
    
    #set up latent positions and right singular vectors
    X <- L[ind,]
    
    tmp <- svd(X) 
    D <- diag(sign(diag(tmp$v)))
    V <- tmp$v 
    Sigma <- diag(tmp$d)
    U <- tmp$u 
    
    #set up Xtil 
    Xtil <- rbind(U %*% Sigma %*% sqrt(C1), 
                  U %*% Sigma %*% sqrt(C2), 
                  U %*% Sigma %*% sqrt(C3))
    
    #calculate theta
    theta.here <- asin((D %*% V)[2,1])
    
    #-----------------------------------------------------------
    #                 Set up P matrices
    #-----------------------------------------------------------
    
    #set up P matrices
    P1_constant <- tcrossprod(X %*% V %*% (expansion(C1, theta.here, expan.size)), X %*% V)
    P1_linear <- tcrossprod(X %*% V %*% (expansion(C1, theta.here, expan.size+1)), X %*% V)
    P1V <- tcrossprod(X %*% C1, X)
    
    P2_constant<- tcrossprod(X %*% V %*% (expansion(C2, theta.here, expan.size)), X %*% V)
    P2_linear <- tcrossprod(X %*% V %*% (expansion(C2, theta.here, expan.size+1)), X %*% V)
    P2V <- tcrossprod(X %*% C2, X)
    
    P3_constant <- tcrossprod(X %*% V %*% (expansion(C3, theta.here, expan.size)), X %*% V)
    P3_linear <- tcrossprod(X %*% V %*% (expansion(C3, theta.here, expan.size+1)), X %*% V)
    P3V <- tcrossprod(X %*% C3, X)
    
    #set up embedding matrices
    P_constant <- make_omni(P1_constant, P2_constant, P3_constant)
    P_linear <- make_omni(P1_linear, P2_linear, P3_linear)
    PV <- make_omni(P1V, P2V, P3V)
    
    #embeddings
    Z0 <- ase(P_constant, 2)
    Z1 <- ase(P_linear, 2)
    Zinf <- ase(PV, 2)
    
    #power method term
    M0 <- pm(P_constant, 2)
    #Minf <- pm(PV, 2)
    
    #align
    Z0 <- procrustes(Z0, Xtil)$X.new
    Z1 <- procrustes(Z1, Xtil)$X.new
    Zinf <- procrustes(Zinf, Xtil)$X.new
    M0 <- procrustes(M0, Xtil)$X.new
    #Minf <- procrustes(Minf, Xtil)$X.new
    
    
    #-----------------------------------------------------------
    #                 Sample & Embedd A matrices
    #-----------------------------------------------------------
    
    #sample A matrices
    A1 <- sampP(P1V)
    A2 <- sampP(P2V)
    A3 <- sampP(P3V)
    
    #Embedd
    Atil <- make_omni(A1, A2, A3)
    Lhat <- ase(Atil, 2)
    
    #align 
    Lhat <- procrustes(Lhat, Xtil)$X.new
    
    #-----------------------------------------------------------
    #                 Set up difference matrices
    #-----------------------------------------------------------
    
    #get scaled matrices
    XC <- rbind(X %*% C1, X %*% C2, X %*% C3) # Bias + Residual 
    XS <- rbind(X %*% S1, X %*% S2, X %*% S3)
    
    #get residual
    R <- Lhat - XS
    
    #-----------------------------------------------------------
    #                 Store Data
    #-----------------------------------------------------------
    
    start <- stop + 1
    stop <- stop + 3*net_size[i]
    
    #store
    df[start:stop,] <- as.matrix(data.frame(Network.Size = rep(net_size[i],3*net_size[i]), #network size
                                            Graph = c(rep("Graph 1", net_size[i]),
                                                      rep("Graph 2", net_size[i]),
                                                      rep("Graph 3", net_size[i])), #which graph
                                            Community = rep(paste("Community", ind), 3), #community assignment
                                            Balance = rep(com.balance[j]/net_size[i], 3*net_size[i]), #MC replicant number
                                            Theta = rep(theta.here, 3 * net_size[i]), #theta,
                                            Det = rep(det(V), 3 * net_size[i]), #det,
                                            Xhat = Lhat[,1], #Lhat
                                            Yhat = Lhat[,2], 
                                            XC1 = XC[,1], #XS
                                            XC2 = XC[,2],
                                            XS1 = XS[,1], #XC
                                            XS2 = XS[,2],
                                            RX = R[,1], #Lhat - XS
                                            RY = R[,2]
    ))
    print(j)
  }
}

plotdf <- data.frame(Network.Size = as.numeric(df[,1]), #network size
                     Graph = df[,2], #which graph
                     Community = df[,3], #community assignment
                     Balance = paste(100*as.numeric(df[,4]), "%", sep = ""), #MC replicant number
                     Theta =  as.numeric(df[,5]), #theta,
                     Det =  as.numeric(df[,6]), #det,
                     Xhat =  as.numeric(df[,7]), #Lhat
                     Yhat =  as.numeric(df[,8]), 
                     XC1 =  as.numeric(df[,9]), # XC
                     XC2 =  as.numeric(df[,10]),
                     XS1 =  as.numeric(df[,11]),# XS 
                     XS2 =  as.numeric(df[,12]),
                     RX =  as.numeric(df[,13]), #R 
                     RY =  as.numeric(df[,14])
)


write.csv(plotdf, "~/Documents/Work/github/BJSE/bias_simulation/data/plotting_data.csv")

