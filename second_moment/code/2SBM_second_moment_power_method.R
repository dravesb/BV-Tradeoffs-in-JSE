#--------------------------------------
#
#       Second Moment Simulations
#         - Power Method Focused
#           
#--------------------------------------

#-----------------------------
#     Preliminaries
#-----------------------------

#load up necessary packages
library(ggplot2)
library(dplyr)
library(reshape2)

#source nesseary functions
source("~/Documents/Work/github/BJSE/second_moment/code/basic_functions.R")

#-----------------------------
#     Set up Base Model
#-----------------------------

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

#-----------------------------
#     Set up C values
#-----------------------------

#converges to ER with p = .3 (1,1) --> (2,0)
C <- function(t){
  diag(c(t + 1, -t + 1))
}

#-----------------------------
#    Set up Bias matrices 
#-----------------------------

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

#-----------------------------
#    Simulation Parameters
#-----------------------------

#set up base parameters 
net_size <- 200
t <- c(0,.25,.5,.75,1)#seq(0, 1, length.out = 11)[-11]
mc_runs <- 100 #number of iterations


#set up storage
here <- 1
df <- matrix(NA, nrow = 2 * mc_runs * sum(net_size) * length(t), ncol = 9)
colnames(df) <- c("Network Size", "iter", "t", "Graph", "Community", "R1X", "R1Y", "R2X", "R2Y")
here <- 1 

#set seed 
set.seed(1985)


for(i in 1:length(net_size)){#iterate over network size
  for(j in 1:length(t)){#iterate over t values
    for(k in 1:mc_runs){#number of MC simulations
      
      #sample group assignments
      samp <- sample(1:2, size = net_size[i], replace = TRUE)
      
      #set up X and P1 and P2
      X <- L[samp,]
      P1 <- tcrossprod(X)
      P2 <- tcrossprod(X, X %*% C(t[j]))
      
      #sample A1 and A2
      A1 <- sampP(P1)
      A2 <- sampP(P2)
      
      #make Omni and embedd for Lhat and Z
      Lhat <- ase(make_omni(list(A1, A2)), 2)
      Z <- ase(make_omni(list(P1, P2)), 2)
      
      #make Ls
      scaling <- S(C(t[j]))
      Ls <- rbind(X %*% scaling[[1]], X %*% scaling[[2]])
      
      #rotate both to Ls
      Lhat.rot <- rotate(Ls, Lhat)
      Z.rot <- rotate(Ls, Z)
      
      #Get Main residual Lhat - Ls  
      R0 <- Lhat.rot - Ls
      
      #power method term
      M <- pm(make_omni(list(P1, P2)), 2)
      PM1 <- make_omni(list(A1 - P1, A2 - P2)) %*% M 
      
      #Ls power method term 
      Sp.inv <- diag(1/eigen(make_omni(list(P1, P2)), only.values = TRUE)$values[1:2])
      PM2 <- make_omni(list(A1 - P1, A2 - P2)) %*% Ls %*% Sp.inv
      
      #rotate terms to Ls
      PM1.rot <- rotate(Ls, PM1)
      PM2.rot <- rotate(Ls, PM2)
      
      #get residual differences
      R1 <- sqrt(net_size[i]) * (R0 - PM1.rot)
      R2 <- sqrt(net_size[i]) * (R0 - PM2.rot)
      
      
      #store residuals in dataframe
      update <- here + 2*net_size[i] - 1
      df[here:update,] <- cbind(rep(net_size[i], 2*net_size[i]), # network size
                                rep(k, 2*net_size[i]), # MC run
                                rep(t[j], 2*net_size[i]), # t
                                c(rep("Graph 1", net_size[i]), rep("Graph 2", net_size[i])), # Graph
                                rep(samp, 2), # Community assignmentl
                                R1, #(Lhat - Ls) - (A - P)M residual
                                R2 #(Lhat - Ls) - (A - P)Ls S^{-1} residual
      )
      here <- update + 1
    }
    print(c(i,j,k))
  }
}


#-----------------------------
#     write out data
#-----------------------------

write.csv(df, "~/Documents/Work/github/BJSE/second_moment/data/power_method_data.csv")
