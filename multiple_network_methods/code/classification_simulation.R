#--------------------------------------
#                                     
#       Multiple Graph Embedding      
#         Clustering Simulation       
#                                     
#--------------------------------------

#-----------------------------
#     Preliminaries
#-----------------------------

#load up necessary packages
library(ggplot2)
library(dplyr)
library(reshape2)
library(matlib)
library(irlba)

#source nesseary functions
source("~/Desktop/paper_figures/multiple_network_methods/code/je_classification.R")
source("~/Desktop/paper_figures/multiple_network_methods/code/mase_classification.R")
source("~/Desktop/paper_figures/multiple_network_methods/code/omni_classification.R")
source("~/Desktop/paper_figures/multiple_network_methods/code/mrdpg_classification.R")

#-----------------------------
#     Set up Model 
#-----------------------------

#set up blocks
B <- matrix(c(.3, .1, .1, 
              .1, .25, .15,
              .1, .15, .25), 
            byrow = T, nrow = 3)

b_ase <- ase(B, 3)

#get latent positions
x1 <- b_ase[1,]; x2 <- b_ase[2,]; x3 <- b_ase[3,]

#set prior probabilities
prior_probs <- c(1/3, 1/3, 1/3)

#get rotation (eigenvectors of Delta)
Delta <- prior_probs[1] * tcrossprod(x1) + 
         prior_probs[2] * tcrossprod(x2) + 
         prior_probs[3] * tcrossprod(x3) 
R <- eigen(Delta)$vectors

#Apply rotation to x1, x2, x3
x1_til <- t(R) %*% x1 
x2_til <- t(R) %*% x2
x3_til <- t(R) %*% x3

#set base latent positions
L <- t(cbind(x1_til, x2_til, x3_til))

#-----------------------------
#     Set up different C 
#        matrices
#-----------------------------

C1 <- function(t) diag(c(1, 1, 1-t)) #Reduces (C2, C3) into a single group
C2 <- function(t) diag(c(1, 1-t, 1)) #(C2,C3) make a 3 group SBM, between group and C1 prob. the same
C3 <- function(t) diag(c(1, 1-t, 1-t)) #Erdos-Reyni

#-----------------------------
#     Simulation parameters 
#-----------------------------

#distance from iid
k <- 20
t <- seq(0, 1, length.out =  k)

#number of mc iterations
mc_runs <- 500

#set network size
#n <- round(10^(seq(log(25, base = 10), log(250, base = 10), length.out = 10))) #network sizes 
n <- 100

#define storage
names <- c("network_size", "t","iter.no", "Omni", "MASE", "JE", "MRDPG")
df <- matrix(NA, ncol = length(names), nrow = k * mc_runs * length(n))
colnames(df) <- names
here <- 1 

#-----------------------------
#     Do simulation
#-----------------------------

for(l in 1:length(n)){#loop over n
  for(i in 1:k){#loop over k 
    for(j in 1:mc_runs){#monte carlo iterations
      
      #sample group assignments
      samp <- sample(1:3, n[l], replace = TRUE)
      X <- L[samp,]
      
      #get different P matrices
      P1 <- tcrossprod(X)
      P2 <- X %*% tcrossprod(C1(t[i]), X)
      P3 <- X %*% tcrossprod(C2(t[i]), X)
      P4 <- X %*% tcrossprod(C3(t[i]), X)
      
      #sample adjaceny matrices
      A1 <- sampP(P1)
      A2 <- sampP(P2)
      A3 <- sampP(P3)
      A4 <- sampP(P4)
  
      #make adjaceny lists
      adj_mats <- list()
      adj_mats[[1]] <- A1
      adj_mats[[2]] <- A2
      adj_mats[[3]] <- A3
      adj_mats[[4]] <- A4
      
      #get classifications from each methods
      omni_estimates <- omni_classes(adj_mats, 3, 3)
      mrdpg_estimates <- mrdpg_classes(adj_mats, 3, 3)
      je_estimates <- je_classes(adj_mats, 3, 3)
      mase_estimates <- mase_classes(adj_mats, 3, 3)
      
      #get MC rates
      omni_mc <- get_mc(omni_estimates, samp)
      mrdpg_mc <- get_mc(mrdpg_estimates, samp)
      je_mc <- get_mc(je_estimates, samp)
      mase_mc <- get_mc(mase_estimates, samp)    
      
      #store data
      df[here, ] <- c(n[l], t[i], j, omni_mc, mase_mc, je_mc, mrdpg_mc)
      
      #update counter 
      here <- here + 1
      
    }
    print(k)
  }
}

#-----------------------------
#     write out simulation
#-----------------------------

write.csv(df, "~/Desktop/paper_figures/multiple_network_methods/data/data.csv")






