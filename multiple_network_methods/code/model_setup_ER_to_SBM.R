#------------------------------------------------------------------------------------------------
#            Simulation set up
#
#  - Define latent positions 
#  - Define simulation parameters
#------------------------------------------------------------------------------------------------

#-----------------------
#simulation paramters
#-----------------------

mc_runs <- 500
net_size <- 100
t <- seq(0, 1, length.out = 26)[-26]
C <- function(t) diag(c(1+t, 1-t), nrow = 2)

#-----------------------
#set up latent positions
#-----------------------

#blocks
B <- matrix(c(.25, .05, .05, .25), byrow = T, nrow = 2)
b_ase <- ase(B, 2)

#get latent positions
x1 <- b_ase[1,]; x2 <- b_ase[2,]

#set prior probabilities
pi <- .5
probs <- c(.5,.5)

#get rotation (eigenvectors of Delta)
Delta <- pi * tcrossprod(x1) + (1 - pi)*tcrossprod(x2) 
R <- eigen(Delta)$vectors

#Apply rotation to x1 and x2
x1_til <- t(R) %*% x1 
x2_til <- t(R) %*% x2

#set base latent positions
L <- rbind(t(x1_til), t(x2_til))

# #set up blocks
# B <- matrix(c(.3, .1, .1, 
#               .1, .25, .15,
#               .1, .15, .25), 
#             byrow = T, nrow = 3)
# 
# b_ase <- ase(B, 3)
# 
# #get latent positions
# x1 <- b_ase[1,]; x2 <- b_ase[2,]; x3 <- b_ase[3,]
# 
# #set prior probabilities
# prior_probs <- c(1/3, 1/3, 1/3)
# 
# #get rotation (eigenvectors of Delta)
# Delta <- prior_probs[1] * tcrossprod(x1) + 
#          prior_probs[2] * tcrossprod(x2) + 
#          prior_probs[3] * tcrossprod(x3) 
# R <- eigen(Delta)$vectors
# 
# #Apply rotation to x1, x2, x3
# x1_til <- t(R) %*% x1 
# x2_til <- t(R) %*% x2
# x3_til <- t(R) %*% x3
# 
# #set base latent positions
# L <- t(cbind(x1_til, x2_til, x3_til))
# 
# #-----------------------------
# #     Set up different C 
# #        matrices
# #-----------------------------
# 
# C1 <- function(t) diag(c(1, 1, 1-t)) #Reduces (C2, C3) into a single group
# C2 <- function(t) diag(c(1, 1-t, 1)) #(C2,C3) make a 3 group SBM, between group and C1 prob. the same
# C3 <- function(t) diag(c(1, 1-t, 1-t)) #Erdos-Reyni
# 
# #-----------------------------
# #     Simulation parameters 
# #-----------------------------
# 
# #distance from iid
# k <- 20
# t <- seq(0, 1, length.out =  k)
# 
# #number of mc iterations
# mc_runs <- 500
# 
# #set network size
# #n <- round(10^(seq(log(25, base = 10), log(250, base = 10), length.out = 10))) #network sizes 
# #n <- 100
# 



