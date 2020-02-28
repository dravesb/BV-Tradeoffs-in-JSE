#------------------------------------------------------------------------------------------------
#            Simulation set up
#
#  - Define latent positions 
#  - Define simulation parameters
#------------------------------------------------------------------------------------------------

#-----------------------
#simulation paramters
#-----------------------

mc_runs <- 1000
net_size <- c(50, 100, 200)
t <- seq(0, 0.5, by = 0.05)
C <- function(t) diag(c(1+t, 1-t), nrow = 2)
alpha <- 0.05

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





