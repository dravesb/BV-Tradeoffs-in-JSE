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
t <- seq(0, 1, length.out = 5)

#-----------------------
#set up latent positions
#-----------------------

# #set up blocks
B <- matrix(c(.3, .1, .1, 
                .1, .25, .15,
                .1, .15, .25), 
              byrow = T, nrow = 3)


#B <- matrix(c(0.42, 0.30, 0.08, 
#              0.30, 0.42, 0.20,
#              0.08, 0.20, 0.42), 
#            byrow = T, nrow = 3)

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
 



