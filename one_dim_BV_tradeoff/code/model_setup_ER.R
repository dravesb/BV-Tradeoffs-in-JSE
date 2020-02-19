#------------------------------------------------------------------------------------------------
#            Simulation set up
#
#  - Define latent positions 
#  - Define simulation parameters
#------------------------------------------------------------------------------------------------

#latent position matrix is just sqrt(p)1_n
p <- .5
L <- sqrt(p)

#simulation paramters
mc_runs <- 1000
net_size <- 100 
C <- seq(0, sqrt(2), length.out = 26)[-1]



