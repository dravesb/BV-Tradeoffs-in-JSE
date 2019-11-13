plotdf <- read.csv("~/Documents/Work/github/BJSE/bias_simulation/data/plotting_data.csv", header = TRUE)[,-1]
#------------------------------------------------
#
#         Cirlce Funtions
#
#------------------------------------------------
library(dplyr)

circleFun <- function(center = c(0,0), r = 1, npoints = 100){
  tt <- seq(0,2*3.14,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

nsizes <- unique(plotdf$Network.Size)
nsizes <- unique(plotdf$Network.Size)

f1 <- function(x) log(3*x)/sqrt(x)
f2 <- function(x) sqrt(log(3*x)/x)

circ_dat1 <- lapply(nsizes, function(x) circleFun(r = f1(x))) %>% bind_rows()
circ_dat1 <- circ_dat1 %>% mutate(Network.Size = rep(nsizes, each = 100)) 

circ_dat2 <- lapply(nsizes, function(x) circleFun(r = f2(x))) %>% bind_rows()
circ_dat2 <- circ_dat2 %>% mutate(Network.Size = rep(nsizes, each = 100)) 

#------------------------------------------------
#
#         Confidence Covariances
#
#------------------------------------------------

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

#Define three different C matrices
C1 <- diag(c(.75, .5))
C2 <- diag(c(.5, .75))
C3 <- diag(c(1, 0))
C_list <- list(C1, C2, C3)

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
S_list <- list(S1, S2, S3)

#calculate sigma tildes
#Variance matrices
sigma_2 <- function(y, C, k){
  as.numeric(crossprod(y, C %*% L[k,]) - crossprod(y, C %*% L[k,])^2)
}
Sigma_tilde <- function(y, g){
  
  #set up preliminaries
  K <- nrow(L)
  d <- ncol(L)
  
  #get sum
  tot <- matrix(0,nrow = d, ncol = d)
  for(k in 1:K){
    tot <- tot + probs[k] * sigma_2(y, C_list[[g]], k) * tcrossprod(L[k,])
    
  }
  return(tot)
  
}
S2D_inv <- solve(.5 * crossprod(L) %*%  Reduce("+", lapply(S_list, function(x) x^2)))
Sigma <- function(y,g){
  
  #preliminary values
  d <- ncol(L)
  m <- length(C_list)
  
  #first summand
  tot1 <- (S_list[[g]] + Reduce("+",S_list)) %*% Sigma_tilde(y, g) %*% (S_list[[g]] + Reduce("+",S_list))
  
  #second summard
  tot2 <- matrix(0, nrow = d, ncol = d)
  for(k in (1:m)[-g]){
    tot2 <- tot2 + S_list[[k]] %*% Sigma_tilde(y, k) %*% S_list[[k]]
  }
  
  #return covariance matrix
  return(.25 * S2D_inv %*% (tot1 + tot2) %*% S2D_inv)
}

#create confidence ellipses
ellips <- function(center = c(0,0), c = qnorm(.975), Sigma, npoints = 1000-1){
  angles <- (0:npoints) * 2 * 3.14/npoints
  unit.circle <- cbind(cos(angles), sin(angles))
  Q <- chol(Sigma, pivot = TRUE)
  order <- order(attr(Q, "pivot"))
  df <- t(center + c * t(unit.circle %*% Q[, order]))
  colnames(df) <- c("x", "y")
  return(df)
}

ell_df <- bind_rows(
  list(
    as.data.frame(do.call(rbind, lapply(nsizes, function(x) ellips(as.numeric(S1 %*% L[1,]), Sigma = Sigma(L[1,], 1)/x)))) %>% 
    mutate(Network.Size = rep(nsizes, each = 1000),Graph = rep("Graph 1", 3000),Community = rep("Community 1", 3000)),
    as.data.frame(do.call(rbind, lapply(nsizes, function(x) ellips(as.numeric(S1 %*% L[2,]), Sigma = Sigma(L[2,], 1)/x)))) %>% 
      mutate(Network.Size = rep(nsizes, each = 1000),Graph = rep("Graph 1", 3000),Community = rep("Community 2", 3000)),
    
    as.data.frame(do.call(rbind, lapply(nsizes, function(x) ellips(as.numeric(S2 %*% L[1,]), Sigma = Sigma(L[1,], 2)/x)))) %>% 
      mutate(Network.Size = rep(nsizes, each = 1000),Graph = rep("Graph 2", 3000),Community = rep("Community 1", 3000)),
    as.data.frame(do.call(rbind, lapply(nsizes, function(x) ellips(as.numeric(S2 %*% L[2,]), Sigma = Sigma(L[2,], 2)/x)))) %>% 
      mutate(Network.Size = rep(nsizes, each = 1000),Graph = rep("Graph 2", 3000),Community = rep("Community 2", 3000)),
    
    as.data.frame(do.call(rbind, lapply(nsizes, function(x) ellips(as.numeric(S3 %*% L[1,]), Sigma = Sigma(L[1,], 3)/x)))) %>% 
      mutate(Network.Size = rep(nsizes, each = 1000),Graph = rep("Graph 3", 3000),Community = rep("Community 1", 3000)),
    as.data.frame(do.call(rbind, lapply(nsizes, function(x) ellips(as.numeric(S3 %*% L[2,]), Sigma = Sigma(L[2,], 3)/x)))) %>% 
      mutate(Network.Size = rep(nsizes, each = 1000),Graph = rep("Graph 3", 3000),Community = rep("Community 2", 3000))
))

colnames(ell_df)[c(1:2)] <- c("X", "Y")


#-------------------------------------------------
#
#        Bias Bounds
#
#-------------------------------------------------
library(ggplot2)

ggplot() +
  geom_point(aes(Xhat, Yhat, col = Community), plotdf, alpha = 0) +
  geom_point(aes(XS1, XS2),  shape = 8, size = 1, plotdf) +
  geom_point(aes(XC1, XC2),  shape = 3, size = 1, plotdf) +
  facet_grid(rows = vars(Graph),
             cols = vars(Network.Size))+ 
  geom_path(aes(X, Y, col = Community), ell_df)+
  scale_color_manual(values=c("red", "blue"))+
  theme_bw()+
  labs(x = "X", y = "Y",
       col = "Community")

ggsave("analytic_bias.jpeg", 
       width = 7, height = 6, units = "in", 
       device = "jpeg", 
       path = "~/Documents/Work/github/BJSE/bias_simulation/figures/")

ggplot() +
  geom_point(aes(Xhat, Yhat, col = Community), plotdf,
             size = .1,
             alpha = .05
             )+
  geom_point(aes(XS1, XS2),  shape = 4, size = 1, plotdf) +
  geom_point(aes(XC1, XC2),  shape = 3, size = 1, plotdf) +
  geom_path(aes(X, Y, col = Community), ell_df, lty = 2)+
  facet_grid(rows = vars(Graph),
             cols = vars(Network.Size))+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_manual(values=c("red", "blue"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1)))+
  labs(x = "X", y = "Y",
       col = "Community")

ggsave("observed_bias.jpeg", 
       width = 6, height = 6, units = "in", 
       device = "jpeg", 
       path = "~/Documents/Work/github/BJSE/bias_simulation/figures/")

#-------------------------------------------------
#
#        Residual Bounds
#
#-------------------------------------------------

ggplot() +
  geom_point(aes(RX, RY, col = Community), plotdf, alpha = 0) +
  geom_path(aes(x, y), data = circ_dat1, linetype = "dashed")+
  facet_grid(rows = vars(Graph),
             cols = vars(Network.Size))+
  scale_color_manual(values=c("red", "blue"))+
  theme_bw()+
  labs(x = "X", y = "Y",
       col = "Community")

ggsave("residual_bounds.jpeg", 
       width = 7, height = 6, units = "in", 
       device = "jpeg", 
       path = "~/Documents/Work/github/BJSE/bias_simulation/figures/")

ggplot() +
  geom_point(aes(RX, RY, col = Community), plotdf,size = .01, alpha = .05) +
  geom_path(aes(x, y), data = circ_dat1, linetype = "dashed")+
  facet_grid(rows = vars(Graph),
             cols = vars(Network.Size))+
  scale_color_manual(values=c("red", "blue"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1)))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = "X", y = "Y",
       col = "Community")

ggsave("residuals.jpeg", 
       width = 6, height = 6, units = "in", 
       device = "jpeg", 
       path = "~/Documents/Work/github/BJSE/bias_simulation/figures/")


