#-----------------------------------------
#
#    Set up of Base Model
#
#-----------------------------------------

#source files
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/basic_functions.R")

#set up blocks
B <- matrix(c(.25, .05, .05, .25), byrow = T, nrow = 2)
b_ase <- ase(B, 2)

#get latent positions
x1 <- b_ase[1,]; x2 <- b_ase[2,]

#set prior probabilities
probs <- c(.5, .5)

#get rotation (eigenvectors of Delta)
Delta <- probs[1] * tcrossprod(x1) + probs[2]*tcrossprod(x2) 
R <- eigen(Delta)$vectors

#Apply rotation to x1 and x2
x1_til <- t(R) %*% x1 
x2_til <- t(R) %*% x2

#set base latent positions
L <- rbind(t(x1_til), t(x2_til))

#-----------------------------------------
#
#     Set up C/S values
#
#-----------------------------------------

#converges to ER with p = .3 (1,1) --> (2,0)
C <- function(t){
  diag(c(t + 1, 1-t))
}

S <- function(x){
  #bias matrices
  v1 <- c(1, C(x)[1,1])
  v2 <- c(1, C(x)[2,2])
  
  #get embeddings
  a1 <- ase(H1(v1), 1)[,1]
  a2 <- ase(H1(v2), 1)[,1]  
  
  #define scaling matrices
  S1 <- abs(diag(c(a1[1], a2[1])))
  S2 <- abs(diag(c(a1[2], a2[2])))
  
  #return results
  return(list(S1, S2))
}

#-----------------------------------------
#
#     Get Variance Values
#
#-----------------------------------------

#set up base parameters 
net_size <- 250
t <- seq(0, 1, length.out = 11)[-11]

#Variance matrices
sigma_2 <- function(y, C, k){
  p <- crossprod(y, C %*% L[k,])
  as.numeric(p - p^2)
}
Sigma_tilde <- function(y, g, C_list){
  
  #set up preliminaries
  K <- nrow(L)
  d <- ncol(L)
  
  #get sum
  tot <- matrix(0, nrow = d, ncol = d)
  for(k in 1:K){
    tot <- tot + (probs[k] * sigma_2(y, C_list[[g]], k) * tcrossprod(L[k,]))
    
  }
  return(tot)
  
}
Sigma <- function(y,g,C_list, S_list){
  
  #get S2D_inv
  S2D_inv <- solve(.5 * crossprod(L) %*%  Reduce("+", lapply(S_list, function(x) x^2)))
  
  #preliminary values
  d <- ncol(L)
  m <- length(C_list)
  
  #first summand
  tot1 <- (S_list[[g]] + Reduce("+",S_list)) %*% Sigma_tilde(y, g, C_list) %*% (S_list[[g]] + Reduce("+",S_list))
  
  #second summard
  tot2 <- matrix(0, nrow = d, ncol = d)
  for(k in (1:m)[-g]){
    tot2 <- tot2 + S_list[[k]] %*% Sigma_tilde(y, k, C_list) %*% S_list[[k]]
  }
  
  #return covariance matrix
  return(.25 * S2D_inv %*% (tot1 + tot2) %*% S2D_inv)
}

#-----------------------------------------
#
#     Make apriori line data
#
#-----------------------------------------


#make apriori lines
apriori <- as.data.frame(matrix(NA, ncol = 7, nrow = 4 * length(t)))
here <-  1
colnames(apriori)[4:7] <- c("Abar", "ASE", "Omni", "Omnibar")


for(j in 1:length(t)){#drift parameter
  #get S_list/C_list
  S_list <- S(t[j])
  C_list <- list(diag(2), C(t[j]))
  
  for(i in 1:2){#network loop
    for(k in 1:2){#community loop
      
      #ASE MSE 
      conj <- solve(.5 *crossprod(L) %*% C_list[[i]])
      ase_bias <- 0
      ase_var <- conj %*% Sigma_tilde(L[k,], i, C_list) %*% conj
      ase_mse_here <- ase_bias + sum(diag(ase_var)) / (net_size + net_size/2)
      
      #ABAR MSE
      Cbar <- Reduce("+", C_list)/length(C_list)
      conj <-  solve(.5 *crossprod(L) %*% Cbar)
      
      abar_bias <- (sqrt(Cbar) - sqrt(C_list[[i]])) %*% L[k,]
      abar_var <- conj %*% Sigma_tilde(L[k,], i, C_list) %*% conj
      
      abar_mse_here <- norm2(abar_bias) + (sum(diag(abar_var)) / (length(S_list)*(net_size + net_size/2)))
      
      #Omni MSE
      omni_bias <- (S_list[[i]] - sqrt(C_list[[i]])) %*% L[k,]
      omni_var <- Sigma(L[k,], i, C_list, S_list) 
      omni_mse_here <- norm2(omni_bias) + (sum(diag(omni_var)) / (net_size + net_size/2))
      
      #Omnibar
      m <- length(S_list)
      Sbar <- Reduce("+", S_list)/m
      S2D.inv <- solve(.5 * crossprod(L) %*%  Reduce("+", lapply(S_list, function(x) x^2)))
      
      omnibar_bias <- (Sbar - sqrt(C_list[[i]])) %*% L[k,]
      omnibar_var <- .25 * S2D.inv %*% Reduce("+",  lapply(1:m, function(ind) (Sbar + m*S_list[[ind]]) %*% Sigma_tilde(L[k,], ind, C_list) %*% (Sbar + m*S_list[[ind]]))) %*%  S2D.inv
      omnibar_mse_here <- norm2(omnibar_bias)+ (sum(diag(omnibar_var)) / (net_size + net_size/2))
        
      #store
      apriori[here,] <- c(paste("Graph", i), #graph
                          paste("Community", k), #community
                          t[j], #dirft parameter
                          #MSE's
                          abar_mse_here,ase_mse_here, omni_mse_here,omnibar_mse_here
                          #norm2(abar_bias), norm2(ase_bias), norm2(omni_bias), norm2(omnibar_bias)
                          )
      
      #update pointer 
      here <- here + 1
    }
    
  }
}

#-----------------------------------------
#
#	     Plot Results
#
#-----------------------------------------
library(ggplot2); library(dplyr); library(reshape2)


apriori_mse <- apriori %>% melt(id.vars = 1:3)
colnames(apriori_mse) <- c("graph", "community", "t", "Method", "MSE")
apriori_mse$MSE <- as.numeric(apriori_mse$MSE)
apriori_mse$t<- as.numeric(apriori_mse$t)

ggplot(apriori_mse, aes(t, MSE, col = Method)) +
  geom_point(alpha = .5)+
  geom_line()+
  facet_grid(rows = vars(graph),
             cols = vars(community))+
  theme_bw()+
  scale_y_log10()+
  labs(y = expression(paste('log'[10], "(MSE)")),
       x = "Deviation from SBM (t = 0) to ER (t = 1)") 

ggsave("analytic_mse.jpeg",
       width = 7, height = 5, units = "in", 
       device = "jpeg", 
       path = "~/Documents/Work/github/BJSE/two_dim_BV_tradeoff/figures/")


#-----------------------------------------
#
#	     Save Data
#
#-----------------------------------------

write.csv(apriori_mse, "~/Documents/Work/github/BJSE/two_dim_BV_tradeoff/data/apriori_mse.csv")
