#--------------------------------------
#
#       Second Moment Simulations
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
#    Set up Parameteric 
#     test statistic 
#-----------------------------

#set up SigmaD unknown
sigma2 <- function(x, y) as.numeric((crossprod(x,y) - crossprod(x,y)^2)) * tcrossprod(y, y)
Sigma_D_unknown <- function(Lhat,n,i){
  #set m 
  m <- nrow(Lhat)/n
  
  #set up scaling parameters
  S2D_hat_inv <- n * m * solve(crossprod(Lhat)) 
  
  #make omnibar
  omnibar <- 1/2 * (Lhat[1:n,] + Lhat[(n+1):(2*n),])
  
  #set up inner sum
  inner <- matrix(0, ncol = ncol(Lhat), nrow = ncol(Lhat))
  for(j in (1:n)[!(1:n %in% c(i))]){
    inner <- inner + sigma2(omnibar[i,], omnibar[j,])
  }
  
  #calculate Sigma
  Sigma <- 1/2 * S2D_hat_inv %*% (1/(n-1) * inner) %*% S2D_hat_inv
  
  #return result
  return(Sigma)
  
}

#gets single test statistic
test_stat <- function(X1, X2, Sigma, n){
  
  #create test statistic - numerically stable matrix inversion
  C <- chol(Sigma)
  test.stat <- n * crossprod(backsolve(C, X1 - X2, trans = TRUE))

  #return test statistics
  return(test.stat)
}

#get all test statistics  
get_test_stats <- function(Lhat, n, X1, X2){
  
  #get SigmaD for each vertex
  SigmaD_list <- lapply(1:n, function(x) Sigma_D_unknown(Lhat, n, x))
  
  #get test statistics
  test_stat_list <- as.numeric(lapply(1:n, function(x) test_stat(X1[x, ], X2[x,], SigmaD_list[[x]], n)))
  
  #return statisitcs
  return(test_stat_list)
}

#-----------------------------
#    Set up Semi Parameteric 
#     threshold
#-----------------------------

get_T_threshold <- function(X){
  
  T.stat <- numeric(500)
    
  for(i in 1:500){
    #sample A1, A2 
    A1 <- sampP(X)
    A2 <- sampP(X)
    
    #sample A1 and A2
    A1 <- sampP(P1)
    A2 <- sampP(P2)
    
    #make Omni and embedd for Lhat and Z
    Lhat <- ase(make_omni(list(A1, A2)), 2)
    
    #make X1 and X2
    X1 <- Lhat[1:(net_size[i]),]
    X2 <- Lhat[(net_size[i] + 1):(2 * net_size[i]),]
    
    #get T value
    T.stat[i] <- sum(apply(X1- X2, 1, norm2))
  }
  
  threshold <- unname(quantile(1:100, probs = .95))
  return(threshold)
}

#-----------------------------
#    Simulation Parameters
#-----------------------------

#set up base parameters 
net_size <- 200
t <- seq(0, 1, length.out = 11)[-11]
mc_runs <- 25 #number of iterations

#set up storage
here <- 1
df <- matrix(NA, nrow = mc_runs * net_size * length(t), ncol = 7)
colnames(df) <- c("network_size", "iter", "t", "test_stat", "omni_stat","reject", "vertex_name")
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
      
      #make X1 and X2
      X1 <- Lhat[1:(net_size[i]),]
      X2 <- Lhat[(net_size[i] + 1):(2 * net_size[i]),]
      
      #Get test statistics
      test_stats <- get_test_stats(Lhat, net_size[i], X1, X2)
      
      #get rejections yes/no
      test_rejs <- ifelse(test_stats > qchisq(.95, 2), 1, 0)
      
      #get omni test statistic
      omni_stats <- apply(X1 - X2, 1, norm2)
      
      #get yes/no based on T value
      thres <- get
        test_rejs <- ifelse(test_stats > qchisq(.95, 2), 1, 0)
      
      #store results
      df[(here:(here + net_size[i] - 1)), ] <- cbind(rep(net_size[i], net_size[i]),#network sizes
                                                 rep(k, net_size[i]), # ieration number
                                                 rep(t[j], net_size[i]), # t value
                                                 test_stats, # store test statistics
                                                 omni_stats, # store omni statistics
                                                 rejs, #store rejects yes/no
                                                 1:net_size[i] #store vertex label
                                                 )
      #update pointer 
      here <- here + net_size[i]
      
    }
  }
}

#-----------------------------
#
#    Visualize Results
#
#-----------------------------

#-----------------------------
#    Behavior of Test Stats
#-----------------------------
library(dplyr); library(reshape2)

#Vertex Level
plotdf <- as.data.frame(df) %>%
  group_by(t, vertex_name) %>% 
  summarize(t_stat= mean(test_stat),
            omni_stat= mean(omni_stat),
            reject = mean(reject)) %>%
  melt(id.vars = c(1, 2, 5))


ggplot(plotdf %>% filter(variable == 't_stat'), aes(t, value, 
                   col = variable, group = interaction(vertex_name, variable)))+
  geom_point(alpha = 1, size = .1)+
  geom_line(alpha = .05)+
  geom_hline(yintercept = qchisq(.95, 2), lty = 2)+
  theme_bw()+
  labs(x = "t", y = "Test Statistic", title = "Vertex Level: Test Statistic Average")

#Full Graph Level
plotdf <- as.data.frame(df) %>%
  group_by(t, iter) %>% 
  summarize(t_stat= mean(test_stat), 
            omni_stat= mean(omni_stat), 
            reject = mean(reject)) %>%
  melt(id.vars = c(1, 2, 5))


ggplot(plotdf %>% filter(variable == 't_stat'), aes(t, value, 
                   col = variable, group = interaction(iter, variable)))+
  geom_point(alpha = 1, size = .1)+
  geom_line(alpha = .05)+
  geom_hline(yintercept = qchisq(.95, 2), lty = 2)+
  theme_bw()+
  labs(x = "t", y = "Test Statistic", title = "Graph Level: Test Statistic Average")


#-----------------------------
#    Power curves
#-----------------------------
plotdf <- as.data.frame(df) %>%
  group_by(t, iter) %>% 
  summarize(power = mean(reject),
            std = sd(reject))

ggplot(plotdf, aes(t, power, group = iter, col = iter))+
  geom_point(alpha = 1, size = .1)+
  geom_line(alpha = .05)+
  theme_bw()+
  labs(x = "t", y = "Empirical Power", title = "Full Graph Hypothesis Testing Power Curve")




