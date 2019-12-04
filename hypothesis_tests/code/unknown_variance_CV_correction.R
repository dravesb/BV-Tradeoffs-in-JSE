#--------------------------------------
#
#      Determine Critical Values
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
probs <- c(.5, .5)

#get rotation (eigenvectors of Delta)
Delta <- probs[1] * tcrossprod(x1) + probs[2]*tcrossprod(x2) 
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
  S1 <- abs(diag(c(a1[1], a2[1])))
  S2 <- abs(diag(c(a1[2], a2[2])))
  
  #return results
  return(list(S1, S2))
}

#--------------------------------------
#    Set up helpful matrix functions
#--------------------------------------
library(MASS)
psd_proj <- function(M){
  #decomp 
  eigen_decomp <- eigen(M)
  U <- eigen_decomp$vectors
  S <- eigen_decomp$values
  
  #threshold eigenvalues
  S.thres <- diag(sapply(S, function(x) max(c(0,x))))
  
  #return
  return(tcrossprod(U %*% S.thres, U))
  
}

#-----------------------------
#    Set up variance 
#     matrices
#-----------------------------
sigma2 <- function(x, y) as.numeric((crossprod(x,y) - crossprod(x,y)^2)) * tcrossprod(y, y)
sigma_2 <- function(y, C, k){
  as.numeric(crossprod(y, C %*% L[k,]) - crossprod(y, C %*% L[k,])^2)
}
Sigma_tilde <- function(y, g, C_list){
  
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


#set up SigmaD unknown
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
  
  #project onto psd cone
  Sigma_proj <- psd_proj(Sigma)
  
  #return result
  return(Sigma_proj)
}
Sigma_D_known <- function(n, samp, L, S_list, C_list){
  
  #set S_sum
  S_sum <- Reduce("+", S_list)
  #set Delta
  Del.inv <- solve(.5 * crossprod(L) %*%  Reduce("+", lapply(S_list, function(x) x^2)))
  
  #get SigmaD for each vertex
  SigmaD_list <- list()
  for(i in 1:n){
    if(samp[i] == 1){
      SigmaD_list[[i]] <- Del.inv %*% S_sum %*% (Sigma_tilde(L[1,], 1, C_list) + Sigma_tilde(L[1,], 2, C_list) ) %*% S_sum %*% Del.inv / 4 
    }else{
      SigmaD_list[[i]] <- Del.inv %*% S_sum %*% (Sigma_tilde(L[2,], 1, C_list) +  Sigma_tilde(L[2,], 2, C_list) ) %*% S_sum %*% Del.inv / 4 
    }
  }
  
  return(SigmaD_list)
  
}
Sigma_D_known2 <- function(n, samp, L, S_list, C_list){
  
  #set Delta
  Del.inv <- solve(.5 * crossprod(L))
  
  #get SigmaD for each vertex
  SigmaD_list <- lapply(1:n, function(i) Del.inv %*% Sigma_tilde(L[samp[i],], g = 1, list(diag(2), diag(2))) %*% Del.inv / 2)
  
  return(SigmaD_list)
  
}

#-----------------------------
#    Set up Parameteric 
#     test statistic 
#-----------------------------

#gets single test statistic
test_stat <- function(X1, X2, Sigma, n){
  
  #create test statistic - numerically stable matrix inversion
  C <- chol(Sigma)
  test.stat <- n * crossprod(backsolve(C, X1 - X2, trans = TRUE))
  
  #return test statistics
  return(test.stat)
}
test_stat_unstable <- function(x1, x2, Sigma, n){
  
  #create test statistic - numerically unstable matrix inversion - generalized inverse
  test.stat <- n * crossprod(x1 - x2, ginv(Sigma) %*% (x1 - x2))
  #return test statistics
  return(test.stat)
}

#get all test statistics  
get_test_stats <- function(Lhat, n, X1, X2){
  
  #get SigmaD for each vertex
  SigmaD_list <- lapply(1:n, function(x) Sigma_D_unknown(Lhat, n, x))
  
  #get test statistics
  test_stat_list <- as.numeric(lapply(1:n, function(x) test_stat_unstable(X1[x, ], X2[x,], SigmaD_list[[x]], n)))
  
  #return statisitcs
  return(test_stat_list)
}
get_test_stats_known_variance <- function(samp, L, X1, X2, S_list, C_list){
  
  #set n 
  n <- length(samp)
  
  #get variances
  Sigma_D_list <- Sigma_D_known(n, samp, L, S_list, C_list)
  
  #get test statistics
  test_stat_list <- as.numeric(lapply(1:n, function(x) test_stat_unstable(X1[x, ], X2[x,], Sigma_D_list[[x]], n)))
  
  #return statisitcs
  return(test_stat_list)
}

#-----------------------------
#    Set up Semi Parameteric 
#     threshold
#-----------------------------
get_T_threshold <- function(X_known){
  
  no.iters <- 50
  T.stat <- numeric(no.iters)
  
  for(i in 1:no.iters){
    #sample A1, A2 
    A1 <- sampP(tcrossprod(X_known))
    A2 <- sampP(tcrossprod(X_known))
    
    #make Omni and embedd for Lhat and Z
    Lhat <- ase(make_omni(list(A1, A2)), 2)
    
    #make X1 and X2
    X1 <- Lhat[1:(nrow(X_known)),]
    X2 <- Lhat[(nrow(X_known) + 1):(2 * nrow(X_known)),]
    
    #get T value
    T.stat[i] <- sum(apply(X1- X2, 1, norm2))
  }
  
  threshold <- unname(quantile(T.stat, probs = .95))
  return(threshold)
}

#-----------------------------
#    Simulation Parameters
#-----------------------------
library(MCMCpack)

#set up base parameters 
net_size <- round(exp(seq(log(100), log(1000), length.out = 20)))
t <- 0 #c(0, .25, .5, .75, 1)
mc_runs <- 100 #number of iterations

#set up storage
df <- matrix(NA, nrow = length(net_size) * mc_runs * length(t) , ncol = 4)
colnames(df) <- c("network_size", "iter", "t", "T_unknown")
here <- 1 

#set seed 
set.seed(1985)

for(i in 1:length(net_size)){#iterate over network size
  
  for(j in 1:length(t)){#iterate over t values
    
    #set S_list and C_list
    S_list <- S(C(t[j]))
    C_list <- list(diag(2), C(t[j]))
    
    for(k in 1:mc_runs){#number of MC simulations
      
      #sample group assignments
      samp <- sample(1:2, size = net_size[i], replace = TRUE)
      
      #set up X and P1 and P2
      X <- L[samp,]
      P1 <- tcrossprod(X)
      P2 <- tcrossprod(X, X %*% C(t[j]))
      
      #define Xtil 
      Xtil <- rbind(X, X %*% sqrt(C(t[j])))
      
      #sample A1 and A2
      A1 <- sampP(P1)
      A2 <- sampP(P2)
      
      #make Omni and embedd for Lhat and Z
      Lhat <- ase(make_omni(list(A1, A2)), 2)
      
      #procrustes!
      Lhat <- procrustes(Lhat, Xtil)$X.new
      
      #make X1 and X2
      X1 <- Lhat[1:(net_size[i]),]
      X2 <- Lhat[(net_size[i] + 1):(2 * net_size[i]),]
      
      #---------------------
      #Get test statistics
      #---------------------
      #unknown variance
      T_star1 <- sum(get_test_stats(Lhat, net_size[i], X1, X2)) #statistic
  
      #store results
      df[here, ] <- c(net_size[i], k, t[j], T_star1)
      
      #update pointer 
      here <- here + 1
    }
    print(j)
  }
}

#-----------------------------
#    Save Data
#-----------------------------

#write.csv(df, "~/Documents/Work/github/BJSE/hypothesis_tests/data/determine_critical_value.csv")
dat <- read.csv("~/Documents/Work/github/BJSE/hypothesis_tests/data/determine_critical_value.csv")[,-1]

#-----------------------------
#
#    Visualize Results
#
#-----------------------------

#-----------------------------
#    Comparing Cutoffs
#-----------------------------
library(dplyr); library(reshape2)

plotdf <- as.data.frame(dat) %>%
  group_by(t, network_size) %>% 
  summarize(cutoff = quantile(T_unknown, probs = c(.95))
  ) %>% 
  mutate(chisq_value = sapply(net_size, function(x) qchisq(.95, 2*x)))

ggplot(plotdf, aes(network_size, cutoff - chisq_value))+
  geom_point(size = 2)+
  geom_line(alpha = .2)+
  theme_bw()+
  labs(y = "Observed minus Theoretical Quantile", 
       x = "Network Size")

ggsave("obs_min_theo.jpeg",
       height = 6, width = 6,
       units = "in",
       path = "~/Documents/Work/github/BJSE/hypothesis_tests/figures/")

#bootstrap quantiles
boot <- function(x, B){replicate(B, {
    samp <- sample(1:length(x), size= length(x), replace = TRUE)
    quant <- unname(quantile(x[samp], probs = .95))
    return(quant)
  })}

B <- 1000 
quants <- numeric()
for(i in 1:length(net_size)){
  quants[((i-1)*B +1):(i*B)] <- boot(dat[dat$network_size == net_size[i], 4], B) - qchisq(.95, 2*net_size[i])
}
plotdf <- data.frame(quants, network_size = rep(net_size, each = B))

ggplot(plotdf, aes(network_size, quants))+
  geom_jitter(size = .2, alpha = .1, height = 5, width = 5)+
  theme_bw()+
  labs(y = "Observed minus Theoretical Quantile", 
       x = "Network Size")

ggsave("obs_min_theo_boot.jpeg",
       height = 6, width = 6,
       units = "in",
       path = "~/Documents/Work/github/BJSE/hypothesis_tests/figures/")

#-----------------------------
#     Visualize Power
#-----------------------------

powerdf <- as.data.frame(dat) %>%
  mutate(rej = ifelse(dat$T_unknown > qchisq(.95, 2 * dat$network_size) , 1, 0)) %>% 
  group_by(network_size) %>%
  summarize(power = mean(rej))
  
ggplot(powerdf, aes(network_size, power)) +
  geom_point(size = 2)+
  geom_line(col = "grey")+
  theme_bw()

ggsave("power_unknown_null.jpeg",
       height = 6, width = 6,
       units = "in",
       path = "~/Documents/Work/github/BJSE/hypothesis_tests/figures/")


#bootstrap quantiles
B <- 1000 
quants <- numeric()
for(i in 1:length(net_size)){
  quants[((i-1)*B +1):(i*B)] <- boot(dat[dat$network_size == net_size[i], 4], B)
}

quantdf <- data.frame(quants, network_size = rep(net_size, each = B)) %>%
  mutate(rej = ifelse(quants > qchisq(.95, 2 * network_size) , 1, 0)) %>% 
  group_by(network_size) %>%
  summarize(quant = mean(rej))


ggplot(quantdf, aes(network_size, quant))+
  geom_point(size = 2)+
  geom_line(col = "grey")+
  geom_hline(yintercept = .05, col = "grey", linetype = "dashed")+
  labs(x = "Network Size", 
       y = "Proportion of Quantiles\n Above Cuttoff", 
       title = "Bootstrapped Cuttoff Values")+
  theme_bw()

ggsave("cutoff_unknown_null_boot.jpeg",
       height = 6, width = 6,
       units = "in",
       path = "~/Documents/Work/github/BJSE/hypothesis_tests/figures/")

#-----------------------------
#     Visualize Test Stat
#-----------------------------

ggplot(dat, aes((T_unknown - 2*network_size)/sqrt(4*network_size))) + 
  geom_histogram(fill = "red", col = "black", alpha = .5, bins = 50)+
  facet_grid(rows = vars(network_size))+
  geom_vline(xintercept = qnorm(.95))+
  theme_bw()+
  labs(x = "Wald Statistic")
  
test <- (dat$T_unknown - 2*dat$network_size)/sqrt(4*dat$network_size)

level <- numeric(length(net_size))
for(i in 1:length(net_size)){
  level[i] <- length(which(test[dat$network_size == net_size[i]] > qnorm(.95)))/net_size[i]
}
level
plot(net_size, level, type = "l")









