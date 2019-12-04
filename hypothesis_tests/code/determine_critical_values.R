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
net_size <- c(100, 250, 500)
t <- 0 #c(0, .25, .5, .75, 1)
mc_runs <- 100 #number of iterations

#set up storage
df <- matrix(NA, nrow = length(net_size) * mc_runs * length(t) , ncol = 9)
colnames(df) <- c("network_size", "iter", "t",
                  "T_known", "T_unknown", "T_Omni",
                  "t_known_rej", "test_rej", "omni_rej")
here <- 1 

#set seed 
set.seed(1985)

for(i in 1:length(net_size)){#iterate over network size
  
  #set up threshold on equal group sized X
  X_known <- L[c(rep(1, net_size[i]/2), rep(2, net_size[i]/2)),]
  thres <- get_T_threshold(X_known)
  
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
      
      #known variance
      T_star0 <- sum(get_test_stats_known_variance(samp, L, X1, X2, S_list, C_list)) #statistic
      test_rej_known <- ifelse(T_star0 > qchisq(.95, 2 * net_size[i]), 1, 0) #get rejections yes/no
      
      #unknown variance
      T_star1 <- sum(get_test_stats(Lhat, net_size[i], X1, X2)) #statistic
      test_rej <- ifelse(T_star1 > qchisq(.95, 2 * net_size[i]), 1, 0) #get rejections yes/no
      
      #get omni test statistic & threshold
      T_star2 <- sum(apply(X1 - X2, 1, norm2)) #get statistic
      omni_rej <- ifelse(T_star2 > thres, 1, 0) #get yes/no based on threshold value
      
      #store results
      df[here, ] <- c(net_size[i], k, t[j], 
                      T_star0,T_star1,T_star2, 
                      test_rej_known, test_rej, omni_rej)
      
      #update pointer 
      here <- here + 1
    }
    print(j)
  }
}

#-----------------------------
#
#    Visualize Results
#
#-----------------------------

#-----------------------------
#    Power Curves
#-----------------------------
library(dplyr); library(reshape2)

#calculate standard errors
se <- as.data.frame(reject_df) %>%
  group_by(t, network_size) %>% 
  summarize(Par_Known_se = 2 * sd(t_known_rej) / sqrt(n()),
            Par_Unknown_se = 2 * sd(test_rej) / sqrt(n()),
            SemiPar_se = 2 * sd(omni_rej) / sqrt(n())
  ) %>% 
  melt(id.vars = c(1,2))

#calculate means
plotdf <- as.data.frame(reject_df) %>%
  group_by(t, network_size) %>% 
  summarize(Par_Known = mean(t_known_rej),
            Par_Unknown = mean(test_rej),
            SemiPar = mean(omni_rej)
  ) %>% 
  melt(id.vars = c(1,2)) %>%
  mutate(se = as.numeric(se[,4]))

#Change variable name

#visualize
ggplot(plotdf, aes(t, value, col = variable, group = variable))+
  geom_point(alpha = 1, size = .1)+
  geom_line(alpha = 1)+
  geom_ribbon(aes(ymin = value - se, ymax = value + se), 
              fill = "grey70", alpha = .2, colour = NA)+ 
  geom_hline(yintercept = .05, col = "grey", linetype = "dashed")+
  facet_grid(rows = vars(network_size))+
  theme_bw()+
  labs(x = "t", y = "Empirical Power", 
       title = "Full Graph Hypothesis Testing Power Curve",
       color = "Test\n Statistic"
  )
#-----------------------------
#    Save Dataframe
#-----------------------------
write.csv(df, "~/Documents/Work/github/BJSE/hypothesis_tests/data/check_critical_value.csv")

df <- read.csv("~/Documents/Work/github/BJSE/hypothesis_tests/data/check_critical_value.csv")[,-1]
#-----------------------------
#    Visualize Test Stat
#-----------------------------

plotdf <- as.data.frame(df[,1:6]) %>%
  melt(id.vars = 1:3)

ggplot(as.data.frame(plotdf) %>% filter(variable != "T_Omni"), 
       aes(value, col = variable))+
  geom_histogram(fill = "white")+
  facet_grid(rows = vars(variable),
             cols = vars(network_size),
             scales="free")+
  theme_bw()

ggsave("test_stat_under_null.jpeg",
       height = 6, width = 6,
       units = "in",
       path = "~/Documents/Work/github/BJSE/hypothesis_tests/figures/")


#-----------------------------
#    Determine Critical Value
#-----------------------------
size <- net_size[3]

t_known <- quantile((plotdf %>% filter(variable == "T_known", network_size == size))[,"value"], probs = .95)
t_unknown <- quantile((plotdf %>% filter(variable == "T_unknown", network_size == size))[,"value"], probs = .95)
CV <- qchisq(.95, 2*size)

c(t_known, t_unknown, CV)

sum((plotdf %>% filter(variable == "T_known", network_size == size))[,"value"] < qchisq(.95, 2*size))/size
sum((plotdf %>% filter(variable == "T_unknown", network_size == size))[,"value"] < qchisq(.95, 2*size))/size


#Take away: the unknown variance estimator is overinflated so we're testing based on 
#           an improper reference distribution
#Next up: Do a simulation that relates the reference distribution to the observed 
#         95% intervals and try to tease out the correction rate

#-----------------------------
#    Visualize Power
#-----------------------------

powerdf <- df %>% 
  group_by(network_size, t) %>% 
  summarize(known_power = mean(t_known_rej), 
            unknown_power = mean(test_rej),
            omni_power = mean(omni_rej)) %>% 
  melt(id.vars = 1:2)

ggplot(powerdf, aes(network_size, value, col = variable))+
  geom_point(size = 2)+
  geom_line()+
  geom_hline(yintercept = .05, col = "grey", linetype = "dashed")+
  labs(x = "Network Size", 
       y = "Power", 
       color = "Method")+
  theme_bw()

ggsave("null_power_few.jpeg",
       height = 6, width = 6,
       units = "in",
       path = "~/Documents/Work/github/BJSE/hypothesis_tests/figures/")








