#--------------------------------------
#
#       Full Graph Hypothesis 
#       Testing Simulation
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

#set base latent positions
L <- rbind(x1, x2)
#-----------------------------
#     Set up C values
#-----------------------------

#converges to ER with p = .5 (1,1) --> (1,0)
C <- function(t){
  diag(c(1, 1-t))
}

#-----------------------------
#    Set up Bias matrices 
#-----------------------------
get_S_list <- function(C_list){
  m <- length(C_list)
  Cm <- 1/m * Reduce('+', lapply(C_list, function(C) C^2))
  Cm_forth <- diag(diag(Cm)^(1/4))
  Cm_forth_inv <- diag(diag(Cm)^(-1/4))
  S_plus_list <- lapply(C_list, function(C) 0.5 * (C %*% Cm_forth_inv + Cm_forth))
  return(S_plus_list)              
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
#    Set up parametric test
#----------------------------
Sigma_D_known <- function(vertex, samp){
  # get m 
  m <- nrow(L)
  
  # Delta inverse
  Delta_inv <- solve(0.5*(tcrossprod(L[1,]) + tcrossprod(L[2,])))
 
  # Sigma tilde
  f <- function(z){
    x1 <- matrix(L[1,], ncol = 1)
    x2 <- matrix(L[2,], ncol = 1)
    
    return( 0.5 * ((as.numeric(crossprod(z, x1)-(crossprod(z, x1))^2)) * tcrossprod(x1) +
                   (as.numeric(crossprod(z, x2)-(crossprod(z, x2))^2)) * tcrossprod(x2)))
  }
  
  Sigma_tilde_1 <- f(matrix(L[1,], ncol = 1))
  Sigma_tilde_2 <- f(matrix(L[2,], ncol = 1))
  
  # return based on group assignment
  if(samp[vertex] == 1){
    return(0.5 * Delta_inv %*% Sigma_tilde_1 %*% Delta_inv)
  }else{
    return(0.5 * Delta_inv %*% Sigma_tilde_2 %*% Delta_inv)
  }
  
}
Sigma_D_unknown <- function(Lhat, n, m, vertex, method = 1){
  # Delta hat inverse 
  Del_hat_inv <- solve(crossprod(Lhat) / (n*m)) 
  
  # Sigma tilde hat
  omnibar <- 0.5 * (Lhat[1:n, ] + Lhat[(n+1):(2*n),])
  
  # Calculate Sigma Tilde Hat
  if(method == 1){
    f1 <- function(i, j){
      x <- as.numeric(crossprod(omnibar[i,], omnibar[j,]))
      return((x - x^2) * tcrossprod(omnibar[j,]))
    }  
    Sigma_tilde_hat <- (1/(n-1)) * Reduce('+', lapply((1:n)[-vertex], function(j) f1(vertex, j)))  
  }else if(method == 2){
    f2 <- function(i, j){
      x <- as.numeric(crossprod(Lhat[i,], Lhat[j,]))
      return((x - x^2) * tcrossprod(Lhat[j,]))
    }  
    Sigma_tilde_hat <- (1/(n*m)) * Reduce('+', lapply(1:(n*m), function(j) f2(vertex, j)))
  }else{
    f3 <- function(i, j){
      x <- as.numeric(crossprod(omnibar[i,], Lhat[j,]))
      return((x - x^2) * tcrossprod(Lhat[j,]))
    }
    Sigma_tilde_hat <- (1/(n*m)) * Reduce('+', lapply(1:(n*m), function(j) f3(vertex, j)))  
  }
  
  # return variance
  Sigma_D_hat <- 0.5 * Del_hat_inv %*% Sigma_tilde_hat %*% Del_hat_inv
  return(Sigma_D_hat)
  
}
get_W_stat <- function(Lhat, n, m, known_cov = FALSE, method = 1){
  
  #fetch variances 
  if(known_cov == TRUE){
    Sigma_D_list <- lapply(1:n, function(i) Sigma_D_known(i, samp))
  }else{
    Sigma_D_list <- lapply(1:n, function(i) Sigma_D_unknown(Lhat, n, m, i, method))
  }
  
  # Create Wald Statistics
  node_differences <- lapply(1:n, function(i) matrix(Lhat[i, ] - Lhat[n + i, ], nrow = 1))
  W_i <- lapply(1:n, function(i) mahalanobis(node_differences[[i]], center = FALSE, cov = Sigma_D_list[[i]]))
  
  #Sigma_D_inv_list <- lapply(Sigma_D_list, function(M) ginv(psd_proj(M)))
  #W_i <- lapply(1:n, function(i) crossprod(node_differences[[i]], Sigma_D_inv_list[[i]]) %*% node_differences[[i]])
  
  # Return aggregate statistic
  return(n*Reduce('+', W_i))
}

#-----------------------------
#    Set up Semi Parameteric 
#     Testing
#-----------------------------
get_T_threshold <- function(X_known, no.iters = 500){
  
  T.stat <- numeric(no.iters)
  
  for(i in 1:no.iters){
    # get P1, P2 
    P <- tcrossprod(X_known); P[which(P < 0)] <- 0; P[which(P > 1)] <- 1
    
    #sample A1, A2 
    A_list <- lapply(1:2, function(x) sampP(P))
    
    #make Omni and embedd for Lhat and Z
    Lhat <- ase(make_omni(A_list), 2)
    
    #make X1 and X2
    X1 <- Lhat[1:(nrow(X_known)),]
    X2 <- Lhat[(nrow(X_known) + 1):(2 * nrow(X_known)),]
    
    #get T value
    T.stat[i] <- sum(apply(X1-X2, 1, function(x) sum(x^2)))
  }
  
  threshold <- unname(quantile(T.stat, probs = .95))
  return(threshold)
}
get_T_stat <- function(Lhat, n){
  Reduce('+', lapply(1:n, function(i) sum((Lhat[i,] - Lhat[n+i,])^2)))
}


#-----------------------------
#    Simulation Parameters
#-----------------------------

#simulation parameters
net_size <- c(50, 100, 200)
t_max <- c(1, 0.75, 0.5)
mc_runs <- 1000

#-----------------------------
#    Parametric simulation
#-----------------------------

#storage 
df_list <- list()
count <- 1
set.seed(1985)
# begin simulation
for(i in 1:length(net_size)){
  # set up t_seq
  n <- net_size[i]
  m <- 2; d <- 2
  t_seq <- seq(0, t_max[i], length.out = 5)
  
  
  for(j in 1:length(t_seq)){
    # set up model
    C_list <- list(diag(2), C(t_seq[j]))
    
    for(k in 1:mc_runs){
      
      # sample and embedd 
      samp <- sample(1:2, size = n, replace = TRUE) 
      X <- L[samp,] 
      A_list <- lapply(1:m, function(g) sampP(X %*% tcrossprod(C_list[[g]], X)))
      Lhat <- ase(make_omni(A_list), d)
      
      # procrustes for known covariance - not needed for estimated
      Lhat_rot <- procrustes(Lhat, rbind(X, X %*% sqrt(C_list[[2]])))$X.new
      
      # fetch statistics
      W_stat <- get_W_stat(Lhat_rot, n, m, known_cov = TRUE) 
      W_hat_stat <- get_W_stat(Lhat, n, m, known_cov = FALSE, method = 3)
      
      # get accept / reject
      W_stat_decision <- ifelse(W_stat > qchisq(0.95, n*d), 1, 0)
      W_hat_stat_decision <- ifelse(W_hat_stat > qchisq(0.95, n*d), 1, 0)
      
      #store results 
      df_list[[count]] <- data.frame(net_size = n, 
                                     m = m, 
                                     iter = k, 
                                     t = t_seq[j],
                                     Test_Statistic = c('W', 'W_hat'),
                                     Data_Dependent = c(FALSE, TRUE),
                                     Statistic = c(W_stat, W_hat_stat), 
                                     Threshold = c(qchisq(0.95, n*d), qchisq(0.95, n*d)), 
                                     Decision = c(W_stat_decision, W_hat_stat_decision))
      count <- count + 1
      
    }
    print(paste('n =', n,'; j = ', j, 'is finished.'))
  }
  print(paste('n =', n, 'is finished.'))
}
df_parm <- Reduce('rbind', df_list) 
write.csv(df_parm, file = '~/Documents/Work/github/BJSE/power_analysis/data/parameteric_test_results.csv')


#-----------------------------
#    Semi-parameteric Simulation
#-----------------------------

#storage 
df_list <- list()
count <- 1
set.seed(1985)
# begin simulation
for(i in 1:length(net_size)){
  # set up t_seq
  n <- net_size[i]
  m <- 2; d <- 2
  t_seq <- seq(0, t_max[i], length.out = 5)
  
  # set up semi-parametric threshold based on even split
  #thres_known <- thres_unknown <-  get_T_threshold(X_known = L[rep(1:2, each = n/2),])
  
  for(j in 1:length(t_seq)){
    # set up model
    C_list <- list(diag(2), C(t_seq[j]))
    
    for(k in 1:mc_runs){
      
      # sample and embedd 
      samp <- sample(1:2, size = n, replace = TRUE) 
      X <- L[samp,] 
      A_list <- lapply(1:m, function(g) sampP(X %*% tcrossprod(C_list[[g]], X)))
      Lhat <- ase(make_omni(A_list), d)
    
      # get unknown threshold based on A1
      thres_unknown <- get_T_threshold(X_known = ase(A_list[[1]], 2), no.iters = 100)
      
      # fetch statistics
      T_stat <- get_T_stat(Lhat, n)
      
      # get accept / reject
      #T_stat_decision <- ifelse(T_stat > thres_known, 1, 0)
      T_hat_stat_decision <- ifelse(T_stat > thres_unknown, 1, 0)
      
      #store results 
      df_list[[count]] <- data.frame(net_size = n, 
                                     m = m, 
                                     iter = k, 
                                     t = t_seq[j],
                                     Test_Statistic = c('T'),
                                     Data_Dependent = c(TRUE),
                                     Statistic = c(T_stat), 
                                     Threshold = c(thres_unknown), 
                                     Decision = c(T_stat_decision))
      count <- count + 1
      
    }
    print(paste('n =', n,'; j = ', j, 'is finished.'))
  }
  print(paste('n =', n, 'is finished.'))
}

df_semi_parm <- Reduce('rbind', df_list) 
write.csv(df_semi_parm, file = '~/Documents/Work/github/BJSE/power_analysis/data/semi_parameteric_unkown_lp_test_results.csv')


#-----------------------------
#    Visualize Results
#-----------------------------


df_semi_parm %>% 
  group_by(net_size, t, Test_Statistic) %>% 
  summarize(Empirical_Power = mean(Decision),
            #lower = mean(Decision) - 2 * sd(Decision)/sqrt(n()), 
            #upper = mean(Decision) + 2 * sd(Decision)/sqrt(n())
            lower = Empirical_Power - 2 *sqrt(Empirical_Power*(1 - Empirical_Power)/sqrt(n())), 
            upper = mean(Decision) + 2 *sqrt(Empirical_Power*(1 - Empirical_Power)/sqrt(n()))
            
            
            ) %>%
  ggplot(aes(x = t, y = Empirical_Power, col = Test_Statistic, group = Test_Statistic)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              fill = 'grey70', alpha = .5, colour = NA) +
  geom_point(alpha = 0.5) + geom_line() + 
  geom_hline(yintercept = 0.05, linetype = 'dashed', col = 'grey') +
  scale_color_manual(labels = c('T',expression(hat('T')), 'W', expression(hat('W'))), 
                     values = c('red','blue','purple', 'darkgreen'))+
  facet_grid(~net_size, scales = 'free', labeller = label_bquote(rows = 'n ='~.(net_size))) +
  labs(x = 't', y = 'Empirical Power', col = 'Statistic') + 
  theme_bw()

df %>% 
  dplyr::select(-Decision) %>% 
  filter(Test_Statistic == 'W_hat', net_size == 50) %>% 
  melt(id.vars = c(1:5)) %>% 
  ggplot(aes(x = value, fill = variable)) + 
    geom_histogram(alpha = 0.8, position = 'dodge') + 
    facet_grid(net_size~t) + 
    theme_bw()



