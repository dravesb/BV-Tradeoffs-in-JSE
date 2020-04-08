#--------------------------------------
#
#       Power Analaysis Simulations
#           
#--------------------------------------

#load packages and source files
library(pacman,lib.loc = "/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
p_load(MCMCpack, ggplot2, dplyr, reshape2, MASS, foreach, doParallel)
setwd('/Users/benjamindraves/Documents/Work/github/BJSE/power_analysis/code/')
source('basic_functions.R')
source('difference_variance_functions.R')
source('model_setup_ER_to_SBM.R')
source('omni_semi_par_threshold.R')
source('test_statistic_functions.R')

#set up storage
here <- 1
df <- matrix(NA, nrow = length(net_size)*length(t)*mc_runs, ncol = 11)
colnames(df) <- c('iter_no','t', 'net_size',
                  'T', 'W', 'W_hat',
                  'T_Thres','W_Thres', 
                  'T_Rej',  'W_Rej', 'W_hat_Rej')

#set up new values for t
t <- rbind(seq(0, 0.5, length.out = 11), 
           seq(0, 0.3, length.out = 11), 
           seq(0, 0.2, length.out = 11))


#set pointer
counter <- 1

set.seed(1000)
#begin simulation
for(i in 1:length(net_size)){ # loop over network sizes
  for(j in 1:ncol(t)){ # loop over t values
    
    #set C list
    C_list <- list()
    C_list[[1]] <- diag(1, ncol = 2, nrow = 2)
    C_list[[2]] <- C(t[i, j])
    
    #get S list
    S_list <- get_Slist(C_list)
    
    #sample rows of X
    comm_ids <- rep(1:2, each = net_size[i]/2)
    X <- L[comm_ids, ]
    Xc <- X %*% sqrt(C(t[i,j]))
    
    #set up scaled latent positions
    X.scaled <- rbind(X, Xc)
    
    #fetch T_threshold
    T.thres <- get_T_threshold_parallel(X, X.scaled)
    
    #set up P matrices
    P1 <- tcrossprod(X)
    P2 <- tcrossprod(Xc)
    
    for(k in 1:mc_runs){
      
      #sample A1 and A2 
      A1 <- sampP(P1)
      A2 <- sampP(P2)
      
      #construct Omni + embed
      Atil <- make_omni(list(A1, A2))
      L_hat.here <- ase(Atil, 2)
      
      #align matrix
      L_hat <- procrustes(L_hat.here, X.scaled)$X.new
      
      #get statistics
      T.here <- get_T_stat(L_hat)
      W.here <- get_W_stat(L_hat, X, S_list, C_list)
      W_hat.here <- get_W_hat_stat(L_hat)
      
      #T and W threshold 
      #Xbar <- 0.5 *(L_hat[1:(net_size[i]),] + L_hat[(net_size[i]+1):(2*net_size[i]),])
      #T.thres <- get_T_threshold(Xbar, X.scaled)
      W.thres <- qchisq(1 - alpha, net_size[i]*ncol(X))
      
      #get reject/fail to reject
      T_rej <- ifelse(T.here > T.thres, 1, 0)
      W_rej <- ifelse(W.here > W.thres, 1, 0)
      W_hat_rej <- ifelse(W_hat.here > W.thres, 1, 0)
      
      
      #store 
      df[counter, ] <- c(k, #iteration
                         t[i,j], #t value
                         net_size[i], #net_size
                         T.here, #T_stat
                         W.here, #W_stat
                         W_hat.here, #W_hat_stat
                         T.thres, #T_thresold
                         W.thres, #W_thresold
                         T_rej, #T Reject 0/1
                         W_rej, #W Reject 0/1
                         W_hat_rej #W_hat Reject 0/1
      )
      #update counter and print update
      counter <- counter + 1
      if(counter %% 1000 == 0){
        print(round(counter/nrow(df), 3))
      }
    }
    }
}


#write out
write.csv(df, "~/Documents/Work/github/BJSE/power_analysis/data/plotting_df_variable_t.csv")

