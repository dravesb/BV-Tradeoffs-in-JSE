#--------------------------------------
#
#       One Dimensional MSE 
#           Simulations
#
#--------------------------------------

#load packages and source files
pacman::p_load(MCMCpack, ggplot2, dplyr, reshape2)
setwd('/Users/benjamindraves/Documents/Work/github/BJSE/one_dim_BV_tradeoff/code/')
source('basic_functions.R')
source('model_setup_ER.R')
source('theoretical_mse_functions.R')

#set up storage
here <- 1
df <- matrix(NA, nrow = 4*2*net_size*length(C)*mc_runs, ncol = 11)
colnames(df) <- c('iter_no','c','node_id', 'graph', 'method', 
                  'estimate', 'bias', 'mse',
                   'theo_bias','theo_var','theo_mse')


#set pointers
chunk_size <- net_size*2
start <- 1 
stop <- chunk_size

#set seed 
set.seed(1985)
for(j in 1:length(C)){
  #set C list
  C_list <- list()
  C_list[[1]] <- diag(1, ncol = 1, nrow = 1)
  C_list[[2]] <- diag(C[j], ncol = 1, nrow = 1)
  
  #get S list
  S_list <- get_Slist(C_list)
  
  #sample rows of X
  X <- matrix(L, nrow = net_size)
  Xc <- X %*% sqrt(C[j])
  
  #set up scaled latent positions
  X.scaled <- rbind(X, Xc)
  
  #set up P matrices
  P1 <- tcrossprod(X)
  P2 <- tcrossprod(Xc)
  
  for(i in 1:mc_runs){
    
    #sample A1 and A2 
    A1 <- sampP(P1)
    A2 <- sampP(P2)
    
    #----------------------
    #   ASE 
    #----------------------
    
    #embedd individually 
    X_hat.here <- ase(A1, 1)
    Y_hat.here <- ase(A2, 1)
    
    #align matrices
    X_hat <- procrustes(X_hat.here, X)$X.new
    Y_hat <- procrustes(Y_hat.here, Xc)$X.new
    
    #rowise bias
    X_bias <- X_hat - X
    Y_bias <- Y_hat - Xc
    
    #rowise mse
    X_mse <- apply(X_hat - X, 1, get_mse)
    Y_mse <- apply(Y_hat - Xc, 1, get_mse)
    
    #store 
    df[start:stop, ] <- cbind(rep(i, chunk_size), #iteration number
                              rep(C[j], chunk_size), #scaling
                              rep(1:net_size, 2), #node id
                              rep(c('Graph 1', 'Graph 2'), each = net_size), #graph id
                              rep('ASE', chunk_size), #Method
                              
                              c(X_hat, Y_hat), #estimate
                              c(X_bias, Y_bias), #observed bias
                              c(X_mse, Y_mse),  #observed MSE
                              
                              ase_bias(X, C_list), #theoretical bias
                              ase_var(X, C_list), #theoretical variance
                              ase_mse(net_size, X, C_list) # theoretical MSE
                              )
    
    
    #update pointers
    start <- stop + 1 
    stop <- start + chunk_size - 1
    
    #----------------------
    #   Abar 
    #----------------------
    
    #embedd A bar
    X_hat.here <- ase((A1 + A2)/2, 1)
    
    #align matrices
    X_hat <- procrustes(X_hat.here, X)$X.new
    
    #get rowwise bias 
    X_bias <- X_hat - X
    Y_bias <- X_hat - Xc
    
    #get rowwise mse
    X_mse <- apply(X_hat - X, 1, get_mse)
    Y_mse <- apply(X_hat - Xc, 1, get_mse)
    
    #store 
    df[start:stop, ] <- cbind(rep(i, chunk_size), #iteration number
                              rep(C[j], chunk_size), #scaling
                              rep(1:net_size, 2), #node id
                              rep(c('Graph 1', 'Graph 2'), each = net_size), #graph id
                              rep('Abar', chunk_size), #Method
                              
                              c(X_hat, X_hat), #estimate
                              c(X_bias, Y_bias), #observed bias
                              c(X_mse, Y_mse),  #observed MSE
                              
                              abar_bias(X, C_list), #theoretical bias
                              abar_var(X, C_list), #theoretical variance
                              abar_mse(net_size, X, C_list) # theoretical MSE
    )
    
    #update pointers
    start <- stop + 1 
    stop <- start + chunk_size - 1
    
    
    #----------------------
    #   Omni 
    #----------------------
    
    #construct Omni + embed
    Atil <- make_omni(list(A1, A2))
    L_hat.here <- ase(Atil, 1)
    
    #align matrices
    L_hat <- procrustes(L_hat.here, X.scaled)$X.new
    X_hat <- L_hat[1:(net_size),]
    Y_hat <- L_hat[(net_size+1):(2*net_size),]
    
    #get rowwise bias
    X_bias <- X_hat - X
    Y_bias <- Y_hat - Xc
    
    #getrowise mse
    X_mse <- apply(X_hat - X, 1, get_mse)
    Y_mse <- apply(Y_hat - Xc, 1, get_mse)
    
    #store 
    df[start:stop, ] <- cbind(rep(i, chunk_size), #iteration number
                              rep(C[j], chunk_size), #scaling
                              rep(1:net_size, 2), #node id
                              rep(c('Graph 1', 'Graph 2'), each = net_size), #graph id
                              rep('Omni', chunk_size), #Method
                              
                              c(X_hat, Y_hat), #estimate
                              c(X_bias, Y_bias), #observed bias
                              c(X_mse, Y_mse),  #observed MSE
                              
                              omni_bias(X, C_list, S_list), #theoretical bias
                              omni_var(X, C_list, S_list), #theoretical variance
                              omni_mse(net_size, X, C_list, S_list) # theoretical MSE
    )
    
    #update pointers
    start <- stop + 1 
    stop <- start + chunk_size - 1
    
    #----------------------
    #   Omnibar 
    #----------------------
    
    #get average
    Lmean <- .5*(X_hat + Y_hat)
    
    #get rowwise bias 
    X_bias <- Lmean - X
    Y_bias <- Lmean - Xc
    
    #get rowwise mse
    X_mse <- apply(Lmean - X, 1, get_mse)
    Y_mse <- apply(Lmean - Xc, 1, get_mse)
    
    #store 
    df[start:stop, ] <- cbind(rep(i, chunk_size), #iteration number
                              rep(C[j], 2*net_size), #scaling
                              rep(1:net_size, 2), #node id
                              rep(c('Graph 1', 'Graph 2'), each = net_size), #graph id
                              rep('Omnibar', chunk_size), #Method
                              
                              c(Lmean, Lmean), #estimate
                              c(X_bias, Y_bias), #observed bias
                              c(X_mse, Y_mse),  #observed MSE
                              
                              omnibar_bias(X, C_list, S_list), #theoretical bias
                              omnibar_var(X, C_list, S_list), #theoretical variance
                              omnibar_mse(net_size, X, C_list, S_list) # theoretical MSE
    )
    
    #update pointers
    start <- stop + 1 
    stop <- start + chunk_size - 1
    
    
    }      
    #print update
    print(j/length(C)) 
  
}


#------------------------------
#
#       Data Prep 
#
#-----------------------------

#restructure dataframe
plotdf <- data.frame(iter_no = as.factor(df[,1]),
                     c = as.numeric(df[,2]),
                     node_id = as.factor(df[,3]),
                     graph = as.factor(df[,4]),
                     method = as.factor(df[,5]),
                     estimate = as.numeric(df[,6]),
                     bias = as.numeric(df[,7]),
                     mse = as.numeric(df[,8]),
                     theo_bias = as.numeric(df[,9]),
                     theo_var = as.numeric(df[,10]),
                     theo_mse = as.numeric(df[,11])) %>% 
  group_by(node_id, graph, method, c) %>%
  summarize(average_mse = mean(mse),
            average_estimate = mean(estimate),
            average_bias = mean(bias),
            empirical_var = var(estimate),
            theo_bias = mean(theo_bias),
            theo_var = mean(theo_var)/net_size,
            theo_mse = mean(theo_mse))

#write out data
write.csv(plotdf, "~/Documents/Work/github/BJSE/one_dim_BV_tradeoff/data/plotting_df.csv")

