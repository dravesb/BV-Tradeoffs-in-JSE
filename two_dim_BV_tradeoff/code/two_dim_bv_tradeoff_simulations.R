#--------------------------------------
#
#       Multidimensional MSE 
#           Simulations
#
#--------------------------------------

#--------------------------------------
#
#       DON'T USE THIS FILE 
#       USE: two_dim_BV_tradeoff.R
#
#--------------------------------------



#set up storage
here <- 1
df <- matrix(NA, nrow = length(net_size)*length(t)*mc_runs, ncol = 19)
colnames(df) <- c("sim_number","net_size", "t",
                  "ase11", "ase21", "ase12", "ase22",
                  "abar11", "abar21", "abar12", "abar22",
                  "omni11", "omni21", "omni12", "omni22",
                  "omnibar11", "omnibar21", "omnibar12", "omnibar22")

#set seed 
set.seed(1985)

for(i in 1:length(net_size)){
  for(j in 1:length(t)){
    
    #sample rows of X
    samp <- rep(c(1, 2), each = net_size[i]/2)#sample(1:2, net_size[i], replace = TRUE) 
    X <- L[samp, ]
    Xc <- X %*% sqrt(C(t[j]))
    
    #set up scaled latent positions
    X.scaled <- rbind(X, Xc)
    
    #set up P matrices
    P1 <- tcrossprod(X)
    P2 <- tcrossprod(Xc)
    
    #iterate over mc_runs
    for(k in 1:mc_runs){
      
      #sample A1 and A2 
      A1 <- sampP(P1)
      A2 <- sampP(P2)
      
      #----------------------
      #   ASE MSE
      #----------------------
      
      #embedd individually 
      X_hat.here <- ase(A1, 2)
      Y_hat.here <- ase(A2, 2)
      
      #align matrices
      X_hat <- procrustes(X_hat.here, X)$X.new
      Y_hat <- procrustes(Y_hat.here, Xc)$X.new
        
      #getrowise mse
      X_mse <- apply(X_hat - X, 1, get_mse)
      Y_mse <- apply(Y_hat - Xc, 1, get_mse)
      
      #get average mse for each group and in graph 1
      ase_comm1_mse1 <- mean(X_mse[samp == 1])
      ase_comm2_mse1 <- mean(X_mse[samp == 2])
  
      #get average mse for each group and in graph 2
      ase_comm1_mse2 <- mean(Y_mse[samp == 1])
      ase_comm2_mse2 <- mean(Y_mse[samp == 2])
      
      #----------------------
      #   Abar MSE
      #----------------------
      
      #embedd A bar
      X_hat.here <- ase((A1 + A2)/2, 2)
      
      #align matrices
      X_hat <- procrustes(X_hat.here, X)$X.new
      
      #get rowwise mse
      X_mse <- apply(X_hat - X, 1, get_mse)
      Y_mse <- apply(X_hat - Xc, 1, get_mse)
      
      #get average mse for each group and in graph 1
      abar_comm1_mse1 <- mean(X_mse[samp == 1])
      abar_comm2_mse1 <- mean(X_mse[samp == 2])
      
      #get average mse for each group and in graph 2
      abar_comm1_mse2 <- Y_mse[1,]mean(Y_mse[samp == 1])
      abar_comm2_mse2 <- mean(Y_mse[samp == 2])
      
      #----------------------
      #   Omni MSE
      #----------------------
      
      #construct Omni + embed
      Atil <- make_omni(A1, A2)
      L_hat.here <- ase(Atil, 2)
      
      #align matrices
      L_hat <- procrustes(L_hat.here, X.scaled)$X.new
      X_hat <- L_hat[1:(net_size[i]),]
      Y_hat <- L_hat[(net_size[i]+1):(2*net_size[i]),]
      
      #getrowise mse
      X_mse <- apply(X_hat - X, 1, get_mse)
      Y_mse <- apply(Y_hat - Xc, 1, get_mse)
      
      #get average mse for each group and in graph 1
      omni_comm1_mse1 <- mean(X_mse[samp == 1])
      omni_comm2_mse1 <- mean(X_mse[samp == 2])
      
      #get average mse for each group and in graph 2
      omni_comm1_mse2 <- mean(Y_mse[samp == 1])
      omni_comm2_mse2 <- mean(Y_mse[samp == 2])
      
      #----------------------
      #   Omnibar MSE
      #----------------------
      
      #get average
      Lmean <- .5*(X_hat + Y_hat)
      
      #get rowwise mse
      X_mse <- apply(Lmean - X, 1, get_mse)
      Y_mse <- apply(Lmean - Xc, 1, get_mse)
      
      #get average mse for each group and in graph 1
      omnibar_comm1_mse1 <- mean(X_mse[samp == 1])
      omnibar_comm2_mse1 <- mean(X_mse[samp == 2])
      
      #get average mse for each group and in graph 2
      omnibar_comm1_mse2 <- mean(Y_mse[samp == 1])
      omnibar_comm2_mse2 <- mean(Y_mse[samp == 2])
      
      #----------------------
      #   Store Results
      #----------------------
      
      df[here, ] <- c(k, #simulation number
                      net_size[i], #network size
                      t[j], #store which C we're on
                      #     Mean Squared Errors
                      ase_comm1_mse1, ase_comm2_mse1, ase_comm1_mse2, ase_comm2_mse2, # ASE 
                      abar_comm1_mse1, abar_comm2_mse1, abar_comm1_mse2, abar_comm2_mse2, # Abar 
                      omni_comm1_mse1, omni_comm2_mse1, omni_comm1_mse2, omni_comm2_mse2, # Omni
                      omnibar_comm1_mse1, omnibar_comm2_mse1,
                      omnibar_comm1_mse2, omnibar_comm2_mse2 # Omnibar 
      )
      
      #update counter 
      here <- here + 1
      
      #print update
      print(k)
    }      
  
    
  }
}


#------------------------------
#
#       Data Prep 
#
#-----------------------------
library(dplyr); library(reshape2)

#melt data frame and add community, graph, and method labels
plotdf <- as.data.frame(df) %>% 
  melt(id.vars = c("sim_number", "net_size", "t")) %>%
  mutate(community = ifelse(variable %in% c("ase11", "ase12","abar11", "abar12","omni11", "omni12","omnibar11", "omnibar12"), "Community 1", "Community 2"),
         graph = ifelse(variable %in% c("ase11", "ase21","abar11", "abar21","omni11", "omni21","omnibar11", "omnibar21"), "Graph 1", "Graph 2"),
         Method = ifelse(variable %in% c("ase11", "ase21", "ase12", "ase22"), "ASE",
                         ifelse(variable %in% c("abar11", "abar21", "abar12", "abar22"), "Abar",
                                ifelse(variable %in% c("omni11", "omni21", "omni12", "omni22"), "Omni", "Omnibar")))) %>%
  group_by(net_size, t, community, graph, Method) %>%
  summarize(mse_se = sd(value),
    average_mse = mean(value)
    #,mse_se = sd(value)/sqrt(mc_runs + net_size/2)
            ) %>% 
  mutate(mse_se = mse_se/sqrt(mc_runs + net_size/2))


#------------------------------
#
#       Read in data 
#             + 
#       format data
#
#-----------------------------
setwd("~/Documents/Work/github/BJSE/two_dim_BV_tradeoff/data/")
#write.csv(plotdf, "plotdf.csv")
plotdf <- read.csv("plotdf.csv")[,-1]

#replace comm1 == 1, graph == 1, method ase with average across t
params1 <- plotdf %>% 
  filter(community == "Community 1", graph == "Graph 1", Method == "ASE") %>% 
  group_by(net_size) %>%
  summarize(m = mean(average_mse), s = mean(mse_se))
params2 <- plotdf %>% 
  filter(community == "Community 2", graph == "Graph 1", Method == "ASE") %>% 
  group_by(net_size) %>%
  summarize(m = mean(average_mse), s = mean(mse_se))

#update graph 1 comm1 parameters
plotdf$average_mse[which(plotdf$community == "Community 1" & plotdf$graph == "Graph 1" & plotdf$Method == "ASE")] <- params1$m
plotdf$mse_se[which(plotdf$community == "Community 1" & plotdf$graph == "Graph 1" &plotdf$Method == "ASE")] <- params1$s

#update graph 1 comm2 parameters
plotdf$average_mse[which(plotdf$community == "Community 2" & plotdf$graph == "Graph 1" & plotdf$Method == "ASE")] <- params2$m
plotdf$mse_se[which(plotdf$community == "Community 2" & plotdf$graph == "Graph 1" &plotdf$Method == "ASE")] <- params2$s


#------------------------------------------------
#
#         Make a priori lines data
#
#------------------------------------------------
#Variance matrices
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

#make apriori lines
apriori <- as.data.frame(matrix(NA, ncol = 7, nrow = 4 * length(t)))
here <-  1
colnames(apriori)[4:7] <- c("Abar", "ASE", "Omni", "Omnibar")


for(j in 1:length(t)){#drift parameter
  #get S_list/C_list
  S_list <- S(C(t[j]))
  C_list <- list(diag(2), C(t[j]))
  
  for(i in 1:2){#network loop
    for(k in 1:2){#community loop
      
      #ASE MSE 
      conj <- solve(.5 *crossprod(L) %*% C_list[[i]])
      ase_bias <- 0
      ase_var <- conj %*% Sigma_tilde(L[k,], i, C_list) %*% conj
      ase_mse_here <- sum(diag(ase_var)) / (net_size^2/2)
      
      #ABAR MSE
      Cbar <- Reduce("+", C_list)/length(C_list)
      conj <-  solve(.5 *crossprod(L) %*% Cbar)
      
      abar_bias <- (sqrt(Cbar) - sqrt(C_list[[i]])) %*% L[k,]
      abar_var <- conj %*% Sigma_tilde(L[k,], i, C_list) %*% conj
      
      abar_mse_here <- norm2(abar_bias) + (sum(diag(abar_var)) / net_size^2/2)
      
      #Omni MSE
      omni_bias <- (S_list[[i]] - sqrt(C_list[[i]])) %*% L[k,]
      omni_var <- Sigma(L[k,], i, C_list, S_list) 
      omni_mse_here <- norm2(omni_bias) + (sum(diag(omni_var)) / net_size^2/2)
      
      
      #Omnibar
      m <- length(S_list)
      Sbar <- Reduce("+", S_list)/m
      S2D.inv <- solve(.5 * crossprod(L) %*%  Reduce("+", lapply(S_list, function(x) x^2)))
      
      omnibar_bias <- (Sbar - sqrt(C_list[[i]])) %*% L[k,]
      omnibar_var <- .25 * S2D.inv %*% Reduce("+",  lapply(1:m, function(ind) (Sbar + m*S_list[[ind]]) %*% Sigma_tilde(L[k,], ind, C_list) %*% (Sbar + m*S_list[[ind]]))) %*%  S2D.inv
      omnibar_mse_here <- norm2(omnibar_bias)+ (sum(diag(omnibar_var)) / 1.5 * net_size)
      
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


apriori_mse <- apriori %>% melt(id.vars = 1:3)
colnames(apriori_mse) <- c("graph", "community", "t", "Method", "MSE")
apriori_mse$MSE <- as.numeric(apriori_mse$MSE)
apriori_mse$t<- as.numeric(apriori_mse$t)


apriori_mse <- read.csv("~/Documents/Work/github/BJSE/two_dim_BV_tradeoff/data/apriori_mse.csv")[,-1]

#------------------------------------------------
#
#         Make Figures
#
#------------------------------------------------

#plot figure
library(ggplot2)
ggplot() +
  geom_point(aes(t, average_mse, col = Method), plotdf, alpha = .5)+
  geom_line(aes(t, average_mse, col = Method), plotdf)+
  geom_line(aes(t, MSE, col = Method), apriori_mse, linetype = "dashed")+
  facet_grid(rows = vars(graph),
             cols = vars(community))+
  geom_ribbon(aes(t, ymin = average_mse - 1.96*mse_se,
                  ymax = average_mse + 1.96*mse_se),
              plotdf, alpha=0.1, linetype = 0)+
  theme_bw()+
  scale_y_log10()+
  labs(y = expression(paste('log'[10], "(MSE)")),
       x = "Deviation from SBM (x = 0) to ER (x = 1)") 


ggsave(filename = "2d_mse_mean.pdf", 
       width = 8, height =8, 
       path = "~/Documents/Work/github/BJSE/two_dim_BV_tradeoff/figures/", 
       units = "in")
  

ggplot(plotdf %>% filter(Method != "Omnibar"),
       aes(t, average_mse, col = Method)) +
  geom_point(alpha = .5, size = 2)+
  geom_line(size = 1)+
  facet_grid(rows = vars(graph),
             cols = vars(community))+
  geom_ribbon(aes(ymin= average_mse - 1.96*mse_se,
                  ymax= average_mse + 1.96*mse_se),
              alpha=0.1, linetype = 0)+
  theme_bw()+
  scale_y_log10()+
  labs(y = expression(paste('log'[10], "(MSE)")),
       x = "t")

ggsave(filename = "2d_mse_3_estimators_median.pdf", 
       width = 6, height = 6, 
       path = "~/Documents/Work/github/BJSE/two_dim_BV_tradeoff/figures/", 
       units = "in")





