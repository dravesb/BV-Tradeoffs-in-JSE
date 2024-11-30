#------------------------------
#
#   Clustering Example 
#    Omni v Abar v ASE
#
#------------------------------

#source libraries + functions
source('/Users/benjamindraves/Desktop/OmniAnalysis/code_base/loadAll.R')
fig_path <- '/Users/benjamindraves/Documents/Work/github/BJSE/Motivating_Ex_Comm_Det/figures'

#---------------------------
# Three block model
#---------------------------

n <- 300 
a <- .2; b <- .1
B <- matrix(b, nrow = 3, ncol = 3)
diag(B) <- a
d <- 3
#helper functions
get_lp <- function(lam, n){
  #set sparsity
  #alpha_n <- log(n)/n^(1/4)
  
  #get latent positions
  Z <- kronecker(diag(3), rep(1, n/3))
  mat <- B + (lam * rbind(c(0, 0, 0),
  						  c(0, 0, a-b),
                          c(0, a-b, 0)))
  
  tmp <- eigen(mat)
  L <- tmp$vectors[, 1:d] %*% diag(sqrt(tmp$values[1:d]))
  return(Z %*% L)
}
get_mc <- function(x, true){
  
  n <- length(x)
  mc <- numeric(6)
  
  #no permutation
  mc[1] <- sum(x != true) / n 
  
  #flip 1 and 2
  y <- x
  y[x == 1] <- 2; y[x == 2] <- 1
  mc[2] <- sum(y != true) / n
  
  #flip 1 and 3
  y <- x
  y[x == 1] <- 3; y[x == 3] <- 1
  mc[3] <- sum(y != true) / n
  
  #flip 2 and 3
  y <- x
  y[x == 2] <- 3; y[x == 3] <- 2
  mc[4] <- sum(y != true) / n
  
  #1 -> 2, 2-> 3, 3 -> 1
  y <- x
  y[x == 1] <- 2; y[x == 2] <- 3; y[x == 3] <- 1;
  mc[5] <- sum(y != true) / n
  
  #1 -> 3, 2-> 1, 3 -> 2
  y <- x
  y[x == 1] <- 3; y[x == 2] <- 1; y[x == 3] <- 2;
  mc[6] <- sum(y != true) / n
  
  return(min(mc))
}

#set up base latent positions
X <- get_lp(0, n)
P_list <- list()
P_list[[1]] <- tcrossprod(X)
labels <- rep(1:3, each = n/3)

#set up storage 
lambda <- seq(0, 1, length.out = 11)
no_iters <- 1000
mat <- matrix(NA, nrow = length(lambda)*no_iters, ncol = 8)
colnames(mat) <- c('iter_no', 'lambda', 'ASE1', 'ASE2', 'ASEbar', 'AASE', 'OMNI', 'MASE', 'RASE')

count1 <- 1

set.seed(1985)
#start simulation 
for(i in 1:length(lambda)){
  #get latent positions & update P_list
  Y <- get_lp(lambda[i],  n)
  P_list[[2]] <- tcrossprod(Y)
  
  for(j in 1:no_iters){
    
    #---------------------------------
    #     Sampling + Embedding
    #---------------------------------
    
    #sample networks
    A_list <- lapply(P_list, sampP)
    
    #store latent positions in order: ASE1, ASE2, ASEBAR, ABAR, OMNI, MASE
    #embedd into R^2
    X_estimates <- lapply(A_list, function(A) ase(A, d)) #ASE1, ASE2
    
    X1.rot <- procrustes(X_estimates[[1]], X_estimates[[2]])$X.new
    X_estimates[[3]] <- 0.5 * (X1.rot +  X_estimates[[2]]) #ASEBar
    
    X_estimates[[4]] <- ase(Reduce('+', A_list)/2, d) #Abar
    
    X_omni <- omni(A_list, d)
    X_estimates[[5]] <- 0.5 * (X_omni[[1]] +  X_omni[[2]]) #Omnibar
    
    X_mase <- mase(A_list, c(d,d), d)
    X1.rot <- procrustes(X_mase[[1]], X_mase[[2]])$X.new
    X_estimates[[6]] <- 0.5 * (X1.rot +  X_mase[[2]])
    
    #---------------------------------
    #     Clustering 
    #---------------------------------
    
    #cluster & get mc rate
    clusters  <- lapply(X_estimates, function(X) kmeans(X, centers = 3)$cluster)
    mc_rates <- unlist(lapply(clusters, function(x) get_mc(x, labels)))
    
    #store 
    mat[count1, ] <- c(j, lambda[i], mc_rates)
    count1 <- count1 + 1
    
  }
  print(lambda[i])
  
}

#------------------------------------
#     plotting results
#------------------------------------
df <- as.data.frame(mat) %>% 
  melt(id.vars = 1:2, value.name = 'mc_rate', variable.name = 'Method') %>% 
  group_by(lambda, Method) %>%
  summarize(MC_Rate = mean(mc_rate),
            MC_sd = sd(mc_rate)/sqrt(n()))

p <- ggplot(df, aes(lambda, MC_Rate, col = Method, group = Method))+
  geom_ribbon(aes(ymin = MC_Rate - MC_sd, ymax = MC_Rate + MC_sd), alpha = 0.2, col = 'lightgrey')+
  geom_point(alpha = .5)+geom_line()+
  scale_y_log10()+ theme_bw()+
  labs(x = expression(lambda), y = expression('log'[10]*'(Misclass. Rate)'))

p
scale <- 2
# ggsave('3_group_sbm_mc_rate.pdf', 
#        path = fig_path, 
#        device = 'pdf',
#        width = scale * 4, height = scale * 3)
# 
# df <- as.data.frame(dist_mat) %>% 
#   melt(id.vars = c(1:2, 9), value.name = 'dist', variable.name = 'Method') %>% 
#   group_by(lambda, Method, Metric) %>%
#   summarize(Distance = mean(dist),
#             Distance_sd = 1.96*sd(dist)/sqrt(n())) %>%
#   mutate(Metric = ifelse(Metric == 0, 'Euclidean', ifelse(Metric == 1, 'Mahalanobis', 'Cosine')))
# 
# ggplot(df, aes(lambda, Distance, col = Method, group = Method))+
#   geom_ribbon(aes(ymin = Distance - Distance_sd, ymax = Distance + Distance_sd), alpha = 0.2, col = 'lightgrey')+
#   geom_point(alpha = .5)+geom_line()+
#   facet_wrap(~Metric, scales = 'free')+
#   theme_bw()+
#   labs(x = expression(lambda), y = 'Ave. Distance b/t Centroids')
# 
# ggsave('3_group_sbm_distance.pdf', 
#        path = fig_path, 
#        device = 'pdf',
#        width = scale * 4, height = scale * 3)

#--------------------------------------------
#     plot networks and adjacency matrices
#--------------------------------------------
pacman::p_load(igraph, gridExtra, reshape2)
lam_plots <- c(0, 0.5, 1)
p_list <- list()
set.seed(1985)
for(i in 1:length(lam_plots)){
  #get P 
  p_base <- melt(tcrossprod(get_lp(lam_plots[i], 3)), 
                      varnames = c('x', 'y')) %>%
    mutate(value = factor(value, levels = c(0.1, 0.15, 0.2))) %>%
    ggplot(aes(x = x, y = -y, 
               fill = value)) + 
    geom_tile() + 
    theme_void() +
    #theme_classic() +
    labs(y = '', title = bquote(lambda*'='*.(lam_plots[i]))) + 
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5))
  
  if(i == 2){
    p_list[[i]] <- p_base + scale_fill_manual(values = c('lightgrey', 'darkgrey', 'black')) 
  } 
  else{
    p_list[[i]] <- p_base + scale_fill_manual(values = c('lightgrey', 'black'))  
  }
  
  
  #get A 
  #P2 <- tcrossprod(get_lp(0, n))
  #p_list[[i+length(lam_plots)]] <- image(Matrix(P2), xlab = '', ylab = '', sub = '', axes = FALSE)
}

#make figure
ind <- length(lam_plots)+1
p_list[[ind]] <- p + 
  labs(#title = 'Misclassification Simulation', 
       x = '') + 
  theme(#plot.margin=unit(c(1,1,-1,1),"cm"), 
        legend.position = 'top'
        )
lay <- rbind(rep(ind, length(lam_plots) + 2),
             rep(ind, length(lam_plots) + 2),
             rep(ind, length(lam_plots) + 2), 
             c(1, NA, 2, NA, 3)
             )
scale1 <- 125
jpeg(file = paste0(fig_path,'/mis_class_sim.jpeg'), 
     width = ncol(lay) * scale1, 
     height = nrow(lay) * scale1
     )
grid.arrange(grobs = p_list, layout_matrix = lay)
dev.off()

scale2 <- 1.5
pdf(file = paste0(fig_path,'/mis_class_sim.pdf'), 
     width = scale2 * ncol(lay), 
     height = scale2 * nrow(lay))
grid.arrange(grobs = p_list, layout_matrix = lay)
dev.off()






