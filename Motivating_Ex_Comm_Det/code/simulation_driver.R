#------------------------------
#
#   Clustering Example 
#    Omni v Abar v ASE
#
#------------------------------

#source libraries + functions
source('/Users/benjamindraves/Desktop/OmniAnalysis/code_base/loadAll.R')

#---------------------------
# Two block model
#---------------------------
#set up model
n <- 100 
a <- .2; b <- .1
B <- matrix(b, nrow = 2, ncol = 2)
diag(B) <- a


#set up model
n <- 100 
a <- .2; b <- .1
B <- rbind(c(a, b),
           c(b, a))
D <- rbind(c(0, a-b),
           c(a-b, 0))

#helper functions
get_lp <- function(a, n){
  Z <- cbind(rep(c(1,0), each = n/2), 
             rep(c(0,1), each = n/2))
  
  tmp <- eigen(B + (a*D))
  L <- tmp$vectors[, 1:2] %*% diag(sqrt(tmp$values[1:2]))
  return(Z %*% L)
}
get_mc <- function(x, true){
  mc1 <- sum(x != true)/length(x)
  mc2 <- sum(ifelse(x == 1, 2, 1) != true)/length(x)
  return(min(mc1, mc2))
}

#set up base latent positions
X <- get_lp(0, n)
P_list <- list()
P_list[[1]] <- tcrossprod(X)
labels <- rep(1:2, each = n/2)

#set up storage 
alpha <- seq(0, 1, length.out = 20)
no_iters <- 1000
mat <- matrix(NA, nrow = length(alpha)*no_iters, ncol = 8)
colnames(mat) <- c('iter_no', 'a', 'ASE1', 'ASE2', 'ASEbar', 'Abar', 'Omni', 'MASE')
count <- 1

#start simulation 
for(i in 1:length(alpha)){
  #get latent positions & update P_list
  Y <- get_lp(alpha[i],  n)
  P_list[[2]] <- tcrossprod(Y)
  
  for(j in 1:no_iters){
    #sample networks
    A_list <- lapply(P_list, sampP)
    
    #embedd into R^2
    X_omni <- omni(A_list, d)
    X_mase <- mase(A_list, c(d,d), d)
    X_bar <- ase(Reduce('+', A_list)/2, d)
    X_ase <- lapply(A_list, function(A) ase(A, d))
    
    #Averages
    X_omni_bar <- 0.5 * (X_omni[[1]] +  X_omni[[2]])
    X1.rot <- procrustes(X_mase[[1]], X_mase[[2]])$X.new
    X_mase_bar <- 0.5 * (X1.rot +  X_mase[[2]])
    X1.rot <- procrustes(X_ase[[1]], X_ase[[2]])$X.new
    X_ase_bar <- 0.5 * (X1.rot +  X_ase[[2]])
    
    #clacualte missclassfication rates
    ase1 <- get_mc(kmeans(X_ase[[1]], centers = 2)$cluster, labels)
    ase2 <- get_mc(kmeans(X_ase[[2]], centers = 2)$cluster, labels)
    asebar <- get_mc(kmeans(X_ase_bar, centers = 2)$cluster, labels)
    abar <- get_mc(kmeans(X_bar, centers = 2)$cluster, labels)
    omni_bar <- get_mc(kmeans(X_omni_bar, centers = 2)$cluster, labels)
    mase_bar <- get_mc(kmeans(X_mase_bar, centers = 2)$cluster, labels)
    
    #store 
    mat[count, ] <- c(j, alpha[i], ase1, ase2, asebar, abar, omni_bar, mase_bar)
    count <- count + 1
  }
  print(alpha[i])
  
}

#------------------------------------
#     plotting results
#------------------------------------
df <- as.data.frame(mat) %>% 
  melt(id.vars = 1:2, value.name = 'mc_rate', variable.name = 'Method') %>% 
  group_by(a, Method) %>%
  summarize(MC_Rate = mean(mc_rate),
            MC_sd = 1.96*sd(mc_rate)/sqrt(n()))

ggplot(df, aes(a, MC_Rate, col = Method, group = Method))+
  geom_ribbon(aes(ymin = MC_Rate - MC_sd, ymax = MC_Rate + MC_sd), alpha = 0.2, col = 'lightgrey')+
  geom_point(alpha = .1)+geom_line()+
  theme_bw()+
  labs(x = expression(alpha), y = 'Misclass. Rate')


#--------------------------------------------------------------------------------------------------------------

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
  #get laten positions
  Z <- kronecker(diag(3), rep(1, n/3))
  if(lam < 0){
    mat <- B + (lam * rbind(c(0, b-a, 0),
                            c(b-a, 0, 0),
                            c(0, 0, 0))) 
  }else{
    mat <- B + (lam * rbind(c(0, 0, 0),
                            c(0, 0, a-b),
                            c(0, a-b, 0)))
  }
  
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

#distance functions
d_euclidean <- function(X, clust){
  #get centers
  c1 <- colMeans(X[clust==1, ])
  c2 <- colMeans(X[clust==2, ])
  c3 <- colMeans(X[clust==3, ])
  
  #get pairwise distances
  d <- dist(rbind(c1, c2, c3))
  
  #return average
  return(mean(d[lower.tri(d, diag = TRUE)], na.rm = TRUE))
  
}

get_mahalanobis <- function(X, Y){
  #get sample sizes
  nx <- nrow(X)
  ny <- nrow(Y)
  
  #get covariances
  CX <- cov(X)
  CY <- cov(Y)
  
  #pool covariances
  Cpool <- ((nx - 1)*CX + (ny - 1)*CY)/(nx + ny -2)
  
  #get column means
  xbar <- colMeans(X)
  ybar <- colMeans(Y)
  
  #calcualte R^{-T}(x - y)
  Rc <- chol(Cpool)
  fact <- forwardsolve(t(Rc), xbar - ybar)
  
  #return distance
  return(sqrt(sum(fact^2)))
  
}
d_mahalanobis <- function(X, clust){
  nk <- nrow(X) / 3
  d1 <- get_mahalanobis(X[clust==1, ], X[clust==2, ])
  d2 <- get_mahalanobis(X[clust==1, ], X[clust==3, ])
  d3 <- get_mahalanobis(X[clust==2, ], X[clust==3, ])
  return(mean(c(d1, d2, d3)))
}

cos_dist <- function(x, y){
  as.numeric(sqrt(2*(1 - crossprod(x, y) / sqrt(crossprod(x) * crossprod(y)))))
}
d_cosine <- function(X, clust){
  #get centers
  c1 <- colMeans(X[clust==1, ])
  c2 <- colMeans(X[clust==2, ])
  c3 <- colMeans(X[clust==3, ])
    
  #get pairwise distances
  d1 <- cos_dist(c1, c2)
  d2 <- cos_dist(c1, c3)
  d3 <- cos_dist(c2, c3)
    
  #return average
  return(mean(c(d1, d2, d3)))
    
}
  

#set up base latent positions
X <- get_lp(0, n)
P_list <- list()
P_list[[1]] <- tcrossprod(X)
labels <- rep(1:3, each = n/3)

#set up storage 
lambda <- seq(-1, 1, length.out = 21)
no_iters <- 1000
mat <- matrix(NA, nrow = length(lambda)*no_iters, ncol = 8)
colnames(mat) <- c('iter_no', 'lambda', 'ASE1', 'ASE2', 'ASEbar', 'Abar', 'Omni', 'MASE')

dist_mat <- matrix(NA, nrow = 3*length(lambda)*no_iters, ncol = 9)
colnames(dist_mat) <- c('iter_no', 'lambda', 'ASE1', 'ASE2', 'ASEbar', 'Abar', 'Omni', 'MASE', 'Metric') 

count1 <- 1
count2 <- 1

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
    X_estimates[[5]] <- 0.5 * (X_omni[[1]] +  X_omni[[2]]) #omnibar
    
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
    
    #---------------------------------
    #     Distance between centroids
    #---------------------------------
    
    #Distances
    euclidean_dist <- sapply(1:length(X_estimates), function(i) d_euclidean(X_estimates[[i]], clusters[[i]]))
    cosine_dist <- sapply(1:length(X_estimates), function(i) d_cosine(X_estimates[[i]], clusters[[i]]))
    mahalanobis_dist <- sapply(1:length(X_estimates), function(i) d_mahalanobis(X_estimates[[i]], clusters[[i]]))
    
    dist_mat[count2, ] <- c(j, lambda[i], euclidean_dist, 0)
    dist_mat[count2+1, ] <- c(j, lambda[i], cosine_dist, 1)
    dist_mat[count2+2, ] <- c(j, lambda[i], mahalanobis_dist, 2)
    count2 <- count2 + 3
    
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
            MC_sd = 1.96*sd(mc_rate)/sqrt(n()))

ggplot(df, aes(lambda, MC_Rate, col = Method, group = Method))+
  geom_ribbon(aes(ymin = MC_Rate - MC_sd, ymax = MC_Rate + MC_sd), alpha = 0.2, col = 'lightgrey')+
  geom_point(alpha = .5)+geom_line()+
  scale_y_log10()+ theme_bw()+
  labs(x = expression(lambda), y = expression('log'[10]*'(Misclass. Rate)'))

scale <- 2
ggsave('3_group_sbm_mc_rate.pdf', 
       path = '/Users/benjamindraves/Documents/Work/github/BJSE/Motivating_Ex_Comm_Det/figures', 
       device = 'pdf',
       width = scale * 4, height = scale * 3)

df <- as.data.frame(dist_mat) %>% 
  melt(id.vars = c(1:2, 9), value.name = 'dist', variable.name = 'Method') %>% 
  group_by(lambda, Method, Metric) %>%
  summarize(Distance = mean(dist),
            Distance_sd = 1.96*sd(dist)/sqrt(n())) %>%
  mutate(Metric = ifelse(Metric == 0, 'Euclidean', ifelse(Metric == 1, 'Mahalanobis', 'Cosine')))

ggplot(df, aes(lambda, Distance, col = Method, group = Method))+
  geom_ribbon(aes(ymin = Distance - Distance_sd, ymax = Distance + Distance_sd), alpha = 0.2, col = 'lightgrey')+
  geom_point(alpha = .5)+geom_line()+
  facet_wrap(~Metric, scales = 'free')+
  theme_bw()+
  labs(x = expression(lambda), y = 'Ave. Distance b/t Centroids')

ggsave('3_group_sbm_distance.pdf', 
       path = '/Users/benjamindraves/Documents/Work/github/BJSE/Motivating_Ex_Comm_Det/figures', 
       device = 'pdf',
       width = scale * 4, height = scale * 3)


