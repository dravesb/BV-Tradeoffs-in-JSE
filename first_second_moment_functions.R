#----------------------------------------------
#
#         Bias _ Variance Functions 
#          for Omnibus 
#
#----------------------------------------------

#------------------------------
#       Bias Functions
#
# returns the bias matrix
# Xhat - X \in (nm x d)
#-----------------------------

ase_bias <- function(X, C_list){
  #get n and m
  n <- nrow(X)
  m <- length(C_list)
  d <- nrow(C_list[[1]])
  
  #return bias
  return(matrix(0, nrow = n*m, ncol = d))
  
}
abar_bias <- function(X, C_list){
  #calculate Cbar
  m <- length(C_list)
  Cbar <- Reduce('+', C_list)/m
  
  #create bias matrix
  return(Reduce(rbind, lapply(C_list, function(C)  X %*% (sqrt(Cbar) - sqrt(C)))))
}
omni_bias <- function(X, C_list){
  #fetch n,m, d
  n <- nrow(X)
  m <- length(C_list)
  d <- ncol(X)
  
  #V (d x m)  matrix with v1, v2, ..., vd in the __rows__ of V
  V <- Reduce('cbind', lapply(C_list, diag))
  
  #Alpha (m x d)  matrix with alpha^1, alpha^2, ..., alpha^d in the __columns__ of Alpha
  Alpha <- apply(V, 1, function(v){
    H <- 0.5* (tcrossprod(rep(1, m), v) + tcrossprod(v, rep(1, m)))
    eigen.system <- eigen(H)
    u <- eigen.system$vectors[, 1]
    u <- u/sqrt(sum(u^2))
    s <- eigen.system$values[1]
    s_half <- sign(s) * sqrt(abs(s))
    alpha <- abs(u * s_half) #chosen to have positive values 
    return(alpha)
  })
  
  #Set up S matrices by taking the diagnol of the __rows__ of Alpha
  S_list <- list()
  for(g in 1:m) S_list[[g]] <- diag(Alpha[g, ], nrow = d, ncol = d)

  #create bias matrix
  return(Reduce(rbind, lapply(1:m, function(g)  X %*% (S_list[[g]] - sqrt(C_list[[g]])))))
}
omnibar_bias <- function(X, C_list){
  #fetch n,m, d
  n <- nrow(X)
  m <- length(C_list)
  d <- ncol(X)
  
  #V (d x m)  matrix with v1, v2, ..., vd in the __rows__ of V
  V <- Reduce('cbind', lapply(C_list, diag))
  
  #Alpha (m x d)  matrix with alpha^1, alpha^2, ..., alpha^d in the __columns__ of Alpha
  Alpha <- apply(V, 1, function(v){
    H <- 0.5* (tcrossprod(rep(1, m), v) + tcrossprod(v, rep(1, m)))
    eigen.system <- eigen(H)
    u <- eigen.system$vectors[,1]
    u <- u/sqrt(sum(u^2))
    s <- eigen.system$values[1]
    s_half <- sign(s) * sqrt(abs(s))
    alpha <- abs(u * s_half) #chosen to have positive values 
    return(alpha)
  })
  
  #Set up Sbar matrix by taking diagnol of the __column means__ of Alpha
  Sbar <- diag(colMeans(Alpha), nrow = d, ncol = d)
  
  #create bias matrix
  return(Reduce(rbind, lapply(1:m, function(g)  X %*% (Sbar - sqrt(C_list[[g]])))))
}

#------------------------------
#       Variance Functions
#
#     For ER Example
# 
#-----------------------------

ase_variance <- function(X, C_list, Sigma_tilde, Delta, graph){
  #set p
  p <- 0.5
  
  #fetch n, m, d
  n <- nrow(X)
  
  #return variance
  return(diag(1/diag(C_list[[graph]] * Delta), nrow = 1, ncol = 1) %*% Sigma_tilde(sqrt(p), C_list[[graph]]) %*% diag(1/diag(C_list[[graph]] * Delta), nrow = 1, ncol = 1)/n)
}
abar_variance <- function(X, C_list, Sigma_tilde, Delta, graph){
  #set p
  p <- 0.5
  
  #fetch n and Cbar
  n <- nrow(X)
  m <- length(C_list)
  Cbar <- Reduce('+', C_list) / m 
  
  #return variance
  return(diag(1/diag(Cbar * Delta), nrow = 1, ncol = 1) %*% Sigma_tilde(sqrt(p), C_list[[graph]]) %*% diag(1/diag(Cbar * Delta), nrow = 1, ncol = 1)/n)
}
omni_variance <- function(X, C_list, Sigma_tilde, Delta, graph){
  #set p 
  p <- 0.5
  
  #fetch n, m
  n <- nrow(X)
  m <- length(C_list)
  d <- ncol(X)
  
  #get S_list
  #V (d x m)  matrix with v1, v2, ..., vd in the __rows__ of V
  V <- Reduce('cbind', lapply(C_list, diag))
  
  #Alpha (m x d)  matrix with alpha^1, alpha^2, ..., alpha^d in the __columns__ of Alpha
  Alpha <- apply(V, 1, function(v){
    H <- 0.5* (tcrossprod(rep(1, m), v) + tcrossprod(v, rep(1, m)))
    eigen.system <- eigen(H)
    u <- eigen.system$vectors[,1]
    u <- u/sqrt(sum(u^2))
    s <- eigen.system$values[1]
    s_half <- sign(s) * sqrt(abs(s))
    alpha <- abs(u * s_half) #chosen to have positive values 
    return(alpha)
  })
  
  #Set up S matrices by taking the diagnol of the __rows__ of Alpha
  S_list <- list()
  for(g in 1:m) S_list[[g]] <- diag(Alpha[g, ], nrow = d, ncol = d)
  
  #set pre and post matrix
  pre_post <- 2 * Reduce('+', lapply(S_list, function(x) x^2)) * Delta
  
  #calculate main main variance term
  V1 <- (S_list[[graph]] + m * Reduce('+', S_list)) %*% Sigma_tilde(sqrt(p), C_list[[graph]]) %*% (S_list[[graph]] + m * Reduce('+', S_list))
  V2 <- matrix(0, nrow = d, ncol = d)
  for(k in (1:m)[-graph]){
    V2 <- V2 + S_list[[k]] %*% Sigma_tilde(sqrt(p), C_list[[k]]) %*% S_list[[k]] 
  }
  
  #return final variance term
  return(diag(1/diag(pre_post), nrow = 1, ncol = 1) %*% (V1 + V2) %*% diag(1/diag(pre_post), nrow = 1, ncol = 1)/n)
  
  
}
omnibar_variance <- function(X, C_list, Sigma_tilde, Delta, graph){
  
  #set p 
  p <- 0.5
  
  #fetch n, m
  n <- nrow(X)
  m <- length(C_list)
  d <- ncol(X)
  
  #get S_list
  #V (d x m)  matrix with v1, v2, ..., vd in the __rows__ of V
  V <- Reduce('cbind', lapply(C_list, diag))
  
  #Alpha (m x d)  matrix with alpha^1, alpha^2, ..., alpha^d in the __columns__ of Alpha
  Alpha <- apply(V, 1, function(v){
    H <- 0.5* (tcrossprod(rep(1, m), v) + tcrossprod(v, rep(1, m)))
    eigen.system <- eigen(H)
    u <- eigen.system$vectors[, 1]
    u <- u/sqrt(sum(u^2))
    s <- eigen.system$values[1]
    s_half <- sign(s) * sqrt(abs(s))
    alpha <- abs(u * s_half) #chosen to have positive values 
    return(alpha)
  })
  
  #Set up S matrices by taking the diagnol of the __rows__ of Alpha
  S_list <- list()
  for(g in 1:m) S_list[[g]] <- diag(Alpha[g, ], nrow = d, ncol = d)
  
  #set pre and post matrix
  pre_post <- 2 * Reduce('+', lapply(S_list, function(x) x^2)) * Delta
  
  #calculate main main variance term
  V1 <- matrix(0, nrow = d, ncol = d)
  for(k in 1:m){
    V1 <- V1 + (m*S_list[[k]] + Reduce('+', S_list)/m) %*% Sigma_tilde(sqrt(p), C_list[[k]]) %*% (m*S_list[[k]] + Reduce('+', S_list)/m)
  }
  
  #return final variance term
  return(diag(1/diag(pre_post), nrow = 1, ncol = 1) %*% V1 %*% diag(1/diag(pre_post), nrow = 1, ncol = 1)/n)
  
  
}

#------------------------------
#  Theoretical Bias Figures
#-----------------------------

#Erdos Reyni Example
p <- 0.5
X <- matrix(sqrt(p), ncol = 1, nrow = 1)
C <- seq(0, 1/sqrt(p), length.out = 25)

plotdf <- matrix(NA, nrow = 2 * 4 * length(C), ncol = 4)
colnames(plotdf) <- c("Method", "c", 'graph', 'bias2')
start <- 1; stop <- 8

for(i in 1:length(C)){
  #make C list
  C_list <- list()
  C_list[[1]] <- diag(1, nrow = 1, ncol = 1)
  C_list[[2]] <- diag(C[i], nrow = 1, ncol = 1)
  
  #fetch bias
  bias.here <- rbind(ase_bias(X, C_list),
                     abar_bias(X, C_list), 
                     omni_bias(X, C_list), 
                     omnibar_bias(X, C_list))
  
  #store
  plotdf[start:stop,] <- cbind(rep(c('ASE', 'Abar', 'Omni', 'Omnibar'), each = 2),
                               rep(C[i], 8),
                               rep(c('Graph 1', 'Graph 2'), 4),
                               bias.here^2)
  start <- stop + 1
  stop <- stop + 8
  
}

#cast as dataframe 
biasdf <- data.frame(Method = plotdf[,1],
                     c = as.numeric(plotdf[,2]),
                     graph = plotdf[,3],
                     bias2 = abs(as.numeric(plotdf[,4])))

ggplot(biasdf, aes(c, bias2, col = Method))+
  geom_point(alpha = .5)+
  geom_line(size = .5)+
  facet_grid(cols = vars(graph))+
  theme_bw()+
  scale_y_sqrt()+
  labs(y = expression('Bias'^2),
       x = "c")   


#-------------------------------
#  Theoretical Variance Figures
#-------------------------------

#Erdos Reyni Example
p <- 0.5
X <- matrix(sqrt(p), ncol = 1, nrow = 100)
C <- seq(0, 1/sqrt(p), length.out = 25)
Delta <- matrix(p, nrow = 1, ncol = 1)
Sigma_tilde <- function(sqrt_p, c){
  matrix((sqrt_p^4*c) - (sqrt_p^6*c^2), nrow = 1, ncol = 1)
}

plotdf <- matrix(NA, nrow = 2 * 4 * length(C), ncol = 4)
colnames(plotdf) <- c("method", "c", 'graph', 'variance')
start <- 1; stop <- 8

for(i in 1:length(C)){
  #make C list
  C_list <- list()
  C_list[[1]] <- diag(1, nrow = 1, ncol = 1)
  C_list[[2]] <- diag(C[i], nrow = 1, ncol = 1)
  
  #fetch variance graph1 
  var.here <- rbind(ase_variance(X, C_list, Sigma_tilde, Delta, graph = 1),
                    ase_variance(X, C_list, Sigma_tilde, Delta, graph = 2),
                    abar_variance(X, C_list, Sigma_tilde, Delta, graph = 1), 
                    abar_variance(X, C_list, Sigma_tilde, Delta, graph = 2), 
                    omni_variance(X, C_list, Sigma_tilde, Delta, graph = 1), 
                    omni_variance(X, C_list, Sigma_tilde, Delta, graph = 2), 
                    omnibar_variance(X, C_list, Sigma_tilde, Delta, graph = 1),
                    omnibar_variance(X, C_list, Sigma_tilde, Delta, graph = 2)
                    )
    
    #store
    plotdf[start:stop,] <- cbind(rep(c('ASE', 'Abar', 'Omni', 'Omnibar'), each = 2),
                                 rep(C[i], 8),
                                 rep(c('Graph 1', 'Graph 2'), 4),
                                 var.here)
    start <- stop + 1
    stop <- stop + 8  
  
}

#cast as dataframe 
vardf <- data.frame(Method = plotdf[,1],
                     c = as.numeric(plotdf[,2]),
                     graph = plotdf[,3],
                     var = as.numeric(plotdf[,4]))

ggplot(vardf, aes(c, var, col = Method))+
  geom_point(alpha = .5)+
  geom_line(size = .5)+
  facet_grid(cols = vars(graph))+
  theme_bw()+
  scale_y_sqrt()+
  labs(y = 'Variance',
       x = "c")   


#-------------------------------
#  Theoretical MSE Figures
#-------------------------------
msedf <- data.frame(biasdf[, c('Method', 'c', 'graph', 'bias2')], 
                    vardf[, 'var'],
                    MSE = biasdf[, 'bias2'] + vardf[, 'var']
)


ggplot(msedf, aes(c, MSE, col = Method))+
  geom_point(alpha = .5)+
  geom_line(size = .5)+
  facet_grid(cols = vars(graph))+
  theme_bw()+
  scale_y_log10()+
  labs(y = 'MSE',
       x = "c")   














