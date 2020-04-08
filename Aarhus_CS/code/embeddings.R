#------------------------------------------------
#
#       Aarhus CS Multiplex            
#
# This multiplex social network consists of 
# five kinds of online and offline relationships 
# (Facebook, Leisure, Work, Co-authorship, Lunch) 
# among 61 employees of the Computer Science 
# Department at Aarhus.
#------------------------------------------------


#load packages and source basic functions
pacman::p_load('igraph', 'RColorBrewer', 'car', 'ggplot2', 'MCMCpack', 'gplots', 'plot.matrix')

#set directory
setwd('~/Documents/Work/github/BJSE/Aarhus_CS/')
source('./code/basic_functions.R')
source('./code/getElbows.R')

#read in data
layer_ids <- read.table('./CS-Aarhus_Multiplex_Social/Dataset/CS-Aarhus_layers.txt', header = TRUE)
node_ids <- read.table('./CS-Aarhus_Multiplex_Social/Dataset/CS-Aarhus_nodes.txt', header = TRUE)

#captilize layer Ids
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

layer_ids$layerLabel <- sapply(layer_ids$layerLabel, function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2)))

#network colors
colors <- c('orange', 'blue', 'green', 'purple', 'red')

#read in adjacency matrices
A1 <- as.matrix(read.table('./adj_mats/A1.txt'))
A2 <- as.matrix(read.table('./adj_mats/A2.txt'))
A3 <- as.matrix(read.table('./adj_mats/A3.txt'))
A4 <- as.matrix(read.table('./adj_mats/A4.txt'))
A5 <- as.matrix(read.table('./adj_mats/A5.txt'))

#make adjacency list
A_list <- list(A1, A2, A3, A4, A5)
Atil <- make_omni(A_list)

#skree plot 
evals <- eigen(Atil)$values 
plot(evals)

#determine d - get elbows
non_zero_evals <- evals[which(evals > 0)]
elbows <- getElbows(non_zero_evals, n = 3)

#set d 
d <- elbows[1]
Lhat <- ase(Atil, d)

#Goodness of Fit Statistics
T_hat <- matrix(NA, nrow = nrow(A1), ncol = length(A_list) * d)
m <- length(A_list)
n <- nrow(A1)

for(i in 1:length(A_list)){
  T_hat[,(d*(i-1) + 1):(d*i)] <- Lhat[(n*(i-1) + 1):(n*i), ]
}

T_hat_SVD <- svd(T_hat)

#singular values show decrease after d
plot(T_hat_SVD$d, col = c(rep('red', d), rep('black', d*(m-1))))

#right singular vectors should be close to diagnol
V <- T_hat_SVD$v[, 1:d]
heatmap.2(V, Rowv=FALSE, Colv=FALSE, dendrogram='none') 
plot(V, key=list(side=3, cex.axis=0.75), main = '')

V.rot <- procrustes(V, rbind(diag(d), matrix(0, nrow = d*(m-1), ncol = d)))$X.new
heatmap.2(V.rot, Rowv=FALSE, Colv=FALSE, dendrogram='none') 
plot(V.rot, key=list(side=3, cex.axis=0.75), main = '')
plot(V.rot[1:5,6:10], key=list(side=3, cex.axis=0.75), main = '')

#Eigen structure of 11^T \otimes I
k <- 3
mat <- kronecker(matrix(1,nrow = k , ncol = k), diag(k))
e_mat <- eigen(mat.rot)$vectors
plot(e_mat, breaks = seq(-1, 1, by = 0.1), key=list(side=3, cex.axis=0.75), main = '')


#plot first two dimensions

#look at eigenvectors of Atil
pairs(eigen(Atil)$vectors[, 1:d], cex = 0.1)
pairs(Lhat, cex = 0.1)

#Dimension 1 is some type of MMSBM? 

#Dimensions 2 - 5 exhibit strong eigenspoke behavior

#Dimensions 6 - 8 Lesser extent eigenspoke behavior 

#which nodes vary the most / least
barplot(apply(T_hat, 1, function(x) sqrt(sum((x - mean(x))^2))))




plot(Lhat[,1], Lhat[,2], col = rep(colors, each = nrow(A1)))
plot(Lhat[,1], Lhat[,2], col = rep(1:61, 5))

plotdf <- data.frame(X1 = Lhat[,1], X2 = Lhat[,2], graph = as.factor(rep(1:length(A_list), each = nrow(A1))), node_id = rep(1:nrow(A1), length(A_list)))

ggplot(plotdf, aes(X1, X2, col = graph)) + geom_point(alpha = 0.5) + theme_bw()

scatter3d(Lhat[,1], Lhat[,2], Lhat[,3], surface = FALSE)
scatter3d(Lhat[,2], Lhat[,3], Lhat[,4], surface = FALSE)

#project onto sphere
Yhat <- diag(1/sqrt(diag(tcrossprod(Lhat)))) %*% Lhat
scatter3d(Yhat[,1], Yhat[,2], Yhat[,3], surface = FALSE)

#calculate Xbar
Xsum <- 0
for(i in 1:5){
	Xsum <- Lhat[(n*(i-1) + 1):(n*i),]
} 
Xbar <- Xsum/5
pairs(Xbar, cex = 0.1)


Xbarhat <- diag(1/sqrt(diag(tcrossprod(Xbar)))) %*% Xbar
plot(Xbarhat[,1], Xbarhat[,2])

