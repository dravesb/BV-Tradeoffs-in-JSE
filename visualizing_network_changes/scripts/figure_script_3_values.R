#--------------------------------------
#
#       Plotting Networks for 
#        different t values
#          (3 values)    
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

#set directory
setwd('~/Documents/Work/github/BJSE/visualizing_network_changes/scripts/')

#-----------------------------
#     Set up Base Model
#-----------------------------

#set up blocks
B <- matrix(c(.25, .05, .05, .25), byrow = T, nrow = 2)
b_ase <- ase(B, 2)

#get latent positions
x1 <- b_ase[1,]; x2 <- b_ase[2,]

#set prior probabilities
pi <- .5

#get rotation (eigenvectors of Delta)
Delta <- pi * tcrossprod(x1) + (1 - pi)*tcrossprod(x2) 
R <- eigen(Delta)$vectors

#Apply rotation to x1 and x2
x1_til <- t(R) %*% x1 
x2_til <- t(R) %*% x2

#set base latent positions
L <- rbind(t(x1_til), t(x2_til))

#-----------------------------
#     Set up Prob matrices
#-----------------------------

#converges to ER with p = .3 (1,1) --> (2,0)
C <- function(t){
  diag(c(1-t/3, -t + 1))
}

#set up latent positions
n <- 50
com <- sample(c(1, 2), size = n, prob = c(pi, 1-pi), replace = TRUE)
X <- L[com,]

#set up P matrices
P1 <- tcrossprod(X, X %*% C(0)) # t = 0 
P2 <- tcrossprod(X, X %*% C(.5)) # t = 0.5
P3 <- tcrossprod(X, X %*% C(1)) # t = 1

#-----------------------------
#     Sample Adjmatrices
#-----------------------------

set.seed(1985)
#sample adjacency matrices
A1 <- sampP(P1)
A2 <- sampP(P2)
A3 <- sampP(P3)

#-----------------------------
#     Visualize
#-----------------------------
library(igraph)

#set up graph
g1 <- graph_from_adjacency_matrix(A1, mode = "undirected")
g2 <- graph_from_adjacency_matrix(A2, mode = "undirected")
g3 <- graph_from_adjacency_matrix(A3, mode = "undirected")

#change node color 
#cols <- c("#ff000088","#0000ff88")
cols <- c("red","blue")
V(g1)$color <- cols[com]
V(g2)$color <- cols[com]
V(g3)$color <- cols[com]

#take off labels
V(g1)$label <- ""
V(g2)$label <- ""
V(g3)$label <- ""

#change edge width
E(g1)$width <- E(g2)$width <- E(g3)$width <- E(g4)$width <- .1

#change edge color
E(g1)$color <- E(g2)$color <- E(g3)$color <- E(g4)$color <- "grey80"

#change vertex size
vs <- 5

#plot networks
pdf("../figures/g1_three_vals.pdf")
plot(g1, vertex.size = vs)
dev.off()

pdf("../figures/g2_three_vals.pdf")
plot(g2, vertex.size = vs)
dev.off()

pdf("../figures/g3_three_vals.pdf")
plot(g3, vertex.size = vs)
dev.off()

#-----------------------------
#     Visualize Change in 
#     Latent Positions
#-----------------------------

setwd("~/Documents/Work/github/BJSE/visualizing_network_changes/")

#set latent positions
X0 <- L %*% C(0) #t = 0
X1 <- L %*% C(0.5) #t = 0.4
X2 <- L %*% C(1) #t = 1

#combine
X_move <- as.data.frame(rbind(X0, X1, X2))
colnames(X_move) <- c("x", "y")
plot(X_move)

#plot together and highlight points for different t values
set.seed(1)
ggplot()+
  geom_point(data = X_move[1:2,], aes(x, y), shape= 23, fill = "white",size = 5)+
  geom_point(data = X_move[1:2,], aes(x, y), col = c("red", "blue"))+
  geom_point(data = X_move[3:4,], aes(x, y), col = rep(c("red", "blue"), 1), alpha = .2, size = 2)+
  geom_jitter(data = X_move[5:6,], aes(x, y), col = c("red", "blue"),
              alpha = .2, size = 2, width = .001, height = 0)+
  labs(x = "", y = "",
       title = )+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("./figures/latent_position_0_three_vals.pdf",device = "pdf", width = 2.8, height = 2.8, units = "in")

set.seed(1)
ggplot()+
  geom_point(data = X_move[3:4,], aes(x, y), shape= 23, fill = "white",size = 5)+
  geom_point(data = X_move[3:4,], aes(x, y), col = c("red", "blue"))+
  geom_point(data = X_move[c(1:2),], aes(x, y), col = rep(c("red", "blue"), 1), alpha = .2, size = 2)+
  geom_jitter(data = X_move[5:6,], aes(x, y), col = c("red", "blue"),
              alpha = .2, size = 2, width = .001, height = 0)+
  labs(x = "", y = "")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("./figures/latent_position_2_three_vals.pdf",device = "pdf", width = 2.8, height = 2.8, units = "in")


set.seed(1)
ggplot()+
  geom_point(data = X_move[5:6,], aes(x, y), shape= 23, fill = "white",size = 5)+
  #geom_point(data = X_move[7:8,], aes(x, y), col = c("red", "blue"))+
  geom_point(data = X_move[1:4,], aes(x, y), col = rep(c("red", "blue"), 2), alpha = .2, size = 2)+
  geom_jitter(data = X_move[5:6,], aes(x, y), col = c("red", "blue"),
              alpha = 1, size = 2, width = .001, height = 0)+
  labs(x = "", y = "")+
  theme_bw()+ 
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("./figures/latent_position_1_three_vals.pdf",device = "pdf",width = 2.8, height = 2.8, units = "in")


#-----------------------------
#     Visualize Omnibus  
#     Embeddings
#-----------------------------
make_omni <- function(mats){
  #H(x) = (1x^T + x1^T)/2
  H <- function(g,m){
    ones <- rep(1, m)
    e <- diag(m)[,g]
    .5 * (tcrossprod(ones,e) + tcrossprod(e,ones))
  }
  
  #sum up each kronecker 
  m <- length(mats)
  Reduce("+", lapply(1:m, function(x) kronecker(H(x,m), mats[[x]])))
}
ase <- function(A,d){
  #function to normalize columns of A
  norm2 <- function(u){
    sqrt(sum(u^2))
  } 
  normalize.cols<- function(A){
    norm.vec <- function(u) u/norm2(u) #define vector normalization func.
    if(ncol(A) == 1) return(norm.vec(A[,1]))
    apply(A, 2, norm.vec) # vectorize 
  } 
  
  #construct ASE 
  E <- eigen(A)
  U <- normalize.cols(as.matrix(E$vectors[,1:d], ncol = d))
  S_sqrt <- diag(x = sign(E$values[1:d]), ncol = d, nrow = d)*diag(sqrt(abs(E$values[1:d])), nrow = d, ncol = d)
  U %*% S_sqrt
}

#make omni and embedd
A_list <- list()
A_list[[1]] <- A1; A_list[[2]] <- A2; A_list[[3]] <- A3
Atil <- make_omni(A_list)
Z <- ase(Atil, d = 2)

#making plotting data frame
plot_df <- data.frame(X = Z[,1], Y = Z[,2], 
                      network = as.factor(rep(1:3, each = n)),
                      community = as.factor(rep(com, 3)))

#plot Embedded points network 1 
ggplot(plot_df %>% filter(network == 1), aes(X, Y, col = community))+
  geom_point(alpha = 0.5)+
  scale_color_manual(values = c('red', 'blue'))+
  theme_bw()+
  theme(legend.position = 'none')
ggsave("./figures/embedded_points_0_three_vals.pdf",device = "pdf",width = 2.8, height = 2.8, units = "in")

#plot Embedded points network 2
ggplot(plot_df %>% filter(network == 2), aes(X, Y, col = community))+
  geom_point(alpha = 0.5)+
  scale_color_manual(values = c('red', 'blue'))+
  theme_bw()+
  theme(legend.position = 'none')
ggsave("./figures/embedded_points_5_three_vals.pdf",device = "pdf",width = 2.8, height = 2.8, units = "in")

#plot Embedded points network 3
ggplot(plot_df %>% filter(network == 3), aes(X, Y, col = community))+
  geom_point(alpha = 0.5)+
  scale_color_manual(values = c('red', 'blue'))+
  theme_bw()+
  theme(legend.position = 'none')
ggsave("./figures/embedded_points_1_three_vals.pdf",device = "pdf",width = 2.8, height = 2.8, units = "in")

#individual embeddings
library(MCMCpack)
Xhat <- lapply(A_list, function(x) ase(x, 2))
Xhat.rot <- lapply(Xhat, function(x) procrustes(x, X)$X.new)

plot_df <- data.frame(X = Xhat.rot[[1]][,1],
                      Y = Xhat.rot[[1]][,2],
                      community = as.factor(com))
ggplot(plot_df, aes(X, Y, col = community))+
  geom_point(alpha = 0.5)+
  scale_color_manual(values = c('red', 'blue'))+
  theme_bw()+
  theme(legend.position = 'none')
ggsave("./figures/ase_0_three_vals.pdf",device = "pdf",width = 2.8, height = 2.8, units = "in")

plot_df <- data.frame(X = Xhat.rot[[2]][,1],
                      Y = Xhat.rot[[2]][,2],
                      community = as.factor(com))
ggplot(plot_df, aes(X, Y, col = community))+
  geom_point(alpha = 0.5)+
  scale_color_manual(values = c('red', 'blue'))+
  theme_bw()+
  theme(legend.position = 'none')
ggsave("./figures/ase_5_three_vals.pdf",device = "pdf",width = 2.8, height = 2.8, units = "in")

plot_df <- data.frame(X = Xhat.rot[[3]][,1],
                      Y = Xhat.rot[[3]][,2],
                      community = as.factor(com))
ggplot(plot_df, aes(X, Y, col = community))+
  geom_point(alpha = 0.5)+
  scale_color_manual(values = c('red', 'blue'))+
  theme_bw()+
  theme(legend.position = 'none')
ggsave("./figures/ase_1_three_vals.pdf",device = "pdf",width = 2.8, height = 2.8, units = "in")

#-----------------------------
#     Visualize Networks 
#     for Beamer
#-----------------------------

#change node color 
#cols <- c("#ff000088","#0000ff88")
cols <- c("red","blue")
V(g1)$color <- cols[com]
V(g2)$color <- cols[com]
V(g3)$color <- cols[com]

#take off labels
V(g1)$label <- ""
V(g2)$label <- ""
V(g3)$label <- ""

#change edge width
E(g1)$width <- E(g2)$width <- E(g3)$width <- E(g4)$width <- .1

#change edge color
E(g1)$color <- E(g2)$color <- E(g3)$color <- E(g4)$color <- "grey40"

#change vertex size
vs <- 8
set.seed(1)
#plot networks
pdf("./figures/g1_beamer.pdf")
plot(g1, vertex.size = vs)
dev.off()

pdf("./figures/g2_beamer.pdf")
plot(g2, vertex.size = vs)
dev.off()

pdf("./figures/g3_beamer.pdf")
plot(g3, vertex.size = vs)
dev.off()


#adjacency matrix visualization
library(Matrix)
perm <- c(which(com == 1), which(com == 2))

pdf("./figures/A1_beamer.pdf")
image(Matrix(A1[perm, perm]), xlab = '', ylab = '', sub = '',
      border.col = NA, col.regions = colorRampPalette(c('white', 'blue'))(30))
dev.off()

pdf("./figures/A2_beamer.pdf")
image(Matrix(A2[perm, perm]), xlab = '', ylab = '', sub = '',
      border.col = NA, col.regions = colorRampPalette(c('white', 'blue'))(30))
dev.off()

pdf("./figures/A3_beamer.pdf")
image(Matrix(A3[perm, perm]), xlab = '', ylab = '', sub = '',
      border.col = NA, col.regions = colorRampPalette(c('white', 'blue'))(30))
dev.off()
