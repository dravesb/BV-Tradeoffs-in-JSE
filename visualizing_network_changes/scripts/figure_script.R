#--------------------------------------
#
#       Plotting Networks for 
#        different t values
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

#set directory
setwd('~/Documents/Work/github/BJSE/visualizing_network_changes/scripts/')

#-----------------------------
#     Set up Base Model
#-----------------------------

#set up blocks
B <- matrix(c(.25, .03, .03, .25), byrow = T, nrow = 2)
b_ase <- ase(B, 2)

#get latent positions
x1 <- b_ase[1,]; x2 <- b_ase[2,]

#set prior probabilities
pi <- .5

#get rotation (eigenvectors of Delta)
Delta <- pi * tcrossprod(x1) + (1 - pi)*tcrossprod(x2) 
R <- eigen(Delta)$vectors

#Apply rotation to x1 and x2
x1_til <- x1 
x2_til <- x2

#set base latent positions
L <- rbind(t(x1_til), t(x2_til))
L <- rbind(x1, x2)

#-----------------------------
#     Set up Prob matrices
#-----------------------------

#converges to ER with p = .3 (1,1) --> (2,0)
#C <- function(t){
#  diag(c(t + 1, -t + 1))
#}

#converges to ER with p = .1 (1,1) --> (2/3,0)
C <- function(t){
  diag(c(1-t/3, -t + 1))
}



#set up latent positions
n <- 50
com <- sample(c(1, 2), size = n, prob = c(pi, 1-pi), replace = TRUE)
X <- L[com,]

#set up P matrices
alpha_n <- log(n) / n^(0.25)

P1 <- tcrossprod(X, X %*% C(0)) # t = 0 
P2 <- tcrossprod(X, X %*% C(.3)) # t = 0.3
P3 <- tcrossprod(X, X %*% C(.6)) # t = 0.6
P4 <- tcrossprod(X, X %*% C(1)) # t = 1

#-----------------------------
#     Sample Adjmatrices
#-----------------------------

set.seed(1995)
#sample adjacency matrices
A1 <- sampP(P1)
A2 <- sampP(P2)
A3 <- sampP(P3)
A4 <- sampP(P4)

#-----------------------------
#     Visualize
#-----------------------------
library(igraph)

#set up graph
g1 <- graph_from_adjacency_matrix(A1, mode = "undirected")
g2 <- graph_from_adjacency_matrix(A2, mode = "undirected")
g3 <- graph_from_adjacency_matrix(A3, mode = "undirected")
g4 <- graph_from_adjacency_matrix(A4, mode = "undirected")

#change node color 
#cols <- c("#ff000088","#0000ff88")
cols <- c("red","blue")
V(g1)$color <- cols[com]
V(g2)$color <- cols[com]
V(g3)$color <- cols[com]
V(g4)$color <- cols[com]

#take off labels
V(g1)$label <- ""
V(g2)$label <- ""
V(g3)$label <- ""
V(g4)$label <- ""

#change edge width
E(g1)$width <- E(g2)$width <- E(g3)$width <- E(g4)$width <- .1

#change edge color
E(g1)$color <- E(g2)$color <- E(g3)$color <- E(g4)$color <- "grey80"

#change vertex size
vs <- 5

#plot networks
pdf("../figures/g1.pdf")
plot(g1, vertex.size = vs)
dev.off()
jpeg("../figures/g1.jpeg")
plot(g1, vertex.size = vs)
dev.off()

pdf("../figures/g2.pdf")
plot(g2, vertex.size = vs)
dev.off()
jpeg("../figures/g2.jpeg")
plot(g2, vertex.size = vs)
dev.off()

pdf("../figures/g3.pdf")
plot(g3, vertex.size = vs)
dev.off()
jpeg("../figures/g3.jpeg")
plot(g3, vertex.size = vs)
dev.off()

pdf("../figures/g4.pdf")
plot(g4, vertex.size = vs)
dev.off()
jpeg("../figures/g4.jpeg")
plot(g4, vertex.size = vs)
dev.off()


#-----------------------------
#     Visualize Change in 
#     Latent Positions
#-----------------------------

setwd("~/Documents/Work/github/BJSE/visualizing_network_changes/")

#set latent positions
X0 <- L %*% sqrt(C(0)) #t = 0
X1 <- L %*% sqrt(C(.3)) #t = 0.6
X2 <- L %*% sqrt(C(.6)) #t = 0.6
X3 <- L %*% sqrt(C(1)) #t = 0

#combine
X_move <- as.data.frame(rbind(X0, X1, X2, X3))
colnames(X_move) <- c("x", "y")
plot(X_move)

#plot together and highlight points for different t values
set.seed(1)
ggplot()+
  geom_point(data = X_move[1:2,], aes(x, y), shape= 23, fill = "white",size = 5)+
  geom_point(data = X_move[1:2,], aes(x, y), col = c("red", "blue"))+
  geom_point(data = X_move[3:6,], aes(x, y), col = rep(c("red", "blue"), 2), alpha = .2, size = 2)+
  geom_jitter(data = X_move[7:8,], aes(x, y), col = c("red", "blue"),
              alpha = .2, size = 2, width = .0008, height = 0)+
  labs(x = "", y = "",
       title = )+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("./figures/latent_position_0.pdf",device = "pdf", width = 2.8, height = 2.8, units = "in")
ggsave("./figures/latent_position_0.jpeg",device = "jpeg", width = 2.8, height = 2.8, units = "in")

set.seed(1)
ggplot()+
  geom_point(data = X_move[3:4,], aes(x, y), shape= 23, fill = "white",size = 5)+
  geom_point(data = X_move[3:4,], aes(x, y), col = c("red", "blue"))+
  geom_point(data = X_move[c(1:2, 5:6),], aes(x, y), col = rep(c("red", "blue"), 2), alpha = .2, size = 2)+
  geom_jitter(data = X_move[7:8,], aes(x, y), col = c("red", "blue"),
              alpha = .2, size = 2, width = .0008, height = 0)+
  labs(x = "", y = "")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("./figures/latent_position_2.pdf",device = "pdf", width = 2.8, height = 2.8, units = "in")
ggsave("./figures/latent_position_2.jpeg",device = "jpeg", width = 2.8, height = 2.8, units = "in")

set.seed(1)
ggplot()+
  geom_point(data = X_move[5:6,], aes(x, y), shape= 23, fill = "white",size = 5)+
  geom_point(data = X_move[5:6,], aes(x, y), col = c("red", "blue"))+
  geom_point(data = X_move[1:4,], aes(x, y), col = rep(c("red", "blue"), 2), alpha = .2, size = 2)+
  geom_jitter(data = X_move[7:8,], aes(x, y), col = c("red", "blue"),
              alpha = .2, size = 2, width = .0008, height = 0)+
  labs(x = "", y = "")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("./figures/latent_position_6.pdf",device = "pdf",width = 2.8, height = 2.8, units = "in")
ggsave("./figures/latent_position_6.jpeg",device = "jpeg",width = 2.8, height = 2.8, units = "in")

set.seed(1)
ggplot()+
  geom_point(data = X_move[7:8,], aes(x, y), shape= 23, fill = "white",size = 5)+
  #geom_point(data = X_move[7:8,], aes(x, y), col = c("red", "blue"))+
  geom_point(data = X_move[1:6,], aes(x, y), col = rep(c("red", "blue"), 3), alpha = .2, size = 2)+
  geom_jitter(data = X_move[7:8,], aes(x, y), col = c("red", "blue"),
              alpha = 1, size = 2, width = .0008, height = 0)+
  labs(x = "", y = "")+
  theme_bw()+ 
  theme(axis.title.x=element_blank(),
        #axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("./figures/latent_position_1.pdf",device = "pdf",width = 2.8, height = 2.8, units = "in")
ggsave("./figures/latent_position_1.jpeg",device = "jpeg",width = 2.8, height = 2.8, units = "in")


