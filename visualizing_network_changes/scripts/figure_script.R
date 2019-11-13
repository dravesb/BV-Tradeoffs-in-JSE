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
  diag(c(t + 1, -t + 1))
}

#set up latent positions
n <- 50
com <- sample(c(1, 2), size = n, prob = c(pi, 1-pi), replace = TRUE)
X <- L[com,]

#set up P matrices
P1 <- tcrossprod(X, X %*% C(0)) # t = 0 
P2 <- tcrossprod(X, X %*% C(.2)) # t = 0.5
P3 <- tcrossprod(X, X %*% C(.6)) # t = .6
P3 <- tcrossprod(X, X %*% C(1)) # t = 1

#-----------------------------
#     Sample Adjmatrices
#-----------------------------

set.seed(1985)
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
g4 <- graph_from_adjacency_matrix(A3, mode = "undirected")

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
E(g1)$width <- 1
E(g2)$width <- 1
E(g3)$width <- 1
E(g4)$width <- 1

#change edge color
E(g1)$color <- "#55555555"
E(g2)$color <- "#55555555"
E(g3)$color <- "#55555555"
E(g4)$color <- "#55555555"

#plot networks
jpeg("../figures/g1.jpeg")
plot(g1, vertex.size = 5)
dev.off()

jpeg("../figures/g2.jpeg")
plot(g2, vertex.size = 5)
dev.off()

jpeg("../figures/g3.jpeg")
plot(g3, vertex.size = 5)
dev.off()

jpeg("../figures/g4.jpeg")
plot(g4, vertex.size = 5)
dev.off()














