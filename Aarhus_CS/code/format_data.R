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

#load packages
pacman::p_load('igraph')

#set directory
setwd('~/Documents/Work/github/BJSE/Aarhus_CS/CS-Aarhus_Multiplex_Social/Dataset/')

#read in data
layer_ids <- read.table('CS-Aarhus_layers.txt', header = TRUE)
edge_list <- read.table('CS-Aarhus_multiplex.edges')
colnames(edge_list) <- c('layer_id', 'node1_id', 'node2_id', 'weight')
node_ids <- read.table('CS-Aarhus_nodes.txt', header = TRUE)

#function to go frome edge list to adjacency matrix
edge_list_to_adj_mat <- function(e_list, n){
  #get edges to place 
  in_edges <- as.matrix(e_list)
  out_edges <- as.matrix(e_list[2:1])
  
  #populate adjacency
  mat <- matrix(0, nrow = n, ncol  = n)
  mat[in_edges] <- mat[out_edges] <- 1
  
  return(mat)
}

#apply across layers
A_list <- list()
for(i in 1:nrow(layer_ids)){
  #fetch layer i's edge list
  layer_i_inds <- which(edge_list$layer_id == i)
  
  #get size of multiplex
  n <- nrow(node_ids)
  
  #fetch adjacency
  A_list[[i]] <- edge_list_to_adj_mat(edge_list[layer_i_inds, 2:3], n)
  
}

#write out each adjacency matrix
write.table(A_list[[1]], '../../adj_mats/A1.txt', row.names = FALSE, col.names = FALSE)
write.table(A_list[[2]], '../../adj_mats/A2.txt', row.names = FALSE, col.names = FALSE)
write.table(A_list[[3]], '../../adj_mats/A3.txt', row.names = FALSE, col.names = FALSE)
write.table(A_list[[4]], '../../adj_mats/A4.txt', row.names = FALSE, col.names = FALSE)
write.table(A_list[[5]], '../../adj_mats/A5.txt', row.names = FALSE, col.names = FALSE)






