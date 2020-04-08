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
pacman::p_load('igraph', 'RColorBrewer')

#set directory
setwd('~/Documents/Work/github/BJSE/Aarhus_CS/')

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

#create igraph objects
g1 <- graph_from_adjacency_matrix(A1, mode = 'undirected')
g2 <- graph_from_adjacency_matrix(A2, mode = 'undirected')
g3 <- graph_from_adjacency_matrix(A3, mode = 'undirected')
g4 <- graph_from_adjacency_matrix(A4, mode = 'undirected')
g5 <- graph_from_adjacency_matrix(A5, mode = 'undirected')

#-------------------------------------------------
#               Visualize Networks
#-------------------------------------------------
set.seed(1985)

jpeg(filename = paste0('./figures/graph_', layer_ids$layerLabel[1], '.jpeg'), width = 512, height = 512)
plot(g1, vertex.label = NA, vertex.size = 6, #layout = lay, 
     vertex.color = adjustcolor( colors[1], alpha.f = 0.5), main = layer_ids$layerLabel[1])
dev.off()

jpeg(filename = paste0('./figures/graph_', layer_ids$layerLabel[2], '.jpeg'), width = 512, height = 512)
plot(g2, vertex.label = NA, vertex.size = 6, #layout = lay, 
     vertex.color = adjustcolor(colors[2], alpha.f = 0.5),  main = layer_ids$layerLabel[2])
dev.off()

jpeg(filename = paste0('./figures/graph_', layer_ids$layerLabel[3], '.jpeg'), width = 512, height = 512)
plot(g3, vertex.label = NA, vertex.size = 6, #layout = lay, 
     vertex.color = adjustcolor( colors[3], alpha.f = 0.5), main = layer_ids$layerLabel[3])
dev.off()

jpeg(filename = paste0('./figures/graph_', layer_ids$layerLabel[4], '.jpeg'), width = 512, height = 512)
plot(g4, vertex.label = NA, vertex.size = 6, #layout = lay, 
     vertex.color = adjustcolor(colors[4], alpha.f = 0.5), main = layer_ids$layerLabel[4])
dev.off()

jpeg(filename = paste0('./figures/graph_', layer_ids$layerLabel[5], '.jpeg'), width = 512, height = 512)
plot(g5, vertex.label = NA, vertex.size = 6, #layout = lay, 
     vertex.color = adjustcolor( colors[5], alpha.f = 0.5),  main = layer_ids$layerLabel[5])
dev.off()

#full network 
edge_list<- read.table('./CS-Aarhus_Multiplex_Social/Dataset/CS-Aarhus_multiplex.edges')
graph <- graph_from_edgelist(as.matrix(edge_list[,2:3]), directed = FALSE)
curve_multiple(graph)
E(graph)$color <- sapply(colors[edge_list[,1]], function(x) adjustcolor(x, alpha.f = 0.5))

jpeg(filename = './figures/aarhus_multiplex.jpeg', width = 700, height = 0.75*700)
plot(graph, vertex.label = NA, vertex.size = 5, #layout = lay, 
     vertex.color = adjustcolor('white', alpha.f = 0.75), main = 'Aarhus Multiplex Network')
dev.off()

#-------------------------------------------------
#          Visualize Networks - Same Layout
#-------------------------------------------------
lay <- layout.fruchterman.reingold(g1)

jpeg(filename = paste0('./figures/graph_same_layout_', layer_ids$layerLabel[1], '.jpeg'), width = 512, height = 512)
plot(g1, vertex.label = NA, vertex.size = 6, layout = lay, 
     vertex.color = adjustcolor( colors[1], alpha.f = 0.5), main = layer_ids$layerLabel[1])
dev.off()

jpeg(filename = paste0('./figures/graph_same_layout_', layer_ids$layerLabel[2], '.jpeg'), width = 512, height = 512)
plot(g2, vertex.label = NA, vertex.size = 6, layout = lay, 
     vertex.color = adjustcolor(colors[2], alpha.f = 0.5),  main = layer_ids$layerLabel[2])
dev.off()

jpeg(filename = paste0('./figures/graph_same_layout_', layer_ids$layerLabel[3], '.jpeg'), width = 512, height = 512)
plot(g3, vertex.label = NA, vertex.size = 6, layout = lay, 
     vertex.color = adjustcolor( colors[3], alpha.f = 0.5), main = layer_ids$layerLabel[3])
dev.off()

jpeg(filename = paste0('./figures/graph_same_layout_', layer_ids$layerLabel[4], '.jpeg'), width = 512, height = 512)
plot(g4, vertex.label = NA, vertex.size = 6, layout = lay, 
     vertex.color = adjustcolor(colors[4], alpha.f = 0.5), main = layer_ids$layerLabel[4])
dev.off()

jpeg(filename = paste0('./figures/graph_same_layout_', layer_ids$layerLabel[5], '.jpeg'), width = 512, height = 512)
plot(g5, vertex.label = NA, vertex.size = 6, layout = lay, 
     vertex.color = adjustcolor( colors[5], alpha.f = 0.5),  main = layer_ids$layerLabel[5])
dev.off()

#full network 
edge_list<- read.table('./CS-Aarhus_Multiplex_Social/Dataset/CS-Aarhus_multiplex.edges')
graph <- graph_from_edgelist(as.matrix(edge_list[,2:3]), directed = FALSE)
curve_multiple(graph)
E(graph)$color <- sapply(colors[edge_list[,1]], function(x) adjustcolor(x, alpha.f = 0.5))

jpeg(filename = './figures/aarhus_multiplex_same_layout.jpeg', width = 700, height = 0.75*700)
plot(graph, vertex.label = NA, vertex.size = 5, layout = lay, 
     vertex.color = adjustcolor('white', alpha.f = 0.75), main = 'Aarhus Multiplex Network')
dev.off()

#-------------------------------------------------
#    Visualize Adjacency Matrices
#-------------------------------------------------

#plot heatmaps of adjacency matrices
jpeg(filename = paste0('./figures/adj_', layer_ids$layerLabel[1], '.jpeg'), width = 512, height = 512)
heatmap(A1, keep.dendro = FALSE, Rowv = NA, Colv = NA, symm = TRUE, 
        labRow = as.character(node_ids[,2]), labCol = as.character(node_ids[,2]),
        col = colorRampPalette(c('white', colors[1]))(100),
        main = layer_ids$layerLabel[1]
        )
dev.off()

jpeg(filename = paste0('./figures/adj_', layer_ids$layerLabel[2], '.jpeg'), width = 512, height = 512)
heatmap(A2, keep.dendro = FALSE, Rowv = NA, Colv = NA, symm = TRUE, 
        labRow = as.character(node_ids[,2]), labCol = as.character(node_ids[,2]),
        col = colorRampPalette(c('white', colors[2]))(100),
        main = layer_ids$layerLabel[2]
)
dev.off()

jpeg(filename = paste0('./figures/adj_', layer_ids$layerLabel[3], '.jpeg'), width = 512, height = 512)
heatmap(A3, keep.dendro = FALSE, Rowv = NA, Colv = NA, symm = TRUE, 
        labRow = as.character(node_ids[,2]), labCol = as.character(node_ids[,2]),
        col = colorRampPalette(c('white', colors[3]))(100),
        main = layer_ids$layerLabel[3]
)
dev.off()

jpeg(filename = paste0('./figures/adj_', layer_ids$layerLabel[4], '.jpeg'), width = 512, height = 512)
heatmap(A4, keep.dendro = FALSE, Rowv = NA, Colv = NA, symm = TRUE, 
        labRow = as.character(node_ids[,2]), labCol = as.character(node_ids[,2]),
        col = colorRampPalette(c('white', colors[4]))(100),
        main = layer_ids$layerLabel[4]
)
dev.off()

jpeg(filename = paste0('./figures/adj_', layer_ids$layerLabel[5], '.jpeg'), width = 512, height = 512)
heatmap(A5, keep.dendro = FALSE, Rowv = NA, Colv = NA, symm = TRUE, 
        labRow = as.character(node_ids[,2]), labCol = as.character(node_ids[,2]),
        col = colorRampPalette(c('white', colors[5]))(100),
        main = layer_ids$layerLabel[5]
)
dev.off()

#-------------------------------------------------
#    Degree Distributions
#-------------------------------------------------
make_dd <- function(graph, i = 1){
  jpeg(filename = paste0('./figures/dd_', layer_ids$layerLabel[i], '.jpeg'), width = 512, height = 512)
  hist(graph.strength(graph), col = colors[i], main = layer_ids$layerLabel[i], xlab = 'Vertex Degree')
  dev.off()
}
make_dd(g1, 1)
make_dd(g2, 2)
make_dd(g3, 3)
make_dd(g4, 4)
make_dd(g5, 5)

make_log_log <- function(graph, i = 1){
  dd <- degree.distribution(graph)
  d <- 0:max(degree(graph))
  ind <- which(dd > 0)[-1]
  jpeg(filename = paste0('./figures/log_log_', layer_ids$layerLabel[i], '.jpeg'), width = 512, height = 512)
  plot(d[ind], dd[ind], col = colors[i], main = paste('Degree Distribution:', layer_ids$layerLabel[i]),
       log = 'xy', xlab = 'Log-Degree', ylab = 'Log - Intensity', pch = 19)
  dev.off()
}
 
make_log_log(g1, 1)
make_log_log(g2, 2)
make_log_log(g3, 3)
make_log_log(g4, 4)
make_log_log(g5, 5)

#-------------------------------------------------
#    Summary Statistics
#-------------------------------------------------

names <- c('density', 'clustering','diameter', 'average_degree')
summary_stat_mat <- matrix(NA, nrow = 5, ncol = length(names))
colnames(summary_stat_mat) <- names 

#fill with summary statistics
g_list <- list(g1, g2, g3, g4, g4)
summary_stat_mat[,1] <- unlist(lapply(g_list, function(x) graph.density(x, loop = FALSE)))
summary_stat_mat[,2] <- unlist(lapply(g_list, function(x) transitivity(x, type = "average")))
summary_stat_mat[,3] <- unlist(lapply(g_list, function(x) diameter(x, directed = FALSE)))
summary_stat_mat[,4] <- unlist(lapply(g_list, function(x) mean(degree(x))))

summary_stat_mat



