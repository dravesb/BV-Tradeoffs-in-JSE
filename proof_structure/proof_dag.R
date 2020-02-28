library('igraph')

# Proof DAG 
names <- c('Z form', 'Matrix\n Berstein','Spectral\n Bound', 'Eigenvalue\n Bound', 'Vtil', 'Hoeffding',
'Commuter', 'Residual\n 2, infty', 'Power\n Method', 'Residual\n Rows', 'Residual\n Exchangeability', 'Concentration', 'CLT')
lemma_names <- c('Lemma 1,2', 'Bernstein', 'Lemma 3', 'Lemma 4', 'Lemma 5', 'Lemma 6',
'Lemma 7', 'Lemma 8', 'Lemma 9', 'Lemma 10', 'Lemma 11', 'Theorem 1B', 'Theorem 2')

mat <- matrix(
c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0,
0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0,
0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
),
ncol = 13, nrow = 13, byrow = TRUE)

coords <- matrix(
c(
0.5,0,
1.5, 0,
0,2,
1,1.5,
2,2,
0,4,
1,3.5,
2,4,
0,6,
1,5.5,
2,6,
0.5,8,
1.5, 8
), 
nrow = 13, ncol  =2, byrow = TRUE)

G <- graph_from_adjacency_matrix(mat, mode = 'directed') %>% set_vertex_attr('label', value = names)
V(G)$label.cex <- 2; V(G)$size <- 40
E(G)$arrow.size <- 1 

jpeg('~/Documents/work/github/BJSE/proof_structure/proof_names.jpeg', height = 1080, width = 1080)
plot(G, vertex.color = adjustcolor('lightblue', alpha.f = .5), layout = coords)
dev.off()


G <- graph_from_adjacency_matrix(mat, mode = 'directed') %>% set_vertex_attr('label', value = lemma_names)
V(G)$label.cex <- 2; V(G)$size <- 40
E(G)$arrow.size <- 1  

jpeg('~/Documents/work/github/BJSE/proof_structure/lemma_names.jpeg', height = 1080, width = 1080)
plot(G, vertex.color = adjustcolor('pink', alpha.f = .5), layout = coords)
dev.off()





