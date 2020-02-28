#-----------------------------------
#
#   2D MSE Visualizations
#
#-----------------------------------

#set up basics
pacman::p_load(ggplot2, dplyr)
setwd('/Users/benjamindraves/Documents/Work/github/BJSE/two_dim_BV_tradeoff/data/')

#read in data
plotdf <- read.csv("plotting_df.csv")[,-1]

#plot MSE figure
ggplot() +
  geom_point(aes(t, average_mse, col = method),
             plotdf,
             alpha = .1, size = .05)+
  geom_line(aes(t, theo_mse,
                col = method),
            plotdf,
            size = .5, alpha = .5)+
  facet_grid(cols = vars(graph),
             rows = vars(comm))+
  theme_bw()+
  scale_y_log10()+
  labs(y = expression(paste('log'[10], "(MSE)")),
       x = "t",
       color = 'Embedding\n Method') 

ggsave(filename = "2d_mse.pdf", 
       width = 6, height = 4, 
       path = "~/Documents/Work/github/BJSE/two_dim_BV_tradeoff/figures/", 
       units = "in")

#plot bias figure
ggplot()+
  geom_line(aes(t, theo_bias2, color = method, group = interaction(node_id, method)), 
            alpha = .1, size = .5,
            plotdf)+
  geom_point(aes(t, average_bias2, color = method, group = interaction(node_id, method)), 
             alpha = .3, size = .1,
             plotdf)+
  facet_grid(cols = vars(graph),
             rows = vars(comm))+
  theme_bw()+
  scale_y_sqrt()+
  labs(x = 't', y = expression(paste('Bias'^2)), 
       color = 'Embedding\n Method')

ggsave(filename = "2d_bias.pdf", 
       width = 8, height = 6, 
       path = "~/Documents/Work/github/BJSE/two_dim_BV_tradeoff/figures/", 
       units = "in")


#plot variance figure
ggplot()+
  geom_line(aes(t, theo_var, color = method, group = interaction(node_id, method)), 
            alpha = .1, size = .5,
            plotdf)+
            #plotdf %>% filter(method %in% c('Omni', 'Omnibar')))+
  geom_point(aes(t, empirical_var, color = method, group = interaction(node_id, method)), 
             alpha = .3, size = .1,
             plotdf)+
             #plotdf %>% filter(method %in% c('Omni', 'Omnibar')))+
  facet_grid(cols = vars(graph),
             rows = vars(comm))+
  theme_bw()+
  scale_y_log10()+
  labs(y = expression(paste('log'[10], "(Variance)")),
       x = "t",
       color = 'Embedding\n Method') 

ggsave(filename = "2d_variance.pdf", 
       width = 8, height = 6, 
       path = "~/Documents/Work/github/BJSE/two_dim_BV_tradeoff/figures/", 
       units = "in")


#-------------------------------------------------

#sample rows of X
Lc <- function(j) L %*% sqrt(C(t[j]))

mat <- matrix(NA, ncol = 5, nrow = 4*length(t))   
colnames(mat) <-  c('X', 'Y', 't', 'comm', 'graph')
for(j in 1:length(t)){
  mat[(4*j-3):(4*j - 2),] <- cbind(Lc(j), c(t[j], t[j]), 1:2, rep('Graph 2', 2))
  mat[(4*j-1):(4*j),] <- cbind(L, c(t[j], t[j]), 1:2, rep('Graph 1', 2))
}

ggplot() +
  geom_point(aes(average_estimate_X, average_estimate_Y, col = t),
             plotdf %>% filter(method == 'ASE', graph == 'Graph 2'),
             alpha = .1)+
  geom_line(aes(X, Y, color = t, group = graph), 
            as.data.frame(mat))+
  #geom_line(aes(theo_bias_X, theo_bias_Y, col = t),
  #           plotdf %>% filter(method == 'Omni'),
  #           alpha = 1, size = .5)+
  theme_bw()+  
  labs(color = 'ER = 1\nSBM = 0 ',
       title = 'ASE Graph 2: Theoretic Lines & Observed Consentration',
       x = 'X', y = 'Y'
       ) 


#plot bias figure
ggplot()+
  geom_point(aes(bias_X, bias_Y, color = method, group = interaction(node_id, method)), 
            alpha = .1, size = .5,
            plotdf)+
  facet_grid(cols = vars(graph),
             rows = vars(comm))+
  theme_bw()+
  #scale_y_sqrt()+
  labs(x = 'Bias_X', y = expression(paste('Bias_Y')), 
       color = 'Embedding\n Method')


mat_df <- as.data.frame(mat)
mat_df[,1] <- as.numeric(as.character(mat_df[,1]))
mat_df[,2] <- as.numeric(as.character(mat_df[,2]))
mat_df[,3] <- as.numeric(as.character(mat_df[,3]))


plotdf2 <- merge(plotdf, as.data.frame(mat), by = c('graph', 'comm', 't'), all.x = TRUE)

#plot bias figure
ggplot(plotdf2)+
  geom_point(aes(average_estimate_X, average_estimate_Y, color = method, group = interaction(node_id, method)), 
             alpha = .1, size = .5)+
  geom_line(aes(as.numeric(as.character(X)), as.numeric(as.character(Y)), group = comm), alpha = .5) +
  geom_point(aes(as.numeric(as.character(X)), as.numeric(as.character(Y))), alpha = .5) +
  facet_grid(cols = vars(graph))+
  theme_bw()+
  #scale_y_sqrt()+
  labs(x = 'X', y = expression(paste('Y')), 
       color = 'Embedding\n Method')




