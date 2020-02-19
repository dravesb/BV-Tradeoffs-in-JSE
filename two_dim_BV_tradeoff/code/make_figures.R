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

ggplot() +
  geom_point(aes(average_estimate_X, average_estimate_Y, col = t),
             plotdf %>% filter(method == 'Omni'),
             alpha = .1, size = .05)+
  #geom_line(aes(theo_bias_X, theo_bias_Y, col = t),
  #           plotdf %>% filter(method == 'Omni'),
  #           alpha = 1, size = .5)+
  facet_grid(#rows = vars(graph),
             cols = vars())+
  theme_bw()+  
  labs(color = 't') 







