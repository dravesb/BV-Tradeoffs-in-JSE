#--------------------------------------
#
#       Make Figures       
#         
#--------------------------------------
pacman::p_load(MCMCpack, ggplot2, dplyr, reshape2)

#------------------------------
#       Read in data 
#-----------------------------
plotdf <- read.csv("~/Documents/Work/github/BJSE/one_dim_BV_tradeoff/data/plotting_df.csv")[,-1]

#plot MSE figure
ggplot() +
  geom_point(aes(c, average_mse, col = method),
             plotdf %>% filter(method != 'Omnibar'),
             alpha = .1, size = .05)+
  geom_line(aes(c, theo_mse,
                col = method),
            plotdf%>% filter(method != 'Omnibar'),
            size = .5, alpha = .5)+
  facet_grid(cols = vars(graph))+
  theme_bw()+
  scale_y_log10()+
  labs(y = expression(paste('log'[10], "(MSE)")),
       x = "c",
       color = 'Embedding\n Method') 


ggsave(filename = "1d_mse.pdf", 
       width = 6, height = 4, 
       path = "~/Documents/Work/github/BJSE/one_dim_BV_tradeoff/figures/", 
       units = "in")

#plot bias figure
ggplot()+
  geom_line(aes(c, theo_bias, color = method, group = interaction(node_id, method)), 
            alpha = .1, size = .5,
            plotdf %>% filter(method != 'Omnibar'))+
  geom_point(aes(c, average_bias, color = method, group = interaction(node_id, method)), 
            alpha = .3, size = .1,
            plotdf %>% filter(method != 'Omnibar'))+
  facet_grid(cols = vars(graph))+
  theme_bw()+
  labs(x = 'c', y = 'Bias', 
       color = 'Embedding\n Method')

ggsave(filename = "1d_bias.pdf", 
       width = 6, height = 4, 
       path = "~/Documents/Work/github/BJSE/one_dim_BV_tradeoff/figures/", 
       units = "in")


#plot variance figure
ggplot()+
  geom_line(aes(c, theo_var, color = method, group = interaction(node_id, method)), 
            alpha = .1, size = .5,
            plotdf %>% filter(method != 'Omnibar'))+
            #plotdf)+
  geom_point(aes(c, empirical_var, color = method, group = interaction(node_id, method)), 
            alpha = .3, size = .1,
            plotdf %>% filter(method != 'Omnibar'))+
            #plotdf)+
  facet_grid(cols = vars(graph))+
  theme_bw()+
  scale_y_log10()+
  labs(y = expression(paste('log'[10], "(Variance)")),
       x = "c",
       color = 'Embedding\n Method') 

ggsave(filename = "1d_variance.pdf", 
       width = 6, height = 4, 
       path = "~/Documents/Work/github/BJSE/one_dim_BV_tradeoff/figures/", 
       units = "in")
