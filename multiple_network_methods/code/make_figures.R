#---------------------------------------
#
#         Make Figures 
#
#---------------------------------------

#-----------------------------
#     read in libraries
#-----------------------------
library(ggplot2)
library(dplyr)
library(reshape2)

#-----------------------------
#     2 graph
#-----------------------------


#-----------------------------
#     format data
#-----------------------------

df <- as.data.frame(read.csv("~/Documents/Work/github/BJSE/multiple_network_methods/data/data.csv")[,-1])

#group by iteration & method
df <- as.data.frame(df) %>% melt(id.vars = 1:3) %>% 
  group_by(network_size, t, variable) %>%  
  summarize(MC_Rate = mean(value, na.rm = TRUE), 
            se_rate = 1.96*sd(value, na.rm = TRUE)/sqrt(n()))
            #se_rate = 1.96*sd(value, na.rm = TRUE))

colnames(df)[3] <- "Method"

#-----------------------------
#     make plots
#-----------------------------
a <- 0.0008
ggplot(df, aes(t, MC_Rate + a, col = Method)) + 
  geom_point(size = 2, alpha = .3)+
  geom_line()+
  #geom_errorbar(aes(ymin = MC_Rate - se_rate, ymax =  MC_Rate + se_rate),alpha = .5,width = .05)+
  geom_ribbon(aes(ymin = MC_Rate + a - se_rate, ymax =  MC_Rate + a + se_rate), alpha = .1, linetype = 0)+
  theme_bw()+
  scale_y_log10()+
  labs(x = "t", 
       y = expression(paste('log'[10], '(Misclassification Rate + ', epsilon, ')')))


#-----------------------------
#     save plot
#-----------------------------

ggsave("multiple_methods_mc_comp.pdf",
       width = 7, height = 5, units = "in", 
       device = "pdf", 
       path = "~/Documents/Work/github/BJSE/multiple_network_methods/figures/")


#-----------------------------
#     4 graph
#-----------------------------


#-----------------------------
#     format data
#-----------------------------

df <- as.data.frame(read.csv("~/Documents/Work/github/BJSE/multiple_network_methods/data/data_4_group.csv")[,-1])
df <- as.data.frame(read.csv("~/Documents/Work/github/BJSE/multiple_network_methods/data/data_3_group_w_rase.csv")[,-1])


#group by iteration & method
df <- df %>% melt(id.vars = 1:3) %>% 
  group_by(network_size, t, variable) %>%  
  summarize(MC_Rate = mean(value, na.rm = TRUE), 
            se_rate = 1.96*sd(value, na.rm = TRUE)/sqrt(n())) %>%
  mutate(Type = rep(c('Multiplex', 'Individual'), c(5, 4)))

colnames(df)[3] <- "Method"

#-----------------------------
#     make plots
#-----------------------------

ggplot(df %>% filter(Method != 'Omnicat'), 
       aes(t, MC_Rate, col = Method, linetype = Type)) + 
  geom_point(size = 2, alpha = .3)+
  geom_line()+
  #geom_errorbar(aes(ymin = MC_Rate - se_rate, ymax =  MC_Rate + se_rate),alpha = .5,width = .05)+
  geom_ribbon(aes(ymin = MC_Rate - se_rate, ymax =  MC_Rate + se_rate),alpha = .1, linetype = 0)+
  theme_bw()+
  scale_y_log10(limits = c(0.05, 0.67), breaks = seq(0.05, 0.7, by = 0.2))+
  labs(x = "t", 
       y = expression(paste('log'[10], '(Misclassification Rate)')))


#-----------------------------
#     save plot
#-----------------------------

ggsave("multiple_methods_mc_comp_new_B_4_graph.pdf",
       width = 7, height = 5, units = "in", 
       device = "pdf", 
       path = "~/Documents/Work/github/BJSE/multiple_network_methods/figures/")
       
       
       
#-----------------------------
#  Mahalanobis Distance Plots
#-----------------------------
#df <- read.csv("~/Documents/Work/github/BJSE/multiple_network_methods/data/data_3_group.csv")[,-1]
df <- as.data.frame(read.csv("~/Documents/Work/github/BJSE/multiple_network_methods/data/data_3_group_w_rase.csv")[,-1])
colnames(df) <- c('net_size', 't', 'iter', 'Method', 'Error', 'D12', 'D13', 'D23')

plotdf <- df %>%
	#filter(Method != 'Omnicat') %>% 
	melt(id.vars = c(1:5), value.name = 'Distance', variable.name = 'Group_Pair') %>% 
	mutate(Group_Pair = factor(Group_Pair, 
							   levels = c('D12', 'D13', 'D23'), 
							   labels = c('Mahalanobis Dist. Comm 1~2', 'Mahalanobis Dist. Comm 1~3', 'Mahalanobis Dist. Comm 2~3'))) %>%
	group_by(net_size, t, Method, Group_Pair) %>% 
	summarize(mean_error = mean(Error), 
		      error_upper = mean_error + 2 * sd(Error)/sqrt(n()), 
		      error_lower = mean_error - 2 * sd(Error)/sqrt(n()), 
		      mean_dist = mean(Distance), 
		      dist_upper = mean_dist + 2 * sd(Distance)/sqrt(n()), 
		      dist_lower = mean_dist - 2 * sd(Distance)/sqrt(n())
		      ) 


plotdf %>% 						   
	ggplot(aes(x = t, y = mean_error, col = Method, group = Method)) + 
	geom_ribbon(aes(ymax = error_upper, ymin = error_lower), 
             		fill = 'grey70', col = NA, alpha = 0.2) +
	geom_point(alpha = 0.8) + geom_line() +
	scale_y_log10(limits = c(0.05, 0.45), breaks = seq(0.05, 0.5, by = 0.1)) + 
	labs(y = expression(paste(log[10], '(Misclassification Rate)'))) + 
	theme_bw()

ggsave("multiple_methods_mc_3_group_w_rase.pdf",
       width = 7, height = 5, units = "in", 
       device = "pdf", 
       path = "~/Documents/Work/github/BJSE/multiple_network_methods/figures/")	
	
plotdf %>% 						   
	ggplot(aes(x = t, y = mean_dist, col = Method, group = Method)) + 
	geom_ribbon(aes(ymax = dist_upper, ymin = dist_lower), 
             		fill = 'grey70', col = NA, alpha = 0.2) +
	geom_point(alpha = 0.8) + geom_line() +
	facet_grid(.~Group_Pair) + 
	#scale_y_log10(limits = c(3.3, 7), breaks = seq(3, 7, by = 1)) + 
	labs(y = 'Mahalanobis Distance') + 
	theme_bw()

ggsave("multiple_methods_mahalanobis_3_group_w_rase.pdf",
       width = 7, height = 5, units = "in", 
       device = "pdf", 
       path = "~/Documents/Work/github/BJSE/multiple_network_methods/figures/")	



p1 <- plotdf %>% mutate(title = 'Misclassification Rate') %>% 
  filter(Method %in% c('ASE1','AASE', 'MASE', 'RASE', 'OMNI')) %>% 
	ggplot(aes(x = t, y = mean_error, col = Method, group = Method)) + 
	geom_ribbon(aes(ymax = error_upper, ymin = error_lower), 
             		fill = 'grey70', col = NA, alpha = 0.2) +
	geom_point(alpha = 0.5) + geom_line(alpha = 0.8, size = 0.5) +
	facet_grid(~title) + 
	scale_y_log10(limits = c(0.05, 0.45), breaks = seq(0.05, 0.5, by = 0.1)) + 
	labs(y = expression(paste(log[10], '(Misclassification Rate)'))) + 
	theme_bw() + 
	theme(legend.position = 'bottom')

	

p2 <- plotdf %>%  				   
  filter(Method %in% c('ASE1','AASE', 'MASE', 'RASE', 'OMNI')) %>% 
	ggplot(aes(x = t, y = mean_dist, col = Method, group = Method)) + 
	geom_ribbon(aes(ymax = dist_upper, ymin = dist_lower), 
             		fill = 'grey90', col = NA, alpha = 0.2) +
	geom_point(alpha = 0.5, size = 0.8 ) + geom_line(alpha = 0.8, size = 0.4) +
	facet_wrap(.~Group_Pair, ncol = 1) + 
	#scale_y_log10(limits = c(3.3, 7), breaks = seq(3, 7, by = 1)) + 
	labs(y = '') + 
	theme_bw() + 
	theme(legend.position = 'none')
	
lay <- cbind(rep(1, 3),rep(1, 3),
rep(2, 3))

pdf('~/Documents/Work/github/BJSE/multiple_network_methods/figures/multiple_methods_3_group_big_plot_w_rase.pdf', 
width = 8, height = 6, pointsize = 1
)
grid.arrange(grobs = list(p1, p2), layout_matrix = lay)
dev.off()	

p1 + labs(x = 'Network Heterogeneity') + theme(legend.position = 'right')
ggsave(filename = '~/Documents/Work/github/Thesis/Simulations/Chapter_3/Comm_Detection/figures/multiple_methods_3_group_big_plot_w_rase_slides.pdf',
       height = 4, width = 6, device = 'pdf', units = 'in')
	
	
	
	