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
df <- df %>% melt(id.vars = 1:3) %>% 
  group_by(network_size, t, variable) %>%  
  summarize(MC_Rate = mean(value, na.rm = TRUE), 
            se_rate = 1.96*sd(value, na.rm = TRUE)/sqrt(n()))
            #se_rate = 1.96*sd(value, na.rm = TRUE))

colnames(df)[3] <- "Method"

#-----------------------------
#     make plots
#-----------------------------

ggplot(df, aes(t, MC_Rate, col = Method)) + 
  geom_point(size = 2, alpha = .3)+
  geom_line()+
  #geom_errorbar(aes(ymin = MC_Rate - se_rate, ymax =  MC_Rate + se_rate),alpha = .5,width = .05)+
  geom_ribbon(aes(ymin = MC_Rate - se_rate, ymax =  MC_Rate + se_rate), alpha = .1, linetype = 0)+
  theme_bw()+
  scale_y_log10()+
  labs(x = "t", 
       y = expression(paste('log'[10], '(Misclassification Rate)')))


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
