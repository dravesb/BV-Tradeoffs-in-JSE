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
#     format data
#-----------------------------

df <- as.data.frame(read.csv("~/Documents/Work/github/BJSE/multiple_network_methods/data/data.csv")[,-1])

#group by iteration & method
df <- df %>% melt(id.vars = 1:3) %>% 
  group_by(network_size, t, variable) %>%  
  summarize(MC_Rate = mean(value), 
            se_rate = sd(value)/sqrt(n()))

colnames(df)[3] <- "Method"

#-----------------------------
#     make plot
#-----------------------------

ggplot(df, aes(t, MC_Rate, col = Method)) + 
  geom_point(size = 2, alpha = .3)+
  geom_line()+
  #geom_errorbar(aes(ymin = MC_Rate - se_rate, ymax =  MC_Rate + se_rate),alpha = .5,width = .05)+
  geom_ribbon(aes(ymin = MC_Rate - se_rate, ymax =  MC_Rate + se_rate),alpha = .1, linetype = 0)+
  theme_bw()+
  labs(x = "t", 
       y = "Misclassification Rate")


#-----------------------------
#     save plot
#-----------------------------

ggsave("multiple_methods_mc_comp.jpeg",
       width = 5, height = 4, units = "in", 
       device = "jpeg", 
       path = "~/Documents/Work/github/BJSE/multiple_network_methods/figures/")
