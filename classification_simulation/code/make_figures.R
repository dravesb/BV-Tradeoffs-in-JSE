
#----------------------------------------
#
#       Plotting figures
#
#----------------------------------------
plotdf <- read.csv("~/Documents/Work/github/BJSE/classification_simulation/data/plotdf.csv", header = TRUE)[,-1]

library(ggplot2)
library(dplyr)

#------------------------------
#  Net Size increasing figures
#-----------------------------

ggplot(plotdf %>% filter(t == .5), aes(net_size, MC_Rate, col = Method))+
  geom_point(alpha = .5)+
  #geom_line(aes(linetype = Graph))+
  geom_line()+facet_grid(~Graph)+
  #geom_errorbar(aes(ymin=MC_Rate - mc_se,ymax=MC_Rate + mc_se), width = .1)+
  geom_ribbon(aes(ymin=MC_Rate - mc_se,ymax=MC_Rate + mc_se), alpha = .1, linetype = 0)+
  scale_x_log10()+
  labs(x = expression(paste('log'[10], "(Network Size)")),
       y = "Misclassification Rate")+
  theme_bw()


#comments: 
#
#       1. Network size increases, MC rate decreases
#       2. Omnibar and Abar looks like theyre clustered (probably can prove this)
#       3. Omni2 reflects changes in A2 but this lets increased performance of Omni1
#       4. Cases where Omni is better than Omnibar and Abar but Omnibar and Abar
#           seem to be the most reliable

ggsave(filename = "mc_rate_by_net_size.pdf", 
       width = 6, height =4, 
       units = "in", 
       path = "~/Documents/Work/github/BJSE/classification_simulation/figures/")

#------------------------------
#  Distance between centers
#     as a function of t
#-----------------------------
ggplot(plotdf %>% filter(net_size == 250), aes(t, Dist, col = Method))+
  geom_point(alpha = .5)+
  #geom_line(aes(linetype = Graph))+
  geom_line()+facet_grid(~Graph)+
  geom_vline(xintercept = .5, linetype = 2, col = "grey")+
  #geom_errorbar(aes(ymin= Dist - dist_se,ymax=Dist + dist_se), width = .03)+
  geom_ribbon(aes(ymin= Dist - dist_se,ymax=Dist + dist_se), alpha = .1, linetype = 0)+
  labs(x = "t", 
       y = "Distance between Centroids")+
  theme_bw()


#comments: 
#
#       1. t increase (more heterogeneous) distances decrease and clustering gets harder
#       2. Omni has the closests centroids in some cases
#       3. Omnibar looks to always have the furthest centers

ggsave(filename = "distance_by_t.pdf", 
       width = 6, height =4, 
       units = "in", 
       path = "~/Documents/Work/github/BJSE/classification_simulation/figures/")

#------------------------------
#  MC Rate by Distance
#-----------------------------
library(ggplot2)

ggplot(plotdf %>% filter(t == 0.5), aes(Dist, MC_Rate, col = Method))+
  geom_point()+
  geom_line()+
  facet_grid(~Graph)+
  scale_x_sqrt()+
  scale_y_sqrt()+
  #geom_ribbon(aes(ymin=MC_Rate - mc_se,
  #                ymax=MC_Rate + mc_se),
  #            alpha=0.1, linetype = 0)+
  labs(x = "Distance between Centroids", y = "Misclassification Rate")+
  theme_bw()

#comments: 
#
#       1. Distance increases, clustering performance improves
#       2. ASE > Omni > Omnibar = Abar have the smallest distances
#       3. But Omni looks to be performing better than other methods 
#           with similar distances - could speak to smaller varianes?

ggsave(filename = "mc_rate_by_distance.pdf", 
       width = 8, height =6, 
       units = "in", 
       path = "../figures")


  ggplot(plotdf %>% filter(t == 0.9), aes(net_size,Dist, col = Method))+
  geom_point()+
  geom_line()+
  facet_grid(~Graph)+
  scale_x_log10()+
  #geom_ribbon(aes(ymin=MC_Rate - mc_se,
  #                ymax=MC_Rate + mc_se),
  #            alpha=0.1, linetype = 0)+
 # labs(x = "Distance between Centroids", y = "Misclassification Rate")+
  theme_bw()

