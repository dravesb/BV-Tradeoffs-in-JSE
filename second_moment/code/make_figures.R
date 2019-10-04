#--------------------------------------
#
#       Second Moment Simulations
#         Visualization
#           
#--------------------------------------

#-----------------------------
#     Preliminaries
#-----------------------------

#load up necessary packages
library(ggplot2)
library(dplyr)
library(reshape2)

#read in data
plotdf <- as.data.frame(read.csv("~/Documents/Work/github/BJSE/second_moment/data/plotting_data.csv")[,-1])

#-----------------------------
#     First Moment Residual
#-----------------------------

ggplot(plotdf, aes(R0X, R0Y, col = t))+
  geom_point(alpha = .1, size = .5)+
  facet_grid(rows = vars(Community),
             cols = vars(Graph))+
  theme_bw()

ggsave(filename = "first_moment_residual.pdf", 
       width = 6, height =6, 
       units = "in", 
       path = "~/Documents/Work/github/BJSE/second_moment/figures/")


#-----------------------------
#     Second Moment Residual
#-----------------------------

ggplot(plotdf, aes(R1X, R1Y, col = t))+
  geom_point(alpha = .1, size = .5)+
  facet_grid(rows = vars(Community),
             cols = vars(Graph))+
  theme_bw()

ggsave(filename = "root_n_Lhat_Ls.pdf", 
       width = 6, height =6, 
       units = "in", 
       path = "~/Documents/Work/github/BJSE/second_moment/figures/")


ggplot(plotdf, aes(R2X, R2Y, col = t))+
  geom_point(alpha = .2, size = .5)+
  facet_grid(rows = vars(Community),
             cols = vars(Graph))+
  theme_bw()

ggsave(filename = "root_n_Lhat_Z.pdf", 
       width = 6, height =6, 
       units = "in", 
       path = "~/Documents/Work/github/BJSE/second_moment/figures/")

ggplot(plotdf, aes(R3X, R3Y, col = t))+
  geom_point(alpha = .2, size = .5)+
  facet_grid(rows = vars(Community),
             cols = vars(Graph))+
  theme_bw()

ggsave(filename = "root_n_Z_Ls.pdf", 
       width = 6, height =6, 
       units = "in", 
       path = "~/Documents/Work/github/BJSE/second_moment/figures/")

#--------------------------------------
#
#       Power Method Visualization
#           
#--------------------------------------

#-----------------------------
#     Preliminaries
#-----------------------------

#load up necessary packages
library(ggplot2)
library(dplyr)
library(reshape2)

#read in data
plotdf <- as.data.frame(read.csv("~/Documents/Work/github/BJSE/second_moment/data/power_method_data.csv")[,-1])

#-----------------------------
#     Power Method Terms
#-----------------------------

ggplot(plotdf, aes(R1X, R1Y, col = t))+
  geom_point(alpha = .1, size = .5)+
  facet_grid(rows = vars(Community),
             cols = vars(Graph))+
  theme_bw()+
  labs("sqrt[n][(Lhat - L) - (A - P)M]")

ggsave(filename = "root_n_pm_M.pdf", 
       width = 6, height =6, 
       units = "in", 
       path = "~/Documents/Work/github/BJSE/second_moment/figures/")


ggplot(plotdf, aes(R2X, R2Y, col = t))+
  geom_point(alpha = .1, size = .5)+
  facet_grid(rows = vars(Community),
             cols = vars(Graph))+
  theme_bw()

ggsave(filename = "root_n_pm_Ls.pdf.pdf", 
       width = 6, height =6, 
       units = "in", 
       path = "~/Documents/Work/github/BJSE/second_moment/figures/")

#--------------------------------------
#
#       ZW - Ls CLT Visualization
#           
#--------------------------------------

#-----------------------------
#     Preliminaries
#-----------------------------

#load up necessary packages
library(ggplot2)
library(dplyr)
library(reshape2)

#read in data
plotdf <- as.data.frame(read.csv("~/Documents/Work/github/BJSE/second_moment/data/data_Z_minus_LS_CLT.csv")[,-1])

#plot residual term
ggplot(plotdf, aes(RX, RY, col = t))+
  geom_point(alpha = .1, size = .5)+
  facet_grid(rows = vars(Community),
             cols = vars(Graph))+
  theme_bw()

#plot clt term
ggplot(plotdf, 
       aes(sqrt(Network.Size) * RX, sqrt(Network.Size) * RY, 
           col = t))+
  geom_point(alpha = .1, size = .5)+
  facet_grid(rows = vars(Community),
             cols = vars(Graph))+
  labs(x = "x - Residual", y = "y - Residual")+
  theme_bw()

#plot residual term
ggplot(plotdf, aes(RX, RY, col = t))+
  geom_point(alpha = .1, size = .5)+
  facet_grid(rows = vars(Community),
             cols = vars(Network.Size))+
  theme_bw()

