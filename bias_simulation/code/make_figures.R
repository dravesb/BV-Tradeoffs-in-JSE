plotdf <- read.csv("~/Documents/Work/github/BJSE/bias_simulation/data/plotting_data.csv", header = TRUE)[,-1]
#------------------------------------------------
#
#         Cirlce Funtions
#
#------------------------------------------------
library(dplyr)

circleFun <- function(center = c(0,0), r = 1, npoints = 100){
  tt <- seq(0,2*3.14,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

nsizes <- unique(plotdf %>% select(Network.Size))[,1]
f1 <- function(x) log(3*x)/sqrt(x)
f2 <- function(x) sqrt(log(3*x)/x)

circ_dat1 <- lapply(nsizes, function(x) circleFun(r = f1(x))) %>% bind_rows()
circ_dat1 <- circ_dat1 %>% mutate(Network.Size = rep(nsizes, each = 100)) 

circ_dat2 <- lapply(nsizes, function(x) circleFun(r = f2(x))) %>% bind_rows()
circ_dat2 <- circ_dat2 %>% mutate(Network.Size = rep(nsizes, each = 100)) 

#-------------------------------------------------
#
#        Bias Bounds
#
#-------------------------------------------------
library(ggplot2)

ggplot() +
  geom_point(aes(Xhat, Yhat, col = Community), plotdf, alpha = 0) +
  geom_point(aes(XS1, XS2),  shape = 8, size = 1, plotdf) +
  geom_point(aes(XC1, XC2),  shape = 3, size = 1, plotdf) +
  facet_grid(rows = vars(Graph),
             cols = vars(Network.Size))+ 
  scale_color_manual(values=c("red", "blue"))+
  theme_bw()+
  labs(x = "X", y = "Y",
       col = "Community")

ggsave("analytic_bias.jpeg", 
       width = 8, height = 6, units = "in", 
       device = "jpeg", 
       path = "~/Documents/Work/github/BJSE/bias_simulation/figures/")

ggplot() +
  geom_point(aes(Xhat, Yhat, col = Community), plotdf,
             size = .01, 
             alpha = .05) +
  geom_point(aes(XS1, XS2),  shape = 4, size = 1, plotdf) +
  geom_point(aes(XC1, XC2),  shape = 3, size = 1, plotdf) +
  facet_grid(rows = vars(Graph),
             cols = vars(Network.Size))+
  theme_bw()+
  theme(legend.position = "none")+
  scale_color_manual(values=c("red", "blue"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1)))+
  labs(x = "X", y = "Y",
       col = "Community")

ggsave("observed_bias.jpeg", 
       width = 8, height = 6, units = "in", 
       device = "jpeg", 
       path = "~/Documents/Work/github/BJSE/bias_simulation/figures/")

#-------------------------------------------------
#
#        Residual Bounds
#
#-------------------------------------------------

ggplot() +
  geom_point(aes(RX, RY, col = Community), plotdf, alpha = 0) +
  geom_path(aes(x, y), data = circ_dat1, linetype = "dashed")+
  facet_grid(rows = vars(Graph),
             cols = vars(Network.Size))+
  scale_color_manual(values=c("red", "blue"))+
  theme_bw()+
  labs(x = "X", y = "Y",
       col = "Community")

ggsave("residual_bounds.jpeg", 
       width = 8, height = 6, units = "in", 
       device = "jpeg", 
       path = "~/Documents/Work/github/BJSE/bias_simulation/figures/")

ggplot() +
  geom_point(aes(RX, RY, col = Community), plotdf,size = .01, alpha = .05) +
  geom_path(aes(x, y), data = circ_dat1, linetype = "dashed")+
  facet_grid(rows = vars(Graph),
             cols = vars(Network.Size))+
  scale_color_manual(values=c("red", "blue"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 1)))+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = "X", y = "Y",
       col = "Community")

ggsave("residuals.jpeg", 
       width = 8, height = 6, units = "in", 
       device = "jpeg", 
       path = "~/Documents/Work/github/BJSE/bias_simulation/figures/")


