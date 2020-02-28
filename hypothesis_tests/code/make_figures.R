#-----------------------------
#
#    Visualize Results
#
#-----------------------------

#-----------------------------
#    Behavior of Test Stats
#-----------------------------
library(dplyr); library(reshape2)

df <- read.csv("~/Documents/Work/github/BJSE/power_analysis/data/plotting_df.csv")[,-1]

#-----------------------------
#    Power curves
#-----------------------------
plotdf <- as.data.frame(df) %>%
  group_by(t, iter) %>% 
  summarize(power = mean(reject),
            std = sd(reject))

ggplot(plotdf, aes(t, power, group = iter, col = iter))+
  geom_point(alpha = 1, size = .1)+
  geom_line(alpha = .05)+
  theme_bw()+
  labs(x = "t", y = "Empirical Power", title = "Full Graph Hypothesis Testing Power Curve")
