#--------------------------------------
#
#       Power Analaysis Figures
#           
#--------------------------------------


#make figures
plotdf <- as.data.frame(df[,c(1:3,9:11)]) %>%
  melt(id.vars = c(1:3)) %>%
  group_by(t, net_size, variable) %>% 
  summarize(power = mean(value), sd_power = 1.96*sd(value)/sqrt(n()))

#plotdf$variable <- factor(plotdf$variable, levels = c('T_Rej', 'W_Rej', 'W_hat_rej'), labels = c('T', 'W', expression(hat(paste('W')))))

ggplot(plotdf, aes(t, power, col = variable))+
  geom_point(alpha = 0.25) + geom_line(alpha = 0.5)+
  geom_hline(yintercept = alpha, linetype = 'dashed', col = 'grey')+ 
  geom_ribbon(aes(ymin = power - sd_power,
                  ymax = power + sd_power),
              fill = 'grey', alpha = 0.25, linetype = 0)+
  facet_grid(rows = vars(net_size))+
  scale_color_manual(labels = c('T', 'W', expression(hat('W'))), values = c('red', 'green', 'blue'))+
  labs(x = 't', y = 'Empirical Power', color = 'Method')+
  theme_bw()



ggplot(as.data.frame(df[,2:6]) %>% melt(id.vars = 1:2) %>% filter(variable == 'W_hat'), 
       aes(value/net_size, col = variable))+
  geom_histogram(fill = 'white')+
  #geom_hline(yintercept = alpha, linetype = 'dashed', col = 'grey')+ 
  facet_grid(rows = vars(t),
             col = vars(net_size))+
  theme_bw()


