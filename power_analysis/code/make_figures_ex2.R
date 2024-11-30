#--------------------------------------
#
#       Power Analaysis Figures
#           
#--------------------------------------
pacman::p_load('lattice', 'ggplot2', 'gridExtra', 'grid', 'dplyr', 'reshape2')

#read in data
df1 <- read.csv('~/Documents/Work/github/BJSE/power_analysis/data/parameteric_test_results.csv')[,-1]
df2 <- read.csv('~/Documents/Work/github/BJSE/power_analysis/data/semi_parameteric_test_results.csv')[,-1]
df2 <- df2 %>% 
		mutate(Test_Statistic = 'T',
			   Data_Dependent = FALSE)

#set level
alpha <- 0.05

#oracle correction 
find_correction <- function(x, net_size){
  correction <- seq(0, 50)
  power <- numeric(length(correction))
  for(i in 1:length(correction)){
    power[i] <- sum(ifelse(x > qchisq(1 - alpha, df = (net_size * 2)  + correction[i]), 1, 0))/length(x)  
  }
  return(min(which(power < 0.05)))
}

#find oracle corrections
net_sizes <- c(50,100,200)
oracle_correction <- numeric(3)
for(i in 1:3){
  x <- df %>% filter(t == 0, net_size == net_sizes[i], Test_Statistic == 'W_hat') %>% .$Statistic
  oracle_correction[i] <- find_correction(x, net_sizes[i])
}

get_correction <- function(ns) ifelse(ns == 50, oracle_correction[1], ifelse(ns == 100, oracle_correction[2], oracle_correction[3]))

df_w_til <- df1 %>%
				filter(Test_Statistic == 'W_hat') %>% 
				rowwise() %>% 
				mutate(Threshold = qchisq(1 - alpha, df = (net_size * 2)  + get_correction(net_size)), 
				       Decision = ifelse(Statistic > Threshold, 1, 0)) %>% 
				ungroup() %>% 
				mutate(Test_Statistic = 'W_til', Data_Dependent = FALSE)
	
#power figure
plotdf <- rbind(df1, df2, df_w_til) %>% 	
	group_by(net_size, t, Test_Statistic, Data_Dependent) %>% 
	summarize(Empirical_Power = mean(Decision),
            lower = mean(Decision) - 2 * sd(Decision)/sqrt(n()), 
            upper = mean(Decision) + 2 * sd(Decision)/sqrt(n())
            #lower = Empirical_Power - 2 *sqrt(Empirical_Power*(1 - Empirical_Power)/sqrt(n())), 
            #upper = mean(Decision) + 2 *sqrt(Empirical_Power*(1 - Empirical_Power)/sqrt(n()))
            ) 
            
plotdf %>%
  ggplot(aes(x = t, y = Empirical_Power, col = Test_Statistic, group = Test_Statistic)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              fill = 'grey90', alpha = .5, colour = NA) +
  geom_line(aes(linetype = Data_Dependent)) + geom_point(alpha = 0.8) + 
  geom_hline(yintercept = 0.05, linetype = 'dashed', col = 'grey') +
  scale_color_manual(labels = c('W', expression(hat('W')), 'T', expression(tilde('W'))), 
                     values = c('red','blue','purple', 'darkgreen'))+
  facet_grid(~net_size, scales = 'free', labeller = label_bquote(rows = 'n ='~.(net_size))) +
  labs(x = 't', y = 'Empirical Power', col = 'Statistic', linetype = 'Data\nDependent') + 
  theme_bw()


ggsave(filename = 'ex2_empirical_power_oracle_correction.pdf', 
       width = 10, height = 5, 
       path = "~/Documents/Work/github/BJSE/power_analysis/figures/", 
       units = "in")

#--------------------------------------
#
#       Variable t Figures
#           
#--------------------------------------

#read in data
df <- read.csv('~/Documents/Work/github/BJSE/power_analysis/data/plotting_df_variable_t.csv')[,-1]
#df <- df %>% mutate(W_hat_Rej_new = ifelse(W_hat > qchisq(1 - alpha, df = net_size * 2  + correction), 1, 0))

#set level
alpha <- 0.05

#power figure
plotdf <- as.data.frame(df[,c(1:3,9:11)]) %>%
  melt(id.vars = c(1:3)) %>%
  group_by(t, net_size, variable) %>% 
  summarize(power = mean(value), sd_power = 1.96*sd(value)/sqrt(n()))

ggplot(plotdf, aes(t, power, col = variable))+
  geom_ribbon(aes(ymin = power - sd_power,
                  ymax = power + sd_power),
              fill = 'grey', alpha = 0.8, linetype = 0)+
  geom_point(alpha = .5) + geom_line(alpha = 0.5)+
  geom_hline(yintercept = alpha, linetype = 'dashed', col = 'grey')+ 
  scale_color_manual(labels = c('T', 'W', expression(hat('W')))
                     , values = c('red', 'darkgreen', 'blue')
  )+
  facet_grid(cols = vars(net_size), scales = "free_x")+
  labs(x = 't', y = 'Empirical Power', color = 'Method')+
  theme_bw()

ggsave(filename = 'empirical_power_variable_t.pdf', 
       width = 10, height = 5, 
       path = "~/Documents/Work/github/BJSE/power_analysis/figures/", 
       units = "in")

#oracle correction visualization
find_correction <- function(x, net_size){
  correction <- seq(0, 50)
  power <- numeric(length(correction))
  for(i in 1:length(correction)){
    power[i] <- sum(ifelse(x > qchisq(1 - alpha, df = net_size * 2  + correction[i]), 1, 0))/length(x)  
  }
  return(min(which(power < 0.05)))
}

#find oracle corrections
net_sizes <- c(50,100,200)
oracle_correction <- numeric(3)
for(i in 1:3){
  x <- (df %>% filter(t == 0, net_size == net_sizes[i]) %>% select(6))[,1]
  oracle_correction[i] <- find_correction(x, net_sizes[i])
}

get_correction <- function(ns) ifelse(ns == 50, oracle_correction[1], ifelse(ns == 100, oracle_correction[2], oracle_correction[3]))

df2 <- df %>% mutate(W_hat_Rej_new = ifelse(W_hat > qchisq(1 - alpha, df = net_size * 2  + get_correction(net_size)), 1, 0))

#power figure
plotdf <- as.data.frame(df2[,c(1:3,9:12)]) %>%
  melt(id.vars = c(1:3)) %>%
  group_by(t, net_size, variable) %>% 
  summarize(power = mean(value), sd_power = 1.96*sd(value)/sqrt(n()))

ggplot(plotdf, aes(t, power, col = variable))+
  geom_ribbon(aes(ymin = power - sd_power,
                  ymax = power + sd_power),
              fill = 'darkgrey', alpha = 0.5, linetype = 0)+
  geom_point(alpha = .5) + geom_line(alpha = 0.5)+
  geom_hline(yintercept = alpha, linetype = 'dashed', col = 'grey')+ 
  scale_color_manual(labels = c('T', 'W', expression(hat('W')), expression(tilde('W'))), 
                     values = c('red','blue','purple', 'darkgreen'))+
  facet_grid(cols = vars(net_size), scales = "free_x")+
  labs(x = 't', y = 'Empirical Power', color = 'Method')+
  theme_bw()

ggsave(filename = 'empirical_power_oracle_correction_variable_t.pdf', 
       width = 10, height = 5, 
       path = "~/Documents/Work/github/BJSE/power_analysis/figures/", 
       units = "in")

#relative power figure
plotdf <- as.data.frame(df2[,c(1:3,9:12)]) %>%
  group_by(t, net_size) %>% 
  summarize(T_power = mean(T_Rej), 
            T_sd_power = 1.96*sd(T_Rej)/sqrt(n()),
            W_power = mean(W_Rej), 
            W_sd_power = 1.96*sd(W_Rej)/sqrt(n()),
            W_hat_power = mean(W_hat_Rej), 
            W_hat_sd_power = 1.96*sd(W_hat_Rej)/sqrt(n()),
            W_til_power = mean(W_hat_Rej_new), 
            W_til_sd_power = 1.96*sd(W_hat_Rej_new)/sqrt(n())
            ) %>%
  melt(id.vars = c(1:4))

#reshape
power_df <- plotdf[plotdf$variable %in% c('W_power', 'W_hat_power', 'W_til_power'),]
sd_df <- plotdf[!(plotdf$variable %in% c('W_power', 'W_hat_power', 'W_til_power')),]
plotdf <- power_df %>% mutate(sd = sd_df$value)
#plotdf$variable <- droplevels(plotdf$variable)

ggplot(plotdf, aes(T_power, value, col = variable))+
   geom_ribbon(aes(ymin = value - sd, ymax = value + sd),
               fill = 'darkgrey', alpha = 0.5, linetype = 0)+
  geom_point(alpha = 0.5)+ geom_line(alpha = 0.5)+
  geom_hline(yintercept = alpha, linetype = 'dashed', col = 'black', alpha = 0.5)+
  geom_abline(slope = 1, linetype = 'dashed', col = 'black', alpha = 0.5)+
  coord_equal()+
  scale_color_manual(labels = c('W', expression(hat('W')), expression(tilde('W'))), 
                     values = c('red','darkgreen','blue'))+
  facet_grid(cols = vars(net_size))+
  labs(x = 'T Empirical Power', y = 'Empirical Power', color = 'Method')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(filename = 'empirical_power_oracle_correction_variable_t.pdf', 
       width = 10, height = 4, 
       path = "~/Documents/Work/github/BJSE/power_analysis/figures/", 
       units = "in")

#power(W)/power(T)
plotdf <- as.data.frame(df2[,c(1:3,9:12)]) %>%
  group_by(t, net_size) %>% 
  summarize(W_power = (mean(W_Rej) - mean(T_Rej))/mean(T_Rej),
            W_hat_power = (mean(W_hat_Rej) - mean(T_Rej))/mean(T_Rej), 
            W_til_power = (mean(W_hat_Rej_new) - mean(T_Rej))/mean(T_Rej)
  ) %>%
  melt(id.vars = c(1:2))

ggplot(plotdf %>% filter(variable == 'W_hat_power'),
       aes(t, value, col = variable))+
  geom_point(alpha = 0.5) + geom_line(alpha = 0.5)+
  #geom_hline(yintercept = 1, linetype = 'dashed', col = 'black', alpha = 0.5)+
  scale_color_manual(labels = c('W', expression(hat('W')), expression(tilde('W'))), 
                     values = c('red','darkgreen','blue'))+
  facet_grid(cols = vars(net_size),scales = "free_x")+
  labs(x = 't', y = 'power(W)/power(T)', color = 'Method')+
  theme_bw()
