#--------------------------------------
#                                     
#       Multiple Graph Embedding      
#         Clustering Simulation       
#           4 Group                          
#--------------------------------------

#-----------------------------
#     Preliminaries
#-----------------------------

#load up necessary packages
library(ggplot2)
library(dplyr)
library(reshape2)
library(matlib)
library(irlba)
library(mclust)

#source nesseary functions
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/ase_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/je_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/mase_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/abar_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/omni_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/mrdpg_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/rase_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/model_setup_4_group.R")


#define storage
df_list <- list()
count <- 1

names <- c("network_size", "t","iter.no",
           'Method', 'MC_rate', 'Distance'
           )

#set seed
set.seed(1985)

#-----------------------------
#     Do simulation
#-----------------------------

for(l in 1:length(net_size)){#loop over n
  for(i in 1:length(t)){#loop over t 
    for(j in 1:mc_runs){#monte carlo iterations
      
      #sample group assignments
      samp <- sample(1:3, net_size[l], replace = TRUE)
      X <- L[samp,]
      
      #get different P matrices
      P_list <- list(tcrossprod(X), 
                     X %*% tcrossprod(C1(t[i]), X), 
                     X %*% tcrossprod(C2(t[i]), X), 
                     X %*% tcrossprod(C3(t[i]), X)
      )
      
      
      #sample adjaceny matrices
      adj_mats <- lapply(P_list, sampP)
      
      #get classifications from each methods
      invisible(capture.output(ase1_estimates <- ase_classes(adj_mats[[1]], 3, 3)))
      invisible(capture.output(abar_estimates <- abar_classes(adj_mats, 3, 3)))
      invisible(capture.output(omni_estimates <- omni_classes(adj_mats, 3, 3)))
      #invisible(capture.output(omnicat_estimates <- omnicat_classes(adj_mats, 3, 3)))
      invisible(capture.output(rase_estimates <- rase_classes(adj_mats, 3, 3)))
      invisible(capture.output(mrdpg_estimates <- mrdpg_classes(adj_mats, 3, 3)))
      invisible(capture.output(je_estimates <- je_classes(adj_mats, 3, 3)))
      invisible(capture.output(mase_estimates <- mase_classes(adj_mats, 3, 3)))
      
      #get MC rates
      abar_mc <- get_mc3(abar_estimates$clusters, samp)
      omni_mc <- get_mc3(omni_estimates$clusters, samp)
      #omnicat_mc <- get_mc3(omnicat_estimates$clusters, samp)
      mrdpg_mc <- get_mc3(mrdpg_estimates$clusters, samp)
      je_mc <- get_mc3(je_estimates$clusters, samp)
      mase_mc <- get_mc3(mase_estimates$clusters, samp)    
      ase1_mc <- get_mc3(ase1_estimates$clusters, samp)    
      rase_mc <- get_mc3(rase_estimates$clusters, samp)    
      #ase2_mc <- get_mc3(ase2_estimates, samp)
      #ase3_mc <- get_mc3(ase3_estimates, samp)    
      #ase4_mc <- get_mc3(ase4_estimates, samp)
      
      
      
      #store data
      df_list[[count]] <- rbind(c(net_size[l], t[i], j, 'AASE', abar_mc, abar_estimates$distances),
                                c(net_size[l], t[i], j, 'OMNI', omni_mc, omni_estimates$distances),
                                #c(net_size[l], t[i], j, 'Omnicat', omnicat_mc, omnicat_estimates$distances),
                                c(net_size[l], t[i], j, 'ALS', mrdpg_mc, mrdpg_estimates$distances),
                                c(net_size[l], t[i], j, 'JE', je_mc, je_estimates$distances),
                                c(net_size[l], t[i], j, 'RASE', rase_mc, rase_estimates$distances),
                                c(net_size[l], t[i], j, 'MASE', mase_mc, mase_estimates$distances),
                                c(net_size[l], t[i], j, 'ASE1', ase1_mc, ase1_estimates$distances)
                                )
                      
      
      #update counter 
      count <- count + 1
      
      }
    #print update
    print(round(i/length(t), 2))
    
  }
}

#-----------------------------
#     write out simulation
#-----------------------------
df <- Reduce('rbind', df_list)
colnames(df) <- c('network_size','t', 'iter.no', 'Method', 'MC_Rate', 'D12', 'D13', 'D23')

write.csv(df, "~/Documents/Work/github/BJSE/multiple_network_methods/data/data_3_group_w_rase.csv")
#write.csv(df, "~/Documents/Work/github/BJSE/multiple_network_methods/data/data_3_group.csv")

as.data.frame(df) %>% 
  group_by(t, Method) %>% 
  summarize(MC_Rate_mean = mean(as.numeric(as.character(MC_Rate)))) %>% 
  ggplot(aes(t, MC_Rate_mean, col = Method, group = Method)) + 
  geom_point() + geom_line() + 
  theme_bw()




