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
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/je_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/mase_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/omni_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/mrdpg_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/omnicat_classification.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/model_setup_4_group.R")


#define storage
names <- c("network_size", "t","iter.no", 
           "Omnibar",'Omnicat', "MASE", "JE", "MRDPG",
           'ASE1', 'ASE2', 'ASE3', 'ASE4')
df <- matrix(NA, ncol = length(names), nrow = length(t) * mc_runs * length(net_size))
colnames(df) <- names
here <- 1 

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
      P1 <- tcrossprod(X)
      P2 <- X %*% tcrossprod(C1(t[i]), X)
      P3 <- X %*% tcrossprod(C2(t[i]), X)
      P4 <- X %*% tcrossprod(C3(t[i]), X)
      
      #sample adjaceny matrices
      A1 <- sampP(P1)
      A2 <- sampP(P2)
      A3 <- sampP(P3)
      A4 <- sampP(P4)
  
      #make adjaceny lists
      adj_mats <- list()
      adj_mats[[1]] <- A1
      adj_mats[[2]] <- A2
      adj_mats[[3]] <- A3
      adj_mats[[4]] <- A4
      
      #get classifications from each methods
      invisible(capture.output(omni_estimates <- omni_classes(adj_mats, 3, 3)))
      invisible(capture.output(omnicat_estimates <- omnicat_classes(adj_mats, 3, 3)))
      invisible(capture.output(mrdpg_estimates <- mrdpg_classes(adj_mats, 3, 3)))
      invisible(capture.output(je_estimates <- je_classes(adj_mats, 3, 3)))
      invisible(capture.output(mase_estimates <- mase_classes(adj_mats, 3, 3)))
      
      #get classifications from individual ASEs
      X1 <- ase(A1, 3); X2 <- ase(A2, 3)
      X3 <- ase(A3, 3); X4 <- ase(A4, 3)
      
      invisible(capture.output(ase1_estimates <- Mclust(X1, G = 3, modelNames = "VVV")$classification))
      invisible(capture.output(ase2_estimates <- Mclust(X2, G = 3, modelNames = "VVV")$classification))
      invisible(capture.output(ase3_estimates <- Mclust(X3, G = 3, modelNames = "VVV")$classification))
      invisible(capture.output(ase4_estimates <- Mclust(X4, G = 3, modelNames = "VVV")$classification))
      
      
      #get MC rates
      omni_mc <- get_mc3(omni_estimates, samp)
      omnicat_mc <- get_mc3(omnicat_estimates, samp)
      mrdpg_mc <- get_mc3(mrdpg_estimates, samp)
      je_mc <- get_mc3(je_estimates, samp)
      mase_mc <- get_mc3(mase_estimates, samp)    
      ase1_mc <- get_mc3(ase1_estimates, samp)    
      ase2_mc <- get_mc3(ase2_estimates, samp)
      ase3_mc <- get_mc3(ase3_estimates, samp)    
      ase4_mc <- get_mc3(ase4_estimates, samp)
      
      #store data
      df[here, ] <- c(net_size[l], t[i], j, omni_mc, omnicat_mc, mase_mc, je_mc, mrdpg_mc,
                      ase1_mc, ase2_mc, ase3_mc,ase4_mc)
      
      #update counter 
      here <- here + 1
      
      #print update
      if(here %% 100 == 0){
        print(round(here/nrow(df), 2))
      }
      
    }
  }
}

#-----------------------------
#     write out simulation
#-----------------------------

write.csv(df, "~/Documents/Work/github/BJSE/multiple_network_methods/data/data_4_group.csv")






