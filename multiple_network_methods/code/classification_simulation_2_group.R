#--------------------------------------
#                                     
#       Multiple Graph Embedding      
#         Clustering Simulation       
#            2 Group                            
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
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/model_setup_ER_to_SBM.R")
source("~/Documents/Work/github/BJSE/multiple_network_methods/code/omnicat_classification.R")


#define storage
names <- c("network_size", "t","iter.no", "Omnibar", "MASE", "JE", "MRDPG", 'ASE1', 'ASE2')
df <- matrix(NA, ncol = length(names), nrow = length(t) * mc_runs * length(net_size))
colnames(df) <- names
here <- 1 

set.seed(1985)

#-----------------------------
#     Do simulation
#-----------------------------

for(l in 1:length(net_size)){#loop over net size
  for(i in 1:length(t)){#loop over t values
    for(j in 1:mc_runs){#monte carlo iterations
      
      #sample group assignments
      samp <- sample(1:2, net_size[l], replace = TRUE)
      X <- L[samp,]
      Xc <- X %*% sqrt(C(t[i]))
      
      #get different P matrices
      P1 <- tcrossprod(X)
      P2 <- tcrossprod(Xc)
      
      #sample adjaceny matrices
      A1 <- sampP(P1)
      A2 <- sampP(P2)
  
      #make adjaceny lists
      adj_mats <- list()
      adj_mats[[1]] <- A1
      adj_mats[[2]] <- A2
      
      #get classifications from multiplex methods
      invisible(capture.output(omni_estimates <- omni_classes(adj_mats, 2, 2)))
      #omnicat_estimates <- omnicat_classes(adj_mats, 2, 2)
      invisible(capture.output(mrdpg_estimates <- mrdpg_classes(adj_mats, 2, 2)))
      invisible(capture.output(je_estimates <- je_classes(adj_mats, 2, 2)))
      invisible(capture.output(mase_estimates <- mase_classes(adj_mats, 2, 2)))
      
      #get classifications from individual ASEs
      X1 <- ase(A1, 2); X2 <- ase(A2, 2)
      invisible(capture.output(ase1_estimates <- Mclust(X1, G = 2, modelNames = "VVV")$classification))
      invisible(capture.output(ase2_estimates <- Mclust(X2, G = 2, modelNames = "VVV")$classification))
            
      #get MC rates
      omni_mc <- get_mc(omni_estimates, samp)
      #omnicat_mc <- get_mc(omnicat_estimates, samp)
      mrdpg_mc <- get_mc(mrdpg_estimates, samp)
      je_mc <- get_mc(je_estimates, samp)
      mase_mc <- get_mc(mase_estimates, samp)    
      ase1_mc <- get_mc(ase1_estimates, samp)    
      ase2_mc <- get_mc(ase2_estimates, samp)    
      
      
      #store data
      #df[here, ] <- c(net_size[l], t[i], j, omni_mc, omnicat_mc, mase_mc, je_mc, mrdpg_mc, ase1_mc, ase2_mc)
      df[here, ] <- c(net_size[l], t[i], j, omni_mc, mase_mc, je_mc, mrdpg_mc, ase1_mc, ase2_mc)
      
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

write.csv(df, "~/Documents/Work/github/BJSE/multiple_network_methods/data/data.csv")






