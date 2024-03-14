library(tidyverse)
library(RNifti)
library(oro.nifti)
library(neurobase)
library(tidyverse)
library(igraph)
library(ggplot2)
library(vroom)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
setwd("~/Documents/group-CCA/")

source('experiments/sparse_CCA/experiment_functions.R')
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source("elena/missing/helper.R")
source("elena/missing/evaluation.R")
source("elena/missing/original_CCA_impute.R")
source("elena/gradient_descent.r")
source("elena/iterative_cca.R")
source("elena/graph_utils.R")
source("elena/reduced_rank_regression.R")
source("elena/graph_reduced_rank_regression.R")
source("elena/group_reduced_rank_regression.R")
#store all values above diagonal of connectome matrices in matrix c
library(vroom)
args <- commandArgs(trailingOnly=TRUE)
part <- as.numeric(args[1])

load("~/Downloads/folds_neuro.RData")
positions = read_csv("~/Downloads/parcellation_coordinates.csv")
real_labels <- readxl::read_xlsx("~/Downloads/BNA_subregions.xlsx") %>% 
  filter(`Label ID.L`<211)
real_labels["name"] = sapply(real_labels$...6, function(x){str_split(x, ",")[[1]][1]})
real_labels["full_label"] = sapply(real_labels$...6, function(x){str_split(x, ",")[[1]][2]})

positions["name"] = sapply(positions$Region, function(x){return(str_split(x, "_")[[1]][2])})
positions  = positions %>% left_join(real_labels %>% dplyr::select(`...6`,
                                                                   name,
                                                                   `Gyrus...10`),
                                     by = c("name"))

positions = positions %>%
  dplyr::rename( "Gyrus" = `Gyrus...10`)
positions = positions %>%
  dplyr::select(-c("...6"))
positions$Gyrus[1:19] = positions$name[1:19] 
positions$Gyrus[15] = "Brain Stem" 
positions$name[15] = "Brain Stem" 

distanceMatrix <- as.matrix(dist(positions[, c("x", "y", "z")]))

# Find the 6 nearest neighbors for each point
numNeighbors <- 4
edges <- c()
for (i in 1:nrow(distanceMatrix)){
  x = distanceMatrix[i,] 
  edges = rbind(edges,
                t(rbind(rep(i, numNeighbors), order(x)[1:numNeighbors+1])))
}
edges = data.frame(edges)
colnames(edges) <- c("Source", "Target")
# Create a graph from these edges
g <- graph_from_edgelist(as.matrix(edges), directed = FALSE)

# Optionally, plot the graph
#plot(g)
# Find the connected components
comp <- components(g)
# Get the number of connected components
num_components <- comp$no
print(num_components)
Gamma <- get_edge_incidence(g, weight = 1)
Gamma_dagger = pinv(Gamma)


edges2 <- c()
for (i in 1:nrow(positions)){
  neighbours= which(positions$Gyrus == positions$Gyrus[i])
  edges2 = rbind(edges2,
                 t(rbind(rep(i, length(neighbours)), neighbours)))
}
edges2 = data.frame(edges2)
colnames(edges2) <- c("Source", "Target")
edges2 = edges2 %>%
  filter(Source < Target)
# Create a graph from these edges
g2 <- graph_from_edgelist(as.matrix(edges2), directed = FALSE)

# Optionally, plot the graph
#plot(g2)
# Find the connected components
comp <- components(g2)
# Get the number of connected components
num_components <- comp$no
print(num_components)
Gamma2 <- get_edge_incidence(g2, weight = 1)
Gamma_dagger2 = pinv(Gamma2)

unique_regions = unique(positions$Gyrus)
groups <- sapply(1:length(unique_regions),
                 function(i){which(positions$Gyrus == unique_regions[i])}
)
##### Split into different 
print("here yeah")
n = nrow(X)
p = ncol(X)
q = ncol(Y)
do.scale = T
if(q >n){
  X_temp <- X
  X <- Y
  Y <- X_temp
}
if (do.scale){
  X <- scale(X)
  Y <- scale(Y)
}
#write_csv(data.frame(X), "data/activations_X_preprocessed.csv")
#write_csv(data.frame(Y), "data/activations_Y_preprocessed.csv")

#write_csv(t(data.frame(folds)), "data/folds.csv")
r = 3
##### Split into different 
### Let's do a cross validation setting
correlation<-c()
print('Starting')
for (i in  part:min(length(folds), (part + 10))){
  print(i)
  test_fold <- i
  val_fold <- ifelse(test_fold < length(length(folds)), test_fold + 1, 1)
  index_test = as.numeric(folds[[test_fold]])
  index_val = as.numeric(folds[[val_fold]])
  index_train = (1:n)[-c(index_test, index_val)]
  for (lambda in 10^seq(-3,2, length.out = 20)){
    final = CCA_rrr(as.matrix(X)[index_train,], 
                    as.matrix(Y)[index_train,],  
                    lambda = lambda, 
                    Kx=NULL, r=r,
                    rho=1, niter=2 * 1e4,
                    do.scale = FALSE, lambda_Kx=0,
                    thresh=1e-5,
                    solver= "ADMM",
                    LW_Sy = TRUE)
    
    correlation <- rbind(
      correlation,
      c("RRR",
        lambda,
        test_fold,
        val_fold,
        diag(cov(as.matrix(X)[index_train,] %*% final$U,
                 as.matrix(Y)[index_train,] %*%  final$V)),
        apply(((as.matrix(X)[index_train,] %*% final$U) -
                 (as.matrix(Y)[index_train,] %*%  final$V))^2, 2, mean),
        diag(t(as.matrix(X)[index_test,] %*% final$U) %*%
               as.matrix(Y)[index_test,] %*% final$V),
        apply(((as.matrix(X)[index_test,] %*% final$U) -
                 (as.matrix(Y)[index_test,] %*%  final$V))^2, 2, mean),
        diag(cor(as.matrix(X)[index_test,] %*% final$U,
                 as.matrix(Y)[index_test,] %*% final$V)),
        diag(t(as.matrix(X)[index_val,] %*% final$U) %*%
               as.matrix(Y)[index_val,] %*% final$V),
        apply(((as.matrix(X)[index_val,] %*% final$U) -
                 (as.matrix(Y)[index_val,] %*%  final$V))^2, 2, mean),
        diag(cor(as.matrix(X)[index_val,] %*% final$U,
                 as.matrix(Y)[index_val,] %*% final$V))
      )
    )
    
    final = CCA_graph_rrr(as.matrix(X)[index_train,], 
                          as.matrix(Y)[index_train,],  
                          Gamma =Gamma,
                          lambda = lambda, 
                          Kx=NULL, r=r,
                          rho=1, niter=1 * 1e4,
                          do.scale = FALSE, lambda_Kx=0,
                          thresh=1e-5,
                          LW_Sy = TRUE)
    
    correlation <- rbind(
      correlation,
      c("RRR-graph1",
        lambda,
        test_fold,
        val_fold,
        diag(cov(as.matrix(X)[index_train,] %*% final$U,
                 as.matrix(Y)[index_train,] %*%  final$V)),
        apply(((as.matrix(X)[index_train,] %*% final$U) -
                 (as.matrix(Y)[index_train,] %*%  final$V))^2, 2, mean),
        diag(t(as.matrix(X)[index_test,] %*% final$U) %*%
               as.matrix(Y)[index_test,] %*% final$V),
        apply(((as.matrix(X)[index_test,] %*% final$U) -
                 (as.matrix(Y)[index_test,] %*%  final$V))^2, 2, mean),
        diag(cor(as.matrix(X)[index_test,] %*% final$U,
                 as.matrix(Y)[index_test,] %*% final$V)),
        diag(t(as.matrix(X)[index_val,] %*% final$U) %*%
               as.matrix(Y)[index_val,] %*% final$V),
        apply(((as.matrix(X)[index_val,] %*% final$U) -
                 (as.matrix(Y)[index_val,] %*%  final$V))^2, 2, mean),
        diag(cor(as.matrix(X)[index_val,] %*% final$U,
                 as.matrix(Y)[index_val,] %*% final$V))
      )
    )
    
    final = CCA_graph_rrr(as.matrix(X)[index_train,], 
                          as.matrix(Y)[index_train,],  
                          Gamma =Gamma2,
                          lambda = lambda, 
                          Kx=NULL, r=r,
                          rho=1, niter=1 * 1e4,
                          do.scale = FALSE, lambda_Kx=0,
                          thresh=1e-5,
                          LW_Sy = TRUE)
    
    correlation <- rbind(
      correlation,
      c("RRR-graph2",
        lambda,
        test_fold,
        val_fold,
        diag(cov(as.matrix(X)[index_train,] %*% final$U,
                 as.matrix(Y)[index_train,] %*%  final$V)),
        apply(((as.matrix(X)[index_train,] %*% final$U) -
                 (as.matrix(Y)[index_train,] %*%  final$V))^2, 2, mean),
        diag(t(as.matrix(X)[index_test,] %*% final$U) %*%
               as.matrix(Y)[index_test,] %*% final$V),
        apply(((as.matrix(X)[index_test,] %*% final$U) -
                 (as.matrix(Y)[index_test,] %*%  final$V))^2, 2, mean),
        diag(cor(as.matrix(X)[index_test,] %*% final$U,
                 as.matrix(Y)[index_test,] %*% final$V)),
        diag(t(as.matrix(X)[index_val,] %*% final$U) %*%
               as.matrix(Y)[index_val,] %*% final$V),
        apply(((as.matrix(X)[index_val,] %*% final$U) -
                 (as.matrix(Y)[index_val,] %*%  final$V))^2, 2, mean),
        diag(cor(as.matrix(X)[index_val,] %*% final$U,
                 as.matrix(Y)[index_val,] %*% final$V))
      )
    )
    
    final = CCA_group_rrr(as.matrix(X)[index_train,], 
                          as.matrix(Y)[index_train,],  
                          groups =groups,
                          lambda = lambda, 
                          Kx=NULL, r=r,
                          rho=1, niter=1 * 1e4,
                          do.scale = FALSE, lambda_Kx=0,
                          thresh=1e-5,
                          LW_Sy = TRUE)
    
    correlation <- rbind(
      correlation,
      c("RRR-group",
        lambda,
        test_fold,
        val_fold,
        diag(cov(as.matrix(X)[index_train,] %*% final$U,
                 as.matrix(Y)[index_train,] %*%  final$V)),
        apply(((as.matrix(X)[index_train,] %*% final$U) -
                 (as.matrix(Y)[index_train,] %*%  final$V))^2, 2, mean),
        diag(t(as.matrix(X)[index_test,] %*% final$U) %*%
               as.matrix(Y)[index_test,] %*% final$V),
        apply(((as.matrix(X)[index_test,] %*% final$U) -
                 (as.matrix(Y)[index_test,] %*%  final$V))^2, 2, mean),
        diag(cor(as.matrix(X)[index_test,] %*% final$U,
                 as.matrix(Y)[index_test,] %*% final$V)),
        diag(t(as.matrix(X)[index_val,] %*% final$U) %*%
               as.matrix(Y)[index_val,] %*% final$V),
        apply(((as.matrix(X)[index_val,] %*% final$U) -
                 (as.matrix(Y)[index_val,] %*%  final$V))^2, 2, mean),
        diag(cor(as.matrix(X)[index_val,] %*% final$U,
                 as.matrix(Y)[index_val,] %*% final$V))
      )
    )
    
  }
  
  for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                   "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                   "SCCA_Parkhomenko")){
    
    print(paste0("Starting ", method))
    tryCatch({
      test1<-additional_checks(as.matrix(X)[c(index_train,
                                              index_test),], 
                               as.matrix(Y)[c(index_train,
                                              index_test),],
                               S=NULL, 
                               rank=r, kfolds=length(folds)-1, 
                               method.type = method,
                               lambdax= 10^seq(-3,2, length.out = 20),
                               lambday = c(0))
      correlation <- rbind(
        correlation,
        c(method,
          0,
          test_fold,
          val_fold,
          diag(cov(as.matrix(X)[index_train,] %*% test1$u,
                   as.matrix(Y)[index_train,] %*%  test1$v)),
          apply(((as.matrix(X)[index_train,] %*% test1$u) -
                   (as.matrix(Y)[index_train,] %*%  test1$v))^2, 2, mean),
          diag(t(as.matrix(X)[index_test,] %*% test1$u) %*%
                 as.matrix(Y)[index_test,] %*% test1$v),
          apply(((as.matrix(X)[index_test,] %*% test1$u) -
                   (as.matrix(Y)[index_test,] %*%  test1$v))^2, 2, mean),
          diag(cor(as.matrix(X)[index_test,] %*% test1$u,
                   as.matrix(Y)[index_test,] %*% test1$v)),
          diag(t(as.matrix(X)[index_val,] %*% test1$u) %*%
                 as.matrix(Y)[index_val,] %*% test1$v),
          apply(((as.matrix(X)[index_val,] %*% test1$u) -
                   (as.matrix(Y)[index_val,] %*%  test1$v))^2, 2, mean),
          diag(cor(as.matrix(X)[index_val,] %*% test1$u,
                   as.matrix(Y)[index_val,] %*% test1$v))
        )
      )
    }, error = function(e) {
      # Print the error message
      cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
      # Skip to the next iteration
    })
    correlation_df = data.frame(correlation)
    write_csv(correlation_df, file=paste0("~/Downloads/corr_3d_", part, ".csv"))
  }
}






