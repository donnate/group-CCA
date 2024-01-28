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
source("elena/reduced_rank_regression.R")
source("elena/graph_reduced_rank_regression.R")
#store all values above diagonal of connectome matrices in matrix c

args <- commandArgs(trailingOnly=TRUE)
lambda <- as.numeric(args[1])
test_fold <- as.numeric(args[2])
val_fold <- ifelse(test_fold < 16, test_fold + 1, 1)

X <- as.matrix(vroom("data/activations_X_preprocessed.csv"))
behaviour <- read_csv("data/activations_Y_preprocessed.csv")
group_assignment <- readxl::read_xlsx("data/activation_groups.xlsx", col_names = FALSE)
folds = t(read_csv("data/folds.csv"))

colnames(group_assignment) <- "group.id"
index_groups = which(group_assignment$group.id!=0)
groups <- sapply(1:length(unique(group_assignment$group.id[index_groups])),
                 function(g){which(group_assignment$group.id[index_groups] == g)})

Y = as.matrix(behaviour)
n = nrow(X)
p = ncol(X)
q = ncol(Y)
Sy = NULL
LW_Sy = FALSE
rho = 1
niter = 100
r = 2
thresh = 1e-3
##### Split into different 
index_test = as.numeric(folds[test_fold,])
index_val = as.numeric(folds[val_fold,])
index_train = (1:n)[-c(index_test, index_val)]


if (is.null(Sy)){
  ###
  Sy = t(Y[index_train,]) %*% Y[index_train,] /n
  if (LW_Sy){
    lw_cov <- corpcor::cov.shrink(Y)
    Sy <- as.matrix(lw_cov)
  }
}



svd_Sy = svd(Sy)
sqrt_inv_Sy = svd_Sy$u %*% diag(sapply(svd_Sy$d, function(x){ifelse(x > 1e-4, 
                                                                    1/sqrt(x), 0)}))  %*% t(svd_Sy$v)
tilde_Y = Y %*% sqrt_inv_Sy




print("Using CVXR")
print("lambda is")
print(lambda)

# large_groups = c(19, 18, 15, 13, 14,11, 12)
# for (g in large_groups){
#  print(g)
#  Sg = t(X[index_train, groups[[g]]]) %*%X[index_train, ] /n
#  write_csv(data.frame(Sg), paste0("~/Downloads/groups_Sg", g, ".csv" ))
# }


#STOP
svd_X = svd(X[index_train, ]/sqrt(n))
svd_left = (svd_X$v) %*% diag(1/ (svd_X$d + rho))
{

  U = matrix(0, p, q)
  Z = matrix(0, p, q)
  B = matrix(0, p, q)
  
  prod_xy = t(X[index_train,]) %*% tilde_Y[index_train,]/length(index_train)
  for (i in 1:niter){
    Uold = U
    Zold = Z
    Bold = B
    time1<-tic()
    B = svd_left  %*% t(svd_X$v) %*% (prod_xy  + rho * (Z - U))
    time2<-toc()
    print(paste0("iteration ", i, " time = ", time2))
    print(max(B))
    Z = B + U
    norm_col = sapply(1:length(groups), function(g){sqrt(sum(Z[groups[[g]],]^2))})
    for (g in 1:length(groups)){
      if(norm_col[g] < lambda * sqrt(length(groups[[g]]))){
        Z[groups[[g]],] = 0
        print(g)
      }else{
        #print(c(g,  (1- (lambda * sqrt(length(groups[[g]])) /rho)/norm_col[g]) ))
        Z[groups[[g]],] =  (1- (lambda * sqrt(length(groups[[g]])) /rho)/norm_col[g]) * Z[groups[[g]],]
      }
      print(paste("Max of group ", g, " is ", max( Z[groups[[g]],])))
    }
    U = U + B - Z
    print(c("ADMM iter", i, norm(Z - B)/sqrt(p), norm(Zold - Z)/sqrt(p), norm(Uold - U)/sqrt(p)))
    if (max(c(norm(Z - B), norm(Zold - Z)/sqrt(p))) <thresh){
      break
    }
  }
  B_opt = B
}

### groups 3, 4, 7,, 8, 10 , 11, 15, 18,19 should be computed ahead
B_opt[which(abs(B_opt)<1e-5)] = 0
print(B_opt)

sol = svd((svd_X$v) %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, sqrt(x), 0)}))  %*% t(svd_X$v)  %*% B_opt)
#sol = svd(t(B_OLS[I, ]), nu = r, nv=r)
V = sqrt_inv_Sy %*% sol$v[, 1:r]
inv_D = diag(sapply(1:r, FUN=function(x){ifelse(sol$d[x]<1e-4, 0, 1/sol$d[x])}))
U = B_opt %*% sol$v[, 1:r] %*% inv_D ### = U\lambda

correlation <-
  data.frame("CCA_graph_rrr",
             lambda,
              test_fold,
             "val_fold" = val_fold,
             diag(cov(as.matrix(X)[train_index,] %*% U,
                      as.matrix(Y)[train_index,] %*%  V)),
             apply(((as.matrix(X)[train_index,] %*% U) -
                      (as.matrix(Y)[train_index,] %*%  V))^2, 2, mean),
             diag(t(as.matrix(X)[test_index,] %*% U) %*%
                    as.matrix(Y)[test_index,] %*% V),
             apply(((as.matrix(X)[test_index,] %*% U) -
                      (as.matrix(Y)[test_index,] %*%  V))^2, 2, mean),
             diag(cor(as.matrix(X)[test_index,] %*% U),
                    as.matrix(Y)[test_index,] %*% V),
             diag(t(as.matrix(X)[val_index,] %*% U) %*%
                    as.matrix(Y)[val_index,] %*% V),
             apply(((as.matrix(X)[val_index,] %*% U) -
                      (as.matrix(Y)[val_index,] %*%  V))^2, 2, mean),
             diag(cor(as.matrix(X)[val_index,] %*% U),
                  as.matrix(Y)[val_index,] %*% V)
  )

write_csv(correlation, paste0("data/results_l", lambda, "_test_fold", test_fold, ".csv"))


