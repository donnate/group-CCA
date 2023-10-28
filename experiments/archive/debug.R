
setwd("~/Documents/group-CCA")

library(tidyverse)
library(CCA)
library(VGAM)
library(matlib)
library(igraph)
library(recipes)
library(parallel)
library(tictoc)
library(PMA)
library(mvtnorm)
library(glmnet)
library(CCA)
library(pls)
library(igraph)
library(pracma)


seed = 1
name_exp = "debug.csv"
r = 2
N = 200
r_pca = 3

n= 200
p= 317
q=10 
sigma=0.1 
k=3
sigma_noise=0.1
power=1.6
lambdaA1seq= c(0.0001, 0.001, 0.01, 0.1, 1)
lambdaA2seq= c(0.0001, 0.001, 0.01, 0.1, 1)
conv=10^{-3} 
max.iter=200
egonet_size= 2
n.cv=5
power=1
rank=k
lambda1=0.5
lambda2=1
lambda3=0
effect_size =2
type_graph="pa"
probs = list('11'= 0.08, '12'=0.001, '13'=0.001, 
             '22'= 0.07, '23' = 0.001,
             '33' = 0.02)

vfn <- function(x){
  ifelse(x=="x", 1, -1)
}

