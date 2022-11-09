library(tidyverse)
library(CCA)
library(VGAM)

simple_experiment <- function(n, p, sigma, k, sigma_noise=0.1){
  ## This is a simple experiment to see how CCA (and variants) can perform in
  ## different settings.
  ## n      :   number of observations
  ## p      :   number of covariates
  ## sigma  :   noise
  ## k      :   number of canonical vectors
  
  Z = matrix(rnorm(n* k), n, k)
  coeffs_x <- runif(k, min=1, max=10)
  coeffs_y <- runif(k, min=1,max=10)
  X <- matrix(rnorm(n* p, sd=sigma_noise), n, p)
  X[, 1:k] = X[, 1:k] + Z  *  coeffs_x 
  Y <- matrix(rnorm(n* p, sd=sigma_noise), n, p)
  Y[, 1:k] = Y[, 1:k] + Z  *  coeffs_y 
  
  print(sapply(1:p, FUN=function(k){cor(X[,k], Y[,k])}))
  cc_results <- cancor(X,Y)
  return(list(X=X, Y=Y, Z=Z, cc_results=cc_results ))
}


simple_experiment <- function(n, p, sigma, k, sigma_noise=0.1, 
                              G = "", structure = "none"){
  ## This is a simple experiment to see how CCA (and variants) can perform in
  ## different settings.
  ## n      :   number of observations
  ## p      :   number of covariates
  ## sigma  :   noise
  ## k      :   number of canonical vectors
  
  Z = matrix(rnorm(n* k), n, k)  #### these serve has the building blocks.
  if (structure=="coeffs"){
    xcoef= runif(k, min=1, max=10) 
    ycoef= runif(k, min=1, max=10)
    X <- matrix(rnorm(n* p, sd=sigma_noise), n, p)
    X[, 1:k] = X[, 1:k] + Z  *  xcoef
    Y <- matrix(rnorm(n* p, sd=sigma_noise), n, p)
    Y[, 1:k] = Y[, 1:k] + Z  *  ycoef
  }
  
  coeffs = switch(  
    structure,  
    "None"= list(xcoef= drop_zeros(matrix(runif(k * p, min=1, max=10), p, k)), 
                 ycoef= runif(k, min=1, max=10)),  
    "Sparse"= list(xcoef= runif(k, min=1, max=10), 
                   ycoef= runif(k, min=1, max=10)), 
    "Group"= cat("Division = ", val1 / val2),  
    "Shrinked"= cat("Multiplication =", val1 * val2),
    "Cliques"= cat("Modulus =", val1 %% val2),
    "Graph"= cat("Power =", val1 ^ val2)
  )
  
  print(sapply(1:p, FUN=function(k){cor(X[,k], Y[,k])}))
  cc_results <- cancor(X,Y)
  return(list(X=X, Y=Y, Z=Z, cc_results=cc_results ))
}







