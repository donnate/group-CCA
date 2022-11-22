library(tidyverse)
library(CCA)
library(VGAM)
library(matlib)

### Under development

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


simple_experiment <- function(n, p, q, r, rhos=c(), sigma, sigma_noise=0.1, 
                              G = NA, structure = "none", prob=0.1){
  ## This is a simple experiment to see how CCA (and variants) can perform in
  ## different settings.
  ## n      :   number of observations
  ## p      :   number of covariates
  ## sigma  :   noise
  ## k      :   number of canonical vectors

  if (r> min(p,q)){
    r = min(p,q)
  }
  
  if (p<q){
    p0 = p
    p = q
    q <- p0
  }
  if (length(rhos)==0){
    rhos = runif(r, min=0.5, max=1)
  }elif(length(rhos)<r){
    rhos <- c(rhos, runif(r-length(rhos), min=0.5, max=1))
  }else{
    rhos = rhos[1:r]
  }
  
  Z = matrix(rnorm(n* k), n, k)  #### these serve has the building blocks.
  if (structure=="coeffs"){
    xcoef= matrix(runif(r * p, min=1, max=10), p, r)
    #### We need to ensure that the coefficients are orthogonal.
    ycoef= matrix(runif(r * q, min=1, max=10), q, r)
    X0 <- matrix(rnorm(n* p, sd=1), n, p)
    X0 <- X0 - apply(X0, 2, mean)
    #### Make sure that ||Xa||=1
    tilde_xcoef = xcoef/apply(X0 %*% xcoef,2, FUN=function(x){sqrt(sum(x^2))})
    Y0 =  X0 %*% tilde_xcoef %*% diag(rhos) %*%inv(t(ycoef) %*% (ycoef)) 
    #Y0 = Y0/ norm(Y0 %*% ycoef, "F")
    Y= Y0 + matrix(rnorm(n* q, sd=sigma_noise), n, q)
    X= X0 + matrix(rnorm(n* p, sd=sigma_noise), n, p)
    cc_results <- cancor(X0,Y0,xcenter = FALSE, ycenter = FALSE)
    
  }
  if (structure == "sparse"){
    coeffs = sparsify(matrix(runif(k * p, min=1, max=10), p, k), prob)
    
  }
  
  
  print(sapply(1:p, FUN=function(k){cor(X[,k], Y[,k])}))
  cc_results <- cancor(X,Y)
  return(list(X=X, Y=Y, Z=Z, cc_results=cc_results ))
}



sparsify <- function(x, prob){
  mask = matrix(rbinom(dim(x)[1] * dim(x)[2], 1, prob), dim(x)[1], dim(x)[2])
  return(x * mask)
}



