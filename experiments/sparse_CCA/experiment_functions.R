library(MASS)
library(stats)
library(CVXR)
library(geigen)
library(pracma)
library(tidyverse)
library(CCA)
library(VGAM)
library(matlib)
library(PMA)
library(mvtnorm)
library(glmnet)
library(caret)
#Example of TGD on Sparse CCA
#n = 500, p1 = p2 = 100, s_u = s_v = 5
#k = 20, eta = 0.0025, lambda =0.01, T = 12000
wd = getwd()
setwd("experiments/GCA")
source("utils.R")
source("gca_to_cca.R")
source("init_process.R")
source("sgca_init.R")
source("sgca_tgd.R")
source("subdistance.R")
source("adaptive_lasso.R")

setwd(wd)
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source("src/cca_bis.R")

setwd(wd)
generate_example <- function(n, p1, p2,   nnzeros = 5,
                             theta = diag( c(0.9,  0.8)),
                             a = 0, r=2, signal_strength="normal"){
  #n  <- 200;
  #p1 <- 100;
  #p2 <- 100;
  # r <- 2
  p <- p1 + p2;
  pp <- c(p1,p2);
  s  <- sample(1:min(p1,p2),nnzeros);
  print('--------------------------------------');
  print('Generating data ...');
  #a <- 0.3;
  Sigma <- diag(p1+p2)
  
  # generate covariance matrix for X and Y
  T1 = toeplitz(a^(0:(pp[1]-1)));
  T1[which(T1<1e-6)] = 0
  Sigma[1:p1, 1:p1] = T1;
  #T1 = Sigma[1:p1, 1:p1]
  Tss = T1[s,s];
  u = matrix(0, pp[1], r)
  u[s,1:r] <- as.matrix(runif( nnzeros * r,max = 5, min=3), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  u <- u %*%(sqrtm(t(u[s,1:r]) %*% Tss %*% u[s,1:r])$Binv)
  
  T2 = toeplitz(a^(0:(pp[2]-1)));
  Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)] = T2;
  Tss = Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)][s, s];
  v = matrix(0, pp[2], r)
  v[s,1:r] <- as.matrix(runif( nnzeros * r,max = 5, min=3), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  v <- v %*%(sqrtm(t(v[s,1:r]) %*% Tss %*% v[s,1:r])$Binv)
  
  Sigma[(p1+1):(p1+p2), 1:p1] = T2 %*%  v  %*% theta %*% t(u) %*% T1;
  Sigma[1:p1, (p1+1):(p1+p2)] = t(Sigma[(p1+1):(p1+p2), 1:p1])
  
  Sigmax = Sigma[1:p1,1:p1];
  Sigmay = Sigma[(p1+1):p,(p1+1):p];
  
  #Generate Multivariate Normal Data According to Sigma
  Data = mvrnorm(n, rep(0, p), Sigma);
  
  X = Data[,1:p1];
  Y = Data[,(p1+1):(p1+p2)];
  
  print('Data generated.');
  print('--------------------------------------');
  
  Mask = matrix(0, p, p);
  idx1 = 1:pp[1];
  idx2 = (pp[1]+1):(pp[1]+pp[2]);
  Mask[idx1,idx1] <- matrix(1,pp[1],pp[1]);
  Mask[idx2,idx2] <- matrix(1,pp[2],pp[2]);
  Sigma0 = Sigma * Mask;
  
  
  S <- cov(Data)
  sigma0hat <- S * Mask
  
  # Estimate the subspace spanned by the largest eigenvector using convex relaxation and TGD
  # First calculate ground truth
  result = geigen::geigen(Sigma,Sigma0)
  evalues <- result$ values
  evectors <-result$vectors
  evectors <- evectors[,p:1]
  a <- evectors[,1:r]
  #scale <- a %*% sqrtm(diag(r)+t(a) %*% Sigma %*% a/lambda)$B;
  return(list(Sigma=Sigma, Sigma0=Sigma0,
         S = S, sigma0hat =  sigma0hat, Mask= Mask,
         X=X, Y = Y, Data=Data,u=u, v=v, 
         Sigmax=Sigmax, Sigmay=Sigmay,
         a=a
        ))
}

cv_function <- function(X, Y, kfolds=5, initu, initv,
                        lambdax, adaptive=TRUE) {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- numeric(length = kfolds)
  # loop over folds
  for (i in seq_along(folds)) {
    # split data into training and validation sets
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]
    
    # fit model on training data with hyperparameters
    tryCatch(
    {
    model <- adaptive_lasso(X_val, Y_val %*% initv, initu, adaptive=adaptive, 
                         lambdax, 
                         max.iter=5000, 
                         max_k = 10, verbose = FALSE, ZERO_THRESHOLD=1e-5)
    
    # make predictions on validation data
    # compute RMSE on validation data
    rmse[i] <- sum((X_val %*% model$Uhat - Y_val%*% initv)^2)
    },
    error = function(e) {
      # If an error occurs, assign NA to the result
      rmse[i] <- NA
    })
  }
  
  # return mean RMSE across folds
  return(median(rmse))
}


cv_function_tgd <- function(X, Y, Mask, kfolds=5, ainit,
                        lambda, r=2, k=20,  maxiter=1000, eta=0.001,
                        convergence=1e-3) {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- numeric(length = kfolds)
  p1 <- dim(X)[2]
  p2 <- dim(Y)[2]
  p <- p1 + p2;
  n <- nrow(X)
  pp <- c(p1,p2);
  S0 = cov(cbind(X, Y)
           )
  
  #init <- gca_to_cca(ainit, S0, pp)
  # loop over folds
  for (i in seq_along(folds)) {
    # split data into training and validation sets
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]
    
    S = cov(cbind(X_train, Y_train))
    sigma0hat = S * Mask
    
    # fit model on training data with hyperparameters
    tryCatch(
    {
    final = sgca_tgd(A=S, B=sigma0hat,
             r=r,ainit, k, lambda = lambda, eta=eta,
             convergence=convergence,
             maxiter=maxiter, plot = FALSE, 
             scale=NULL)
    final <- gca_to_cca(final, S, pp)
    
    # make predictions on validation data
    # compute RMSE on validation data
    rmse[i] <- sum((X_val %*% final$u - Y_val%*% final$v)^2)
    },
    error = function(e) {
      # If an error occurs, assign NA to the result
      rmse[i] <- sum((X_val %*% final$u - Y_val%*% final$v)^2)
    })
  }
  
  # return mean RMSE across folds
  return(median(rmse))
}

pipeline_adaptive_lasso <- function(Data, Mask, sigma0hat, r, nu=1, Sigmax, 
                                    Sigmay, maxiter=30, lambdax=NULL, lambday=NULL,
                                    adaptive=TRUE, kfolds=5, param1=10^(seq(-4, 2, by = 0.25)),
                                    create_folds=TRUE){
  
  ### data splitting procedure 3 folds
  Data=example$Data;Mask = example$Mask; sigma0hat=example$sigma0hat; r=2; 
  nu=1; Sigmax = example$Sigmax;  Sigmay = example$Sigmay; maxiter=maxiter=100;  lambdax=NULL;
  adaptive=TRUE; kfolds=5;  param1=10^(seq(-5, 1, by = 0.5));
  create_folds=FALSE
  
  p1 <- dim(Sigmax)[1]
  p2 <- dim(Sigmay)[1]
  p <- p1 + p2;
  n <- nrow(Data)
  pp <- c(p1,p2);
  if(create_folds){
    folds <- createFolds(1:nrow(Data), k = 2, list = TRUE, returnTrain = FALSE)
    S1 <- cov(Data[folds[[1]],])
    S3 <- cov(Data[folds[[1]],])
    X = Data[folds[[2]],1:p1]
    Y = Data[folds[[2]],(p2+1):p]
  }else{
    S1 <- cov(Data)
    S3 <- cov(Data)
    X = Data[,1:p1]
    Y = Data[,(p2+1):p]
  }
  sigma0hat1 <- S1 * Mask
  sqx <- sqrtm(Sigmax)$B
  sqy <- sqrtm(Sigmay)$B
  apply(S1[1:p1, (p1+1):p], 1, max)
  ag <- sgca_init(A=S1, B=sigma0hat1, rho=0.5 * sqrt(log(p)/dim(X)[1]),
                  K=r ,nu=nu,trace=FALSE, maxiter = maxiter) ###needs to be changed to be a little more scalable
  
  
  ainit <- init_process(ag$Pi, r) 
  init <- gca_to_cca(ainit, S3, pp)
  print("Init done")
  initu<- init$u
  initv <- init$v


  if (is.null(lambdax)){
    ### do CV
    resultsx <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_function(X, Y, kfolds, initu, initv,
                                                  lambdax = .x, adaptive=adaptive)))

    
    # print best hyperparameters and corresponding RMSE
    best_hyperparams <- resultsx[which.min(resultsx$rmse), ]
    which_lambdax = which(abs(resultsx$rmse-min(resultsx$rmse))/(1e-6 + min(resultsx$rmse)) <0.05)
    lambdax = max(resultsx$param1[which_lambdax])
  }
  if (is.null(lambday)){
    ### do CV
    results <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_function(Y, X, kfolds, initv, initu,
                                                  lambdax = .x, adaptive=adaptive)))
    best_hyperparams <- results[which.min(results$rmse), ]
    which_lambday = which(abs(results$rmse-min(results$rmse))/(1e-6 + min(results$rmse)) <0.05)
    lambday = max(results$param1[which_lambday])
  }
  
  ufinal = adaptive_lasso(X, Y %*% initv, initu, adaptive=adaptive, lambdax, 
                          max.iter=5000, 
                       max_k = 10, verbose = FALSE, ZERO_THRESHOLD=1e-5)
  
  vfinal = adaptive_lasso(Y, X %*% initu, initv, adaptive=adaptive, lambday, max.iter=5000, 
                       max_k = 10, verbose = FALSE, ZERO_THRESHOLD=1e-5)
  a_estimate = rbind(ufinal$Uhat, vfinal$Uhat)
  a_estimate <- gca_to_cca(a_estimate, S3, pp)
  return(list( ufinal = a_estimate$u, vfinal = a_estimate$v,
               initu=initu, initv=initv,
               Uhat= ufinal$Uhat, 
               Vhat=vfinal$Uhat,
               lambdax=lambdax, 
               lambday=lambday,
               resultsx=resultsx,
               resultsy=results
         )) #### Not too bad
}

## Running initialization using convex relaxation

#(Sigma,Sigma0, lambda, rho, eta=0.001, nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE)
pipeline_thresholded_gradient <- function(Data, Mask, sigma0hat, r=2, nu=1, Sigmax, 
                                          Sigmay, maxiter.init=30, 
                                          lambda=NULL, k=NULL, kfolds=5,
                                          maxiter=2000, convergence=1e-3, eta=1e-3,
                                          param1=10^(seq(-4, 1, by = 1)),
                                          param2=c(20, 1000)){
  p1 <- dim(Sigmax)[1]
  p2 <- dim(Sigmay)[1]
  p <- p1 + p2;
  n <- nrow(Data)
  pp <- c(p1,p2);
  S = cov(Data)
  ag <- sgca_init(A=S, B=sigma0hat, rho=0.5 * sqrt(log(p)/n),
                  K=r ,nu=1,trace=FALSE, maxiter = maxiter.init) ###needs to be changed
  print("Done with initialization")
  ainit <- init_process(ag$Pi, r) 
  init <- gca_to_cca(ainit, S, pp)
  
  if (is.null(lambda) | is.null(k)){
    resultsx <- expand.grid(lambda = param1, k = param2) %>%
      mutate(rmse = map2_dbl(lambda, k, ~ cv_function_tgd(Data[, 1:p1], Data[, (p1+1):p], 
                                                          Mask, kfolds=5, ainit,
                                                          lambda = .x,
                                                          k = .y, r=r,
                                                          maxiter=maxiter, eta=eta, convergence=convergence)))
                                                          #X, Y, Mask, kfolds=5, ainit,lambda, r=2, k=20,  
                                                          #maxiter=1000, eta=0.001, convergence=1e-

    #print(resultsx)
    ###### (X, Y, Mask, kfolds=5, ainit, lambda, k=20)
  
    # print best hyperparameters and corresponding RMSE
    best_hyperparams <- resultsx[which.min(resultsx$rmse), ]
    which_lambdax = which(abs(resultsx$rmse-min(resultsx$rmse))/(1e-6  + min(resultsx$rmse)) <0.05)
    lambda = max(resultsx$lambda[which_lambdax])
    k = max(resultsx$k[which_lambdax])
    #print(c("selected", k, lambda))
  }
  final <- sgca_tgd(A=S, B=sigma0hat,
                    r=r, ainit,k=k, lambda = lambda, eta=eta,convergence=convergence,
                    maxiter=maxiter, plot=FALSE)
  a_estimate <- gca_to_cca(final, S, pp)
  return(list( ufinal = a_estimate$u, vfinal = a_estimate$v,
               initu=init$u, initv=init$v,
               final=final,
               lambda=lambda, 
               k=k,
               resultsx=resultsx
  ))
  
}

additional_checks <- function(X_train, Y_train, S=NULL, 
                              rank=2, kfolds=5, method.type="FIT_SAR_BIC"){
  p1 <- dim(X_train)[2]
  p2 <- dim(Y_train)[2]
  p <- p1 + p2;
  n <- nrow(X_train)
  pp <- c(p1,p2);
  if(is.null(S)){
    S = cov(cbind(X_train, Y_train))
  }
  
  if (method.type=="FIT_SAR_BIC"){
    method<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                           lambdaAseq=seq(from=20,to=0.02,length=20),
                           lambdaBseq=seq(from=20,to=0.02,length=20),
                           max.iter=100,conv=10^-2,
                           selection.criterion=1,n.cv=5)
    a_estimate = rbind(method$uhat, method$vhat)
    
  }
  if(method.type=="FIT_SAR_CV"){
    method<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                          lambdaAseq=seq(from=0.2,to=0.02,length=10),
                          lambdaBseq=seq(from=0.2,to=0.02,length=10),
                          max.iter=100,conv=10^-2, selection.criterion=2, n.cv=5)
    a_estimate = rbind(method$uhat, method$vhat)
    
  }
  if (method.type=="Witten_Perm"){
    Witten_Perm <- CCA.permute(x=X_train,z=Y_train,typex="standard",typez="standard", nperms=50)
    method<-CCA(x=X_train, z=Y_train, typex="standard",typez="standard",K=rank,
                         penaltyx=Witten_Perm$bestpenaltyx,penaltyz=Witten_Perm$bestpenaltyz,trace=F)
    a_estimate = rbind(method$u, method$v)
  }
  if(method.type=="Witten.CV"){
    Witten_CV<-Witten.CV(X=X_train,Y=Y_train,n.cv=5,lambdax=matrix(seq(from=0,to=1,length=20),nrow=1),
                         lambday=matrix(seq(from=0,to=1,length=20),nrow=1))
    method <-CCA(x=X_train,z=Y_train,typex="standard",typez="standard",
                 K=rank,penaltyx=Witten_CV$lambdax.opt,penaltyz=Witten_CV$lambday.opt,trace=F)
    a_estimate = rbind(method$u, method$v)
    
  }
  if(method.type=="Waaijenborg-Author"){
    method<-Waaijenborg(X=X_train,Y=Y_train,
                        lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                        lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                        rank=rank,selection.criterion=1)
    a_estimate = rbind(method$vhat, method$uhat)
    
  }
  if(method.type=="Waaijenborg-CV"){
    method<-Waaijenborg(X=X_train,
                        Y=Y_train,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                        lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                        rank=rank, selection.criterion=2)
    a_estimate = rbind(method$vhat, method$uhat)
    
  }
  if(method.type=="SCCA_Parkhomenko"){
    method<- SCCA_Parkhomenko(x.data=X_train, y.data=Y_train, Krank=rank)
    a_estimate = rbind(method$uhat, method$vhat)
    
  }
  if(method.type=="Canonical Ridge-Author"){
    RCC_cv<-estim.regul_crossvalidation(X_train,Y_train,n.cv=5)
    method<-rcc(X_train,Y_train, RCC_cv$lambda1.optim, RCC_cv$lambda2.optim)
    a_estimate = rbind(method$xcoef[,1:rank], method$ycoef[,1:rank])
    
    
  }
  a_estimate <- gca_to_cca(a_estimate, S, pp)
  return(a_estimate)
}


