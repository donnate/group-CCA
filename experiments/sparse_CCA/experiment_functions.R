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
setwd("experiments/alternative_methods/GCA")
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

cv_function <- function(X, Y, 
                        kfolds=10, initu, initv,
                        lambdax, adaptive=TRUE, normalize=FALSE,
                        criterion="prediction") {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- numeric(length = kfolds)
  nnz  <- numeric(length = kfolds)
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
    model <- adaptive_lasso(X_train, Y_train %*% initv, initu, adaptive=adaptive, 
                         lambdax, 
                         max.iter=5000, 
                         max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
    
    # make predictions on validation data
    # compute RMSE on validation data
    ##### Normalize Uhat
    if (normalize == FALSE){
      if (criterion=="prediction"){
        rmse[i] <- sum((X_val %*% model$Uhat - Y_val%*% initv)^2)
        if (norm(model$Uhat) == 0){ ####prevents selecting values that would make everything 0
          rmse[i] <- 1e8
        }
      }else{
        rmse[i] <- sum(abs(cor(X_val %*% model$Uhat, Y_val%*% initv)))
      }
      nnz[i] <- sum(apply(model$Uhat^2, 1, sum) >1e-4)
    }else{
      sol <- gca_to_cca(rbind(model$Uhat, initv), 
                        cov(rbind(X_val, Y_val)), pp)
      if (criterion=="prediction"){
      rmse[i] <- sum((X_val %*% sol$u - Y_val%*% initv)^2)
      }else{
        rmse[i] <- sum(abs(cor(X_val %*%sol$u, Y_val%*% initv)))
      }
      nnz[i] <- sum(apply(sol$u^2, 1, sum) >1e-4)
     }
    },
    error = function(e) {
      # If an error occurs, assign NA to the result
      rmse[i] <- NA
    })
  }
  
  # return mean RMSE across folds
  if (mean(is.na(rmse)) == 1){
      return(1e8)
   }else{
  return(median(rmse, na.rm=TRUE))
   }
}


cv_function_tgd <- function(X, Y, Mask, kfolds=5, ainit,
                        lambda, r=2, k=20,  maxiter=1000, eta=0.001,
                        convergence=1e-3, normalize=FALSE, criterion="prediction") {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- numeric(length = kfolds)
  p1 <- dim(X)[2]
  p2 <- dim(Y)[2]
  p <- p1 + p2;
  n <- nrow(X)
  pp <- c(p1,p2);
  S0 = cov(cbind(X, Y))
  
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
    if (normalize == FALSE){
      if (criterion=="prediction"){
        rmse[i] <- sum((X_val %*% final$u - Y_val%*% initv)^2)
      }else{
        rmse[i] <- sum(abs(cor(X_val %*% final$u, Y_val%*% initv)))
      }
    }else{
      sol <- gca_to_cca(rbind(final$u, initv), 
                        cov(rbind(X_val, Y_val)), pp)
      if (criterion=="prediction"){
        rmse[i] <- sum((X_val %*% sol$u - Y_val%*% initv)^2)
      }else{
        rmse[i] <- sum(abs(cor(X_val %*%sol$u, Y_val%*% initv)))
      }
    }
    },
    error = function(e) {
      # If an error occurs, assign NA to the result
      rmse[i] <- NA
    })
  }
  
  # return mean RMSE across folds
  if (mean(is.na(rmse)) == 1){
      return(1e8)
   }else{
  return(median(rmse, na.rm=TRUE))
   }
}

preselection <-function(Data, CorrelationMat, p1, r, alpha){
  p = ncol(Data)
  n=  nrow(Data)
  p2 =  p-p1
  t = apply(CorrelationMat -diag(diag(CorrelationMat)), 1, 
            function(x){max(x^2)})
  J = order(-t)[1: ceiling(alpha  * n/(r))]
  set_u = J[which(J <= p1)]
  set_v = J[which(J > p1)]
  t=CCA::cc(as.matrix(Data[,set_u]), as.matrix(Data[, set_v]))
  Uhat = matrix(0, p, r)
  Uhat[set_u, ] =  t$xcoef[,1:r]
  Uhat[set_v, ] =  t$ycoef[,1:r]
  return(Uhat)
}


pipeline_adaptive_lasso <- function(Data, Mask, sigma0hat, r, nu=1, Sigmax, 
                                    Sigmay, maxiter=30, lambdax=NULL, lambday=NULL,
                                    adaptive=TRUE, kfolds=5, param1=10^(seq(-4, 2, by = 0.25)),
                                    create_folds=TRUE, init ="Fantope", normalize=FALSE, alpha=0.5,
                                    criterion="prediction",
                                    fantope_solution=NULL){
  
  ### data splitting procedure 3 folds
  #maxiter=100;  lambdax=NULL;
 #adaptive=TRUE; kfolds=5;  param1=10^(seq(-5, 1, by = 0.5));
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
    X1 = Data[folds[[2]],1:p1]
    Y1= Data[folds[[2]],(p2+1):p]
    sigma0hat1 <- S1 * Mask
  }else{
    S1 <- cov(as.matrix(Data))
    S3 <- S1
    X = Data[,1:p1]
    Y = Data[,(p2+1):p]
    X1 = Data[,1:p1]
    Y1 = Data[,(p2+1):p]
    sigma0hat1 <- sigma0hat
  }
 
  if (init == "Fantope"){
    if (is.null(fantope_solution)){
      ag <- sgca_init(A=S1, B=sigma0hat1, rho=0.5 * sqrt(log(p)/dim(X)[1]),
                      K=r ,nu=nu,trace=FALSE, maxiter = maxiter) ###needs to be changed to be a little more scalable
      ainit <- init_process(ag$Pi, r) 
    }else{
      ainit <- fantope_solution
    }
      
  
  }else{
    CorrelationMatrix =  diag(1/sqrt(diag(example$S))) %*% example$S %*% diag(1/sqrt(diag(example$S)))
     ainit= preselection(example$Data, CorrelationMatrix, p1, r, alpha)
  }

  init <- gca_to_cca(ainit, S3, pp)
  print("Init done")
  initu<- init$u
  initv <- init$v


  if (is.null(lambdax)){
    ### do CV
    resultsx <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_function(X, Y, kfolds, initu, initv,
                                                  lambdax = .x, adaptive=adaptive,
                                                   normalize=normalize,
                                                  criterion=criterion)))

    
    # print best hyperparameters and corresponding RMSE
    best_hyperparams <- resultsx[which.min(resultsx$rmse), ]
    which_lambdax = which(abs(resultsx$rmse-min(resultsx$rmse))/(1e-6 + min(resultsx$rmse)) <0.05)
    lambdax = max(resultsx$param1[which_lambdax])
    resultsx =NULL
  }
  if (is.null(lambday)){
    ### do CV
    results <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_function(Y, X, kfolds, initv, initu,
                                                  lambdax = .x, adaptive=adaptive,
                                                  criterion=criterion)))
    best_hyperparams <- results[which.min(results$rmse), ]
    which_lambday = which(abs(results$rmse-min(results$rmse))/(1e-6 + min(results$rmse)) <0.05)
    lambday = max(results$param1[which_lambday])
    results =NULL
  }
  
  ufinal = adaptive_lasso(X, Y %*% initv, initu, adaptive=adaptive, lambdax, 
                          max.iter=5000, 
                       max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
  
  vfinal = adaptive_lasso(Y, X %*% initu, initv, adaptive=adaptive, lambday, max.iter=5000, 
                       max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
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

pipeline_alternating_lasso <- function(Data, Mask, sigma0hat, r, nu=1, Sigmax, 
                                       Sigmay, maxiter=30, lambdax=NULL, lambday=NULL,
                                       adaptive=TRUE, kfolds=5, param1=10^(seq(-4, 2, by = 0.25)),
                                       create_folds=TRUE, init ="Fantope", normalize=FALSE, alpha=0.5,
                                       criterion="prediction",
                                       fantope_solution=NULL){
  
  ### Replaces the initialization step by an rcc
  p1 <- dim(Sigmax)[1]
  p2 <- dim(Sigmay)[1]
  p <- p1 + p2;
  n <- nrow(Data)
  pp <- c(p1,p2);
  RCC_cv<-estim.regul_crossvalidation(Data[,1:p1], Data[,(p2+1):p],
                                      n.cv=5)
  method<-rcc(Data[,1:p1], Data[,(p2+1):p], 
              RCC_cv$lambda1.optim, RCC_cv$lambda2.optim)
  U0 = method$xcoef[,1:r]; ### Initial values
  V0 = method$ycoef[,1:r];
  X = Data[,1:p1]
  Y = Data[,(p2+1):p]
  Uinit = U0
  Vinit = V0
  converged = FALSE
  it = 0
  while(converged==FALSE & it <10){
    resultsx <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_function(X, Y, kfolds, U0, V0,
                                                  lambdax = .x, adaptive=adaptive,
                                                  normalize=normalize,
                                                  criterion=criterion)))
    
    
    # print best hyperparameters and corresponding RMSE
    best_hyperparams <- resultsx[which.min(resultsx$rmse), ]
    which_lambdax = which(abs(resultsx$rmse-min(resultsx$rmse))/(1e-6 + min(resultsx$rmse)) <0.05)
    lambdax = max(resultsx$param1[which_lambdax])
    Unew = adaptive_lasso(X, Y %*% V0, U0, adaptive=adaptive, lambdax, 
                          max.iter=5000, 
                          max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
    
    results <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_function(Y, X, kfolds,  V0, Unew$Uhat,
                                                  lambdax = .x, adaptive=adaptive,
                                                  criterion=criterion)))
    best_hyperparams <- results[which.min(results$rmse), ]
    which_lambday = which(abs(results$rmse-min(results$rmse))/(1e-6 + min(results$rmse)) <0.05)
    lambday = max(results$param1[which_lambday])
    
    Vnew = adaptive_lasso(Y, X %*% Unew$Uhat, V0, adaptive=adaptive, lambday, max.iter=5000, 
                            max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
    
    converged = (mean((Unew$Uhat-U0)^2) <1e-4) & (mean((Vnew$Uhat-V0)^2) <1e-4)
    it = it + 1
    U0  = Unew$Uhat
    V0 = Vnew$Uhat
    print(it)
    
  }
  
  return(list( ufinal = Unew, vfinal = Vnew,
               initu=Uinit, initv=Vinit,
               lambda.x=lambdax, 
               lambda.y=lambday
  ))
  
  
  
}

## Running initialization using convex relaxation

#(Sigma,Sigma0, lambda, rho, eta=0.001, nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE)
pipeline_thresholded_gradient <- function(Data, Mask, sigma0hat, r=2, nu=1, Sigmax, 
                                          Sigmay, maxiter.init=30, 
                                          lambda=NULL, k=NULL, kfolds=5,
                                          maxiter=2000, convergence=1e-3, eta=1e-3,
                                          param1=10^(seq(-4, 1, by = 1)),
                                          param2=c(20, 1000), init="Fantope",
                                          normalize=FALSE,criterion="prediction",
                                          fantope_solution=NULL){
  p1 <- dim(Sigmax)[1]
  p2 <- dim(Sigmay)[1]
  p <- p1 + p2;
  n <- nrow(Data)
  pp <- c(p1,p2);
  S = cov(Data)

  if (init == "Fantope"){
    if (is.null(fantope_solution)){
      ag <- sgca_init(A=S1, B=sigma0hat1, rho=0.5 * sqrt(log(p)/n),
                      K=r ,nu=nu,trace=FALSE, maxiter = maxiter) ###needs to be changed to be a little more scalable
      ainit <- init_process(ag$Pi, r) 
    }else{
      ainit <- fantope_solution
    }
  }else{
      RCC_cv<-estim.regul_crossvalidation(Data[,1:p1], Data[,(p2+1):p],
                                          n.cv=5)
      method<-rcc(Data[,1:p1], Data[,(p2+1):p], 
                    RCC_cv$lambda1.optim, RCC_cv$lambda2.optim)
      ainit= rbind(method$xcoef[,1:r], method$ycoef[,1:r])
  }

  init <- gca_to_cca(ainit, S, pp)
  
  if (is.null(lambda) | is.null(k)){
    resultsx <- expand.grid(lambda = param1, k = param2) %>%
      mutate(rmse = map2_dbl(lambda, k, ~ cv_function_tgd(Data[, 1:p1], Data[, (p1+1):p], 
                                                          Mask, kfolds=5, ainit,
                                                          lambda = .x,
                                                          k = .y, r=r,
                                                          maxiter=maxiter, eta=eta, convergence=convergence,
                                                          normalize=normalize,
                                                          criterion=criterion)))
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


evaluate_results <- function(Uhat, Vhat, example, name_method, overlapping_amount, 
                             thres = 0.0001, lambdax= NULL,lambday = NULL, it=1,
                             normalize_diagonal=TRUE,
                             criterion="prediction"){
  Uhat_tot = rbind(Uhat, Vhat)
  U_tot = rbind(example$u, example$v)
  p1 = ncol(example$X)
  p2 = ncol(example$Y)
  n = nrow(example$X)
  r = ncol(example$u)
  sparsity = length(which(apply(example$u^2, 1, sum)>0))
  silly_benchmark = subdistance(matrix(0, p1, 2), example$u)
  data.frame("method" = name_method,
             "exp" = it,
             "n" = n,
             "nnz" = nnz,
             "p1" = p1,
             "p2" = p2,
             "sparsity" = sparsity,
             "criterion" = criterion,
             "overlapping_amount" = overlapping_amount,
             "zero_benchmark" = silly_benchmark,
             "nb_discoveries" = sum(apply(Uhat_tot^2, 1, sum)>0),
             "nb_real_discoveries" = sum(apply(Uhat_tot^2, 1, sum)>thres),
             "param1" = lambdax,
             "param2" = lambday,
             "normalize_diagonal" = normalize_diagonal,
             "distance_tot" = subdistance(Uhat_tot, U_tot),
             "distance_U" = subdistance(Uhat, example$u),
             "distance_V" = subdistance(Vhat, example$v),
             "sinTheta_tot" = sinTheta(Uhat_tot, U_tot),
             "sinTheta_U" = sinTheta(Uhat, example$u),
             "sinTheta_V" = sinTheta(Vhat, example$v),
             "sinTheta_U_0" = sinTheta(matrix(0, p1, r), example$u),
             "sinTheta_V_0" = sinTheta(matrix(0, p2, r), example$v),
             "prediction_tot" = mean((example$X %*% Uhat - example$Y %*% Vhat)^2),
             "prediction_U" = mean((example$X %*% Uhat - example$X %*% example$u)^2),
             "prediction_V" = mean((example$Y %*% Vhat - example$Y %*% example$v)^2),
             "TPR" =TPR(apply(Uhat_tot^2, 1, sum), apply(U_tot^2, 1, sum), tol=thres),
             "TNR" = TNR(apply(Uhat_tot^2, 1, sum), apply(U_tot^2, 1, sum), tol=thres),
             "FPR" = FPR(apply(Uhat_tot^2, 1, sum), apply(U_tot^2, 1, sum), tol=thres),
             "FNR" = FPR(apply(U_tot^2, 1, sum),apply(Uhat_tot^2, 1, sum), tol=thres)
  )
}
