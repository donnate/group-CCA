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

# LOAD Functions


source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source("src/cca_bis.R")
### Under development

numCores <- detectCores()

vfn <- function(x){
  ifelse(x=="x", 1, -1)
}

library(foreach)    # install.packages('foreach')
library(caret)      # install.packages('caret', dependencies = c("Depends", "Suggests"))
library(doParallel) # install.packages('doParallel')
registerDoParallel(makeCluster(4)) # Use 4 cores for parallel CV




simple_experiment <- function(n, p, q, sigma, k, effect_size = 2,
                              sigma_noise=0.1, power=1,
                              lambdaA1seq= c(0.0001, 0.001, 0.01, 0.1, 1),
                              lambdaA2seq= c(0.0001, 0.001, 0.01, 0.1, 1),
                              conv=10^{-3}, max.iter=200,
                              egonet_size=3, n.cv=3,
                              mode = c("ECOS", "CD")){
  ## This is a simple experiment to see how CCA (and variants) can perform in
  ## different settings.
  ## n      :   number of observations
  ## p      :   number of covariates
  ## sigma  :   noise
  ## k      :   number of canonical vectors
  
  ### Create some sort of graph
  
  ### set 
  # lambdaA2 = lambdaA1/(8 * 12) 
  G <- sample_pa(p, power = power)
  E = data.frame(as_edgelist(G))
  colnames(E)  = c("x", "y")
  E["e"] = 1:nrow(E)
  E = pivot_longer(E, cols=-c("e"))
  E["fill_value"] = sapply(E$name, vfn)
  D = pivot_wider(E, id_cols =c("e"), names_from = "value", values_from  =  fill_value)
  D[is.na(D)] <- 0
  D = as.matrix(D[, 2:(p+1)])
  
  indices <- sample(1:p, k)
  X = matrix(rnorm(p * n, mean=0, sd=sigma_noise), n, p)
  Y = matrix(rnorm(q*n,  mean=0, sd=sigma_noise), n, q)
  source = matrix(0, k, p)
  colors = matrix(0, k, p)
  colors_l = rep(0, p)
  trueA = matrix(0, p, k)
  trueB = matrix(0, q, k)
  for (i in 1:k){
    print(i)
    idx = indices[i]
    subg <- ego(G, order=egonet_size, nodes = idx, 
                mode = "all", mindist = 0)[[1]]
    nodes_in_network <- as.numeric(subg)
    colors_l[nodes_in_network] = i
    colors[i, nodes_in_network]  = i
    source[i, nodes_in_network] = 1
    mean_value <- sapply(1:n, function(i){rnorm(1, mean=effect_size, sd=sigma)})
    X[, nodes_in_network] <- mean_value + X[, nodes_in_network] 
    Y[, i] =Y[, i] + (0.9)^(i) * X[, nodes_in_network] %*% matrix(rep(1, length(nodes_in_network)), nrow=length(nodes_in_network), ncol=1)
  }
  #heatmap(source)
  
  
  df.x <- data.frame(X) %>% mutate_all(~(scale(.) %>% as.vector))
  df.y <- data.frame(Y) %>% mutate_all(~(scale(.) %>% as.vector))
  X <- as.matrix(df.x)
  Y <- as.matrix(df.y)
  
  for ( i in 1:k){
    trueA[source[i, ],i] = 1/sqrt(sum(X[, source[i,]]^2))
    trueB[i,i] = 1/ abs(Y[i,i])
  }
  
  
  # Store results Estimation Accuracy
  Nsim = 1
  MSEa<-matrix(NA,ncol=10,nrow=Nsim)
  colnames(MSEa)<-c("Regular-CCA",  "SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author",
                   "genCCA")
  MSEb<-matrix(NA,ncol=10,nrow=Nsim)
  colnames(MSEb)<-c("Regular-CCA",  "SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author",
                     "genCCA")
  
  # Store results Sparsity Recognition Performance
  TPRa<-matrix(NA,ncol=10,nrow=Nsim)
  colnames(TPRa)<-c("Regular-CCA",  "SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author",
                    "genCCA")
  TPRb<-matrix(NA,ncol=10,nrow=Nsim)
  colnames(TPRb)<-c("Regular-CCA",  "SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author",
                    "genCCA")
  TNRa<-matrix(NA,ncol=10,nrow=Nsim)
  colnames(TNRa)<-c("Regular-CCA",  "SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author",
                    "genCCA")
  TNRb<-matrix(NA,ncol=10,nrow=Nsim)
  colnames(TNRb)<-c("Regular-CCA",  "SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author",
                    "genCCA")
  l2loss<-matrix(NA,ncol=10,nrow=Nsim)
  colnames(l2loss)<-c("Regular-CCA",  "SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author",
                    "genCCA")
  l2loss_test<-matrix(NA,ncol=10,nrow=Nsim)
  colnames(l2loss_test)<-c("Regular-CCA",  "SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author",
                    "genCCA")
  
  
  
  #### divide data into train and test
  train  = sample(1:n, size=as.integer(0.8 *n), replace = FALSE)
  test = setdiff(1:n, train)
  X_train = X[train, ]
  X_test = X[test, ]
  Y_train = Y[train, ]
  Y_test = Y[test, ]
  
  # if (p<n){
  #   cc_results <- cancor(X_train, Y_train, xcenter=FALSE, ycenter=FALSE)
  #   MSEa[1,1]<-principal_angles(trueA, cc_results$xcoef)$angles[1,1]
  #   MSEb[1,1]<-principal_angles(trueB, cc_results$ycoef)$angles[1,1]
  #   cancors <- cc_results$cor
  #   TPRa[1,1]<-TPR(trueA,cc_results$xcoef[,order(cancors,decreasing=T)])
  #   TPRb[1,1]<-TPR(trueB,cc_results$xcoef[,order(cancors,decreasing=T)])
  #   TNRa[1,1]<-TNR(trueA,cc_results$xcoef[,order(cancors,decreasing=T)])
  #   TNRb[1,1]<-TNR(trueB,cc_results$xcoef[,order(cancors,decreasing=T)])
  # }
  # 

  FIT_SAR_BIC<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                         lambdaAseq=seq(from=0.2,to=0.02,length=10),
                         lambdaBseq=seq(from=0.2,to=0.02,length=10),
                         max.iter=100,conv=10^-2,
                         selection.criterion=1,n.cv=5)
  MSEa[1,2]<-principal_angles(trueA,FIT_SAR_BIC$ALPHA)$angles[1,1]
  MSEb[1,2]<-principal_angles(trueB,FIT_SAR_BIC$BETA)$angles[1,1]
  cancors_SAR_BIC <- FIT_SAR_BIC$cancors
  TPRa[1,2]<-TPR(trueA,FIT_SAR_BIC$ALPHA[,order(cancors_SAR_BIC,decreasing=T)])
  TPRb[1,2]<-TPR(trueB,FIT_SAR_BIC$BETA[,order(cancors_SAR_BIC,decreasing=T)])
  TNRa[1,2]<-TNR(trueA,FIT_SAR_BIC$ALPHA[,order(cancors_SAR_BIC,decreasing=T)])
  TNRb[1,2]<-TNR(trueB,FIT_SAR_BIC$BETA[,order(cancors_SAR_BIC,decreasing=T)])
  l2loss[1,2] <- sqrt(mean((X_train %*% FIT_SAR_BIC$ALPHA - Y_train %*% FIT_SAR_BIC$BETA)^2))
  l2loss_test[1,2] <- sqrt(mean((X_test %*% FIT_SAR_BIC$ALPHA - Y_test %*% FIT_SAR_BIC$BETA)^2))
  
  FIT_SAR_CV<-SparseCCA(X=X,Y=Y,rank=rank,
                        lambdaAseq=seq(from=0.2,to=0.02,length=10),
                        lambdaBseq=seq(from=0.2,to=0.02,length=10),
                        max.iter=100,conv=10^-2, selection.criterion=2, n.cv=5)
  MSEa[1,3]<-principal_angles(trueA,FIT_SAR_CV$ALPHA)$angles[1,1]
  MSEb[1,3]<-principal_angles(trueB,FIT_SAR_CV$BETA)$angles[1,1]
  cancors_SAR_CV <- FIT_SAR_CV$cancors
  TPRa[1,3]<-TPR(trueA,FIT_SAR_CV$ALPHA[,order(cancors_SAR_CV,decreasing=T)])
  TPRb[1,3]<-TPR(trueB,FIT_SAR_CV$BETA[,order(cancors_SAR_CV,decreasing=T)])
  TNRa[1,3]<-TNR(trueA,FIT_SAR_CV$ALPHA[,order(cancors_SAR_CV,decreasing=T)])
  TNRb[1,3]<-TNR(trueB,FIT_SAR_CV$BETA[,order(cancors_SAR_CV,decreasing=T)])
  l2loss[1,3] <- sqrt(mean((X_train %*% FIT_SAR_CV$ALPHA - Y_train %*% FIT_SAR_CV$BETA)^2))
  l2loss_test[1, 3] <- sqrt(mean((X_test %*% FIT_SAR_CV$ALPHA - Y_test %*% FIT_SAR_CV$BETA)^2))
  
  
  Witten_Perm<-CCA.permute(x=X,z=Y,typex="standard",typez="standard",nperms=50)
  WittenSCCA_Perm<-CCA(x=X,z=Y,typex="standard",typez="standard",K=rank,penaltyx=Witten_Perm$bestpenaltyx,penaltyz=Witten_Perm$bestpenaltyz,trace=F)
  MSEa[1,4]<-principal_angles(trueA,WittenSCCA_Perm$u)$angles[1,1]
  MSEb[1,4]<-principal_angles(trueB,WittenSCCA_Perm$v)$angles[1,1]
  cancors_witten <- WittenSCCA_Perm$cors
  TPRa[1,4]<-TPR(trueA,WittenSCCA_Perm$u[,order(cancors_witten,decreasing=T)])
  TPRb[1,4]<-TPR(trueB,WittenSCCA_Perm$v[,order(cancors_witten,decreasing=T)])
  TNRa[1,4]<-TNR(trueA,WittenSCCA_Perm$u[,order(cancors_witten,decreasing=T)])
  TNRb[1,4]<-TNR(trueB,WittenSCCA_Perm$v[,order(cancors_witten,decreasing=T)])
  l2loss[1,4] <- sqrt(mean((X_train %*% WittenSCCA_Perm$u - Y_train %*% WittenSCCA_Perm$v)^2))
  l2loss_test[1,4] <- sqrt(mean((X_test %*% WittenSCCA_Perm$u - Y_test %*% WittenSCCA_Perm$v)^2))
  
  Witten_CV<-Witten.CV(X=X,Y=Y,n.cv=5,lambdax=matrix(seq(from=0,to=1,length=20),nrow=1),lambday=matrix(seq(from=0,to=1,length=20),nrow=1))
  WittenSCCA_CV<-CCA(x=X,z=Y,typex="standard",typez="standard",K=rank,penaltyx=Witten_CV$lambdax.opt,penaltyz=Witten_CV$lambday.opt,trace=F)
  MSEa[1,5]<-principal_angles(trueA,WittenSCCA_CV$u)$angles[1,1]
  MSEb[1,5]<-principal_angles(trueB,WittenSCCA_CV$v)$angles[1,1]
  cancors_witten <- WittenSCCA_CV$cors
  TPRa[1,5]<-TPR(trueA,WittenSCCA_CV$u[,order(cancors_witten,decreasing=T)])
  TPRb[1,5]<-TPR(trueB,WittenSCCA_CV$v[,order(cancors_witten,decreasing=T)])
  TNRa[1,5]<-TNR(trueA,WittenSCCA_CV$u[,order(cancors_witten,decreasing=T)])
  TNRb[1,5]<-TNR(trueB,WittenSCCA_CV$v[,order(cancors_witten,decreasing=T)])
  l2loss[1,5] <- sqrt(mean((X_train %*% WittenSCCA_CV$u - Y_train %*% WittenSCCA_CV$v)^2))
  l2loss_test[1,5] <- sqrt(mean((X_test %*% WittenSCCA_CV$u - Y_test %*% WittenSCCA_CV$v)^2))
  
  
  Waaijenborg_Delta<-Waaijenborg(X=X,Y=Y,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),rank=rank,selection.criterion=1)
  MSEa[1,6]<-principal_angles(trueA,Waaijenborg_Delta$vhat)$angles[1,1]
  MSEb[1,6]<-principal_angles(trueB,Waaijenborg_Delta$uhat)$angles[1,1]
  cancors_Waaijenborg <- Waaijenborg_Delta$cancors
  TPRa[1,6]<-TPR(trueA, Waaijenborg_Delta$vhat[,order(cancors_Waaijenborg,decreasing=T)])
  TPRb[1,6]<-TPR(trueB, Waaijenborg_Delta$uhat[,order(cancors_Waaijenborg,decreasing=T)])
  TNRa[1,6]<-TNR(trueA, Waaijenborg_Delta$vhat[,order(cancors_Waaijenborg,decreasing=T)])
  TNRb[1,6]<-TNR(trueB, Waaijenborg_Delta$uhat[,order(cancors_Waaijenborg,decreasing=T)])
  l2loss[1,6] <- sqrt(mean((X_train %*% Waaijenborg_Delta$vhat - Y_train %*% Waaijenborg_Delta$uhat) ^2))
  l2loss_test[1,6] <- sqrt(mean((X_test %*% Waaijenborg_Delta$vhat  - Y_test %*% Waaijenborg_Delta$uhat )^2))
  
  Waaijenborg_Test<-Waaijenborg(X=X,Y=Y,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),rank=rank,selection.criterion=2)
  MSEa[1,7]<-principal_angles(trueA,Waaijenborg_Test$vhat)$angles[1,1]
  MSEb[1,7]<-principal_angles(trueB,Waaijenborg_Test$uhat)$angles[1,1]
  cancors_Waaijenborg <- Waaijenborg_Test$cancors
  TPRa[1,7]<-TPR(trueA, Waaijenborg_Test$vhat[,order(cancors_Waaijenborg,decreasing=T)])
  TPRb[1,7]<-TPR(trueB, Waaijenborg_Test$uhat[,order(cancors_Waaijenborg,decreasing=T)])
  TNRa[1,7]<-TNR(trueA, Waaijenborg_Test$vhat[,order(cancors_Waaijenborg,decreasing=T)])
  TNRb[1,7]<-TNR(trueB, Waaijenborg_Test$uhat[,order(cancors_Waaijenborg,decreasing=T)])
  l2loss[1,7] <- sqrt(mean((X_train %*% Waaijenborg_Test$vhat - Y_train %*% Waaijenborg_Test$uhat) ^2))
  l2loss_test[1,7] <- sqrt(mean((X_test %*% Waaijenborg_Test$vhat  - Y_test %*% Waaijenborg_Test$uhat )^2))
  
  Parkhomenko_SCCA<-SCCA_Parkhomenko(x.data=X,y.data=Y,Krank=rank)
  MSEa[1,8]<-principal_angles(trueA,Parkhomenko_SCCA$a)$angles[1,1]
  MSEb[1,8]<-principal_angles(trueB,Parkhomenko_SCCA$b)$angles[1,1]
  cancors_Parkhomenko <- Parkhomenko_SCCA$cancor
  TPRa[1,8]<-TPR(trueA, Parkhomenko_SCCA$a[,order(cancors_Parkhomenko,decreasing=T)])
  TPRb[1,8]<-TPR(trueB, Parkhomenko_SCCA$b[,order(cancors_Parkhomenko ,decreasing=T)])
  TNRa[1,8]<-TNR(trueA, Parkhomenko_SCCA$a[,order(cancors_Parkhomenko ,decreasing=T)])
  TNRb[1,8]<-TNR(trueB, Parkhomenko_SCCA$b[,order(cancors_Parkhomenko ,decreasing=T)])
  l2loss[1,8] <- sqrt(mean((X_train %*% Parkhomenko_SCCA$a - Y_train %*%Parkhomenko_SCCA$b) ^2))
  l2loss_test[1,8] <- sqrt(mean((X_test %*% Parkhomenko_SCCA$a  - Y_test %*% Parkhomenko_SCCA$b)^2))
  
  RCC_cv<-estim.regul_crossvalidation(X,Y,n.cv=5)
  RCCA<-rcc(X,Y,RCC_cv$lambda1.optim,RCC_cv$lambda2.optim)
  MSEa[1,9]<-principal_angles(trueA,RCCA$xcoef)$angles[1,1]
  MSEb[1,9]<-principal_angles(trueB,RCCA$ycoef)$angles[1,1]
  cancors_RCCA<- RCCA$cor
  TPRa[1,9]<-TPR(trueA, RCCA$xcoef[,order(cancors_RCCA,decreasing=T)[1:rank]])
  TPRb[1,9]<-TPR(trueB, RCCA$ycoef[,order(cancors_RCCA,decreasing=T)[1:rank]])
  TNRa[1,9]<-TNR(trueA, RCCA$xcoef[,order(cancors_RCCA,decreasing=T)[1:rank]])
  TNRb[1,9]<-TNR(trueB, RCCA$ycoef[,order(cancors_RCCA,decreasing=T)[1:rank]])
  l2loss[1,9] <- sqrt(mean((X_train %*% RCCA$xcoef - Y_train %*%RCCA$ycoef) ^2))
  l2loss_test[1,9] <- sqrt(mean((X_test %*%  RCCA$xcoef - Y_test %*% RCCA$ycoef)^2))
  

  genCCA.results <-  genCCA.CV(X, Y, D, rank, n.cv=3,
                      lambda1seq=c(0, 0.01, 0.1, 1, 10, 100),
                      lambda2seq=c(0, 0.01, 0.1, 1, 10, 100))
  MSEa[1,10]<-principal_angles(trueA,genCCA.results$xcoef)$angles[1,1]
  MSEb[1,10]<-principal_angles(trueB,genCCA.results$ycoef)$angles[1,1]
  cancors = diag(t(X %*% genCCA.results$xcoef) %*% (Y %*% genCCA.results$ycoef))
  TPRa[1,10]<-TPR(trueA, genCCA.results$xcoef[,order(cancors,decreasing=T)[1:rank]])
  TPRb[1,10]<-TPR(trueB, genCCA.results$ycoef[,order(cancors,decreasing=T)[1:rank]])
  TNRa[1,10]<-TNR(trueA, genCCA.results$xcoef[,order(cancors,decreasing=T)[1:rank]])
  TNRb[1,10]<-TNR(trueB, genCCA.results$ycoef[,order(cancors,decreasing=T)[1:rank]])
  l2loss[1,10] <- sqrt(mean((X_train %*% genCCA.results$xcoef - Y_train %*%genCCA.results$ycoef) ^2))
  l2loss_test[1,10] <- sqrt(mean((X_test %*%  genCCA.results$xcoef - Y_test %*% genCCA.results$ycoef)^2))

  ###
  
  ### Select the best out of all the different options
  ### selected lambdas = which.min(cv.results)
  #### look at the predicted correlations
  
  #### Assess the quality of the reconstruction


  
  
  return(list(X=X, Y=Y, Z=Z, 
              gencca_results = gencca_results.final,
              cc_results=cc_results, cv.results= cv.results ))
}




process_results <- function(gencca_results.final, source){
  
  ### postprocessing: 
  
  t.test(gencca_results.final$ALPHA[which(source[1,]==1),1], gencca_results.final$ALPHA[which(source[1,]==0),1])
  ggplot(coefs)+
    geom_point(aes(x=index, y=X1, colour= colour1, shape=colour2))
  
  ggplot(coefs)+
    geom_point(aes(x=index, y=X2, colour= colour1, shape=colour2))
  
  df.all = cbind(df.x , df.y)
  colnames(df.all) = c(sapply(1:p, function(x){paste0("X", x)}), sapply(1:q, function(x){paste0("Y", x)}))
  
  
  test = genCCA(as.matrix(df.x),  as.matrix(df.y), D, NULL,
                lambdaA1=0,
                lambdaB1=NULL,
                lambdaA2=0.,
                lambdaB2=0.,
                rank=k,
                A.initial=NULL,B.initial=NULL,
                max.iter=100,conv=10^-4)
  train_test_split <- initial_split(df.all, prop = 0.9)
  df_train <- training(train_test_split)
  df_test <- testing(train_test_split)
  validation_data =  mc_cv(df_train, prop = 0.9, times = 10)
  coefs_long = pivot_longer(coefs, cols=-c("index", sapply(1:k, FUN = function(k){paste0("colour", k)})))
  coefs_long
  
  ggplot(coefs_long %>% filter(name %in% c("X1", "X2", "X3", "X4")))+
    geom_point(aes(x=index, y=value, shape= name))
}




