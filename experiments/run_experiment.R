#!/usr/bin/env Rscript
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source("src/cca_bis.R")
source("src/CCA_alaChao.R")

library(tidyverse)
library(PMA)
library(mvtnorm)
library(CCA)
args = commandArgs(trailingOnly=TRUE)

#### load(data)
name_experiment  = args[1]
type_experiment = args[2]
mysavedir = args[3]
rank =as.numeric( args[4])


load(paste0(mysavedir, name_experiment, '.RData'))

if (type_experiment == "others"){
  methods = c("Regular-CCA",  "SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author")
  MSEa<-matrix(NA,ncol=length(methods),nrow=Nsim)
  colnames(MSEa)<-methods
  MSEb<-matrix(NA,ncol=length(methods),nrow=Nsim)
  colnames(MSEb)<-methods
  # Store results Sparsity Recognition Performance
  TPRa<-matrix(NA,ncol=length(methods),nrow=Nsim)
  colnames(TPRa)<-methods
  TPRb<-matrix(NA,ncol=length(methods),nrow=Nsim)
  colnames(TPRb)<-methods
  TNRa<-matrix(NA,ncol=length(methods),nrow=Nsim)
  colnames(TNRa)<-methods
  TNRb<-matrix(NA,ncol=length(methods),nrow=Nsim)
  colnames(TNRb)<-methods
  l2loss<-matrix(NA,ncol=length(methods),nrow=Nsim)
  colnames(l2loss)<-methods
  l2loss_test<-matrix(NA,ncol=length(methods),nrow=Nsim)
  colnames(l2loss_test)<-methods
  
  results.x <-list()
  results.y <-list()
  recovered_corr <- list()
  
  
  # cc_results <- cancor(X_train, Y_train, xcenter=FALSE, ycenter=FALSE)
  # MSEa[1,1]<-principal_angles(trueA, cc_results$xcoef)$angles[1,1]
  # MSEb[1,1]<-principal_angles(trueB, cc_results$ycoef)$angles[1,1]
  # cancors <- cc_results$cor
  # TPRa[1,1]<-TPR(trueA[, order(true_corr,decreasing=T)],cc_results$xcoef[,order(cancors,decreasing=T)[1:rank]])
  # TPRb[1,1]<-TPR(trueB[, order(true_corr,decreasing=T)],cc_results$ycoef[,order(cancors,decreasing=T)[1:rank]])
  # TNRa[1,1]<-TNR(trueA[, order(true_corr,decreasing=T)],cc_results$xcoef[,order(cancors,decreasing=T)[1:rank]])
  # TNRb[1,1]<-TNR(trueB[, order(true_corr,decreasing=T)],cc_results$ycoef[,order(cancors,decreasing=T)[1:rank]])
  # l2loss[1,1] <- sqrt(mean((X_train %*% cc_results$xcoef[,order(cancors,decreasing=T)[1:rank]]- Y_train %*% cc_results$ycoef[,order(cancors,decreasing=T)[1:rank]])^2))
  # l2loss_test[1, 1] <- sqrt(mean((X_test %*% cc_results$xcoef[,order(cancors,decreasing=T)[1:rank]]- Y_test %*% cc_results$ycoef[,order(cancors,decreasing=T)[1:rank]])^2))
  # results.x[['Regular-CCA']] <- cc_results$xcoef
  # results.y[['Regular-CCA']] <- cc_results$ycoef
  # recovered_corr [['Regular-CCA']] <- cancors
  
  FIT_SAR_BIC<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                         lambdaAseq=seq(from=0.2,to=0.02,length=10),
                         lambdaBseq=seq(from=0.2,to=0.02,length=10),
                         max.iter=100,conv=10^-2,
                         selection.criterion=1,n.cv=5)
  MSEa[1,2]<-principal_angles(trueA,FIT_SAR_BIC$ALPHA)$angles[1,1]
  MSEb[1,2]<-principal_angles(trueB,FIT_SAR_BIC$BETA)$angles[1,1]
  cancors_SAR_BIC <- FIT_SAR_BIC$cancors
  TPRa[1,2]<-TPR(trueA[, order(true_corr,decreasing=T)],FIT_SAR_BIC$ALPHA[,order(cancors_SAR_BIC,decreasing=T)[1:rank]])
  TPRb[1,2]<-TPR(trueB[, order(true_corr,decreasing=T)],FIT_SAR_BIC$BETA[,order(cancors_SAR_BIC,decreasing=T)[1:rank]])
  TNRa[1,2]<-TNR(trueA[, order(true_corr,decreasing=T)],FIT_SAR_BIC$ALPHA[,order(cancors_SAR_BIC,decreasing=T)[1:rank]])
  TNRb[1,2]<-TNR(trueB[, order(true_corr,decreasing=T)],FIT_SAR_BIC$BETA[,order(cancors_SAR_BIC,decreasing=T)[1:rank]])
  l2loss[1,2] <- sqrt(mean((X_train %*% FIT_SAR_BIC$ALPHA - Y_train %*% FIT_SAR_BIC$BETA)^2))
  l2loss_test[1,2] <- sqrt(mean((X_test %*% FIT_SAR_BIC$ALPHA - Y_test %*% FIT_SAR_BIC$BETA)^2))
  results.x[['SAR-Author']] <- FIT_SAR_BIC$ALPHA
  results.y[['SAR-Author']] <- FIT_SAR_BIC$BETA
  recovered_corr [['SAR-Author']] <- cancors_SAR_BIC
  
  
  FIT_SAR_CV<-SparseCCA(X=X,Y=Y,rank=rank,
                        lambdaAseq=seq(from=0.2,to=0.02,length=10),
                        lambdaBseq=seq(from=0.2,to=0.02,length=10),
                        max.iter=100,conv=10^-2, selection.criterion=2, n.cv=5)
  MSEa[1,3]<-principal_angles(trueA,FIT_SAR_CV$ALPHA)$angles[1,1]
  MSEb[1,3]<-principal_angles(trueB,FIT_SAR_CV$BETA)$angles[1,1]
  cancors_SAR_CV <- FIT_SAR_CV$cancors
  TPRa[1,3]<-TPR(trueA[, order(true_corr,decreasing=T)],FIT_SAR_CV$ALPHA[,order(cancors_SAR_CV,decreasing=T)[1:rank]])
  TPRb[1,3]<-TPR(trueB[, order(true_corr,decreasing=T)],FIT_SAR_CV$BETA[,order(cancors_SAR_CV,decreasing=T)[1:rank]])
  TNRa[1,3]<-TNR(trueA[, order(true_corr,decreasing=T)],FIT_SAR_CV$ALPHA[,order(cancors_SAR_CV,decreasing=T)[1:rank]])
  TNRb[1,3]<-TNR(trueB[, order(true_corr,decreasing=T)],FIT_SAR_CV$BETA[,order(cancors_SAR_CV,decreasing=T)[1:rank]])
  l2loss[1,3] <- sqrt(mean((X_train %*% FIT_SAR_CV$ALPHA - Y_train %*% FIT_SAR_CV$BETA)^2))
  l2loss_test[1, 3] <- sqrt(mean((X_test %*% FIT_SAR_CV$ALPHA - Y_test %*% FIT_SAR_CV$BETA)^2))
  results.x[['SAR-CV']] <- FIT_SAR_CV$ALPHA
  results.y[['SAR-CV']] <- FIT_SAR_CV$BETA
  recovered_corr [['SAR-CV']] <- cancors_SAR_CV
  
  
  Witten_Perm<-CCA.permute(x=X,z=Y,typex="standard",typez="standard",nperms=50)
  WittenSCCA_Perm<-CCA(x=X,z=Y,typex="standard",typez="standard",K=rank,penaltyx=Witten_Perm$bestpenaltyx,penaltyz=Witten_Perm$bestpenaltyz,trace=F)
  MSEa[1,4]<-principal_angles(trueA,WittenSCCA_Perm$u)$angles[1,1]
  MSEb[1,4]<-principal_angles(trueB,WittenSCCA_Perm$v)$angles[1,1]
  cancors_witten <- WittenSCCA_Perm$cors
  TPRa[1,4]<-TPR(trueA[, order(true_corr,decreasing=T)], WittenSCCA_Perm$u[,order(cancors_witten,decreasing=T)[1:rank]])
  TPRb[1,4]<-TPR(trueB[, order(true_corr,decreasing=T)], WittenSCCA_Perm$v[,order(cancors_witten,decreasing=T)[1:rank]])
  TNRa[1,4]<-TNR(trueA[, order(true_corr,decreasing=T)], WittenSCCA_Perm$u[,order(cancors_witten,decreasing=T)[1:rank]])
  TNRb[1,4]<-TNR(trueB[, order(true_corr,decreasing=T)], WittenSCCA_Perm$v[,order(cancors_witten,decreasing=T)[1:rank]])
  l2loss[1,4] <- sqrt(mean((X_train %*% WittenSCCA_Perm$u - Y_train %*% WittenSCCA_Perm$v)^2))
  l2loss_test[1,4] <- sqrt(mean((X_test %*% WittenSCCA_Perm$u - Y_test %*% WittenSCCA_Perm$v)^2))
  results.x[['Witten-Author']] <- WittenSCCA_Perm$u
  results.y[['Witten-Author']] <- WittenSCCA_Perm$v
  recovered_corr [['Witten-Author']] <- cancors_witten
  
  Witten_CV<-Witten.CV(X=X,Y=Y,n.cv=5,lambdax=matrix(seq(from=0,to=1,length=20),nrow=1),lambday=matrix(seq(from=0,to=1,length=20),nrow=1))
  WittenSCCA_CV<-CCA(x=X,z=Y,typex="standard",typez="standard",K=rank,penaltyx=Witten_CV$lambdax.opt,penaltyz=Witten_CV$lambday.opt,trace=F)
  MSEa[1,5]<-principal_angles(trueA,WittenSCCA_CV$u)$angles[1,1]
  MSEb[1,5]<-principal_angles(trueB,WittenSCCA_CV$v)$angles[1,1]
  cancors_witten <- WittenSCCA_CV$cors
  TPRa[1,5]<-TPR(trueA[, order(true_corr,decreasing=T)],WittenSCCA_CV$u[,order(cancors_witten,decreasing=T)[1:rank]])
  TPRb[1,5]<-TPR(trueB[, order(true_corr,decreasing=T)],WittenSCCA_CV$v[,order(cancors_witten,decreasing=T)[1:rank]])
  TNRa[1,5]<-TNR(trueA[, order(true_corr,decreasing=T)],WittenSCCA_CV$u[,order(cancors_witten,decreasing=T)[1:rank]])
  TNRb[1,5]<-TNR(trueB[, order(true_corr,decreasing=T)],WittenSCCA_CV$v[,order(cancors_witten,decreasing=T)[1:rank]])
  l2loss[1,5] <- sqrt(mean((X_train %*% WittenSCCA_CV$u - Y_train %*% WittenSCCA_CV$v)^2))
  l2loss_test[1,5] <- sqrt(mean((X_test %*% WittenSCCA_CV$u - Y_test %*% WittenSCCA_CV$v)^2))
  results.x[['Witten-CV']] <- WittenSCCA_Perm$u
  results.y[['Witten-CV']] <- WittenSCCA_Perm$v
  recovered_corr [['Witten-CV']] <- cancors_witten
  
  
  Waaijenborg_Delta<-Waaijenborg(X=X,Y=Y,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),rank=rank,selection.criterion=1)
  MSEa[1,6]<-principal_angles(trueA,Waaijenborg_Delta$vhat)$angles[1,1]
  MSEb[1,6]<-principal_angles(trueB,Waaijenborg_Delta$uhat)$angles[1,1]
  cancors_Waaijenborg <- Waaijenborg_Delta$cancors
  TPRa[1,6]<-TPR(trueA[, order(true_corr,decreasing=T)], Waaijenborg_Delta$vhat[,order(cancors_Waaijenborg,decreasing=T)[1:rank]])
  TPRb[1,6]<-TPR(trueB[, order(true_corr,decreasing=T)], Waaijenborg_Delta$uhat[,order(cancors_Waaijenborg,decreasing=T)[1:rank]])
  TNRa[1,6]<-TNR(trueA[, order(true_corr,decreasing=T)], Waaijenborg_Delta$vhat[,order(cancors_Waaijenborg,decreasing=T)[1:rank]])
  TNRb[1,6]<-TNR(trueB[, order(true_corr,decreasing=T)], Waaijenborg_Delta$uhat[,order(cancors_Waaijenborg,decreasing=T)[1:rank]])
  l2loss[1,6] <- sqrt(mean((X_train %*% Waaijenborg_Delta$vhat - Y_train %*% Waaijenborg_Delta$uhat) ^2))
  l2loss_test[1,6] <- sqrt(mean((X_test %*% Waaijenborg_Delta$vhat  - Y_test %*% Waaijenborg_Delta$uhat )^2))
  results.x[['Waaijenborg-Author']] <- Waaijenborg_Delta$vhat
  results.y[['Waaijenborg-Author']] <- Waaijenborg_Delta$uhat
  recovered_corr [['Waaijenborg-Author']] <- cancors_Waaijenborg
  
  Waaijenborg_Test<-Waaijenborg(X=X,Y=Y,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                                rank=rank, selection.criterion=2)
  MSEa[1,7]<-principal_angles(trueA,Waaijenborg_Test$vhat)$angles[1,1]
  MSEb[1,7]<-principal_angles(trueB,Waaijenborg_Test$uhat)$angles[1,1]
  cancors_Waaijenborg <- Waaijenborg_Test$cancors
  TPRa[1,7]<-TPR(trueA[, order(true_corr,decreasing=T)], Waaijenborg_Test$vhat[,order(cancors_Waaijenborg,decreasing=T)])
  TPRb[1,7]<-TPR(trueB[, order(true_corr,decreasing=T)], Waaijenborg_Test$uhat[,order(cancors_Waaijenborg,decreasing=T)])
  TNRa[1,7]<-TNR(trueA[, order(true_corr,decreasing=T)], Waaijenborg_Test$vhat[,order(cancors_Waaijenborg,decreasing=T)])
  TNRb[1,7]<-TNR(trueB[, order(true_corr,decreasing=T)], Waaijenborg_Test$uhat[,order(cancors_Waaijenborg,decreasing=T)])
  l2loss[1,7] <- sqrt(mean((X_train %*% Waaijenborg_Test$vhat - Y_train %*% Waaijenborg_Test$uhat) ^2))
  l2loss_test[1,7] <- sqrt(mean((X_test %*% Waaijenborg_Test$vhat  - Y_test %*% Waaijenborg_Test$uhat )^2))
  results.x[['Waaijenborg-CV']] <- Waaijenborg_Test$vhat
  results.y[['Waaijenborg-CV']] <- Waaijenborg_Test$uhat
  recovered_corr [['Waaijenborg-CV']] <- cancors_Waaijenborg
  
  Parkhomenko_SCCA<-SCCA_Parkhomenko(x.data=X,y.data=Y, Krank=rank)
  MSEa[1,8]<-principal_angles(trueA,Parkhomenko_SCCA$a)$angles[1,1]
  MSEb[1,8]<-principal_angles(trueB,Parkhomenko_SCCA$b)$angles[1,1]
  cancors_Parkhomenko <- Parkhomenko_SCCA$cancor
  TPRa[1,8]<-TPR(trueA[, order(true_corr,decreasing=T)], Parkhomenko_SCCA$a[,order(cancors_Parkhomenko,decreasing=T)])
  TPRb[1,8]<-TPR(trueB[, order(true_corr,decreasing=T)], Parkhomenko_SCCA$b[,order(cancors_Parkhomenko ,decreasing=T)])
  TNRa[1,8]<-TNR(trueA[, order(true_corr,decreasing=T)], Parkhomenko_SCCA$a[,order(cancors_Parkhomenko ,decreasing=T)])
  TNRb[1,8]<-TNR(trueB[, order(true_corr,decreasing=T)], Parkhomenko_SCCA$b[,order(cancors_Parkhomenko ,decreasing=T)])
  l2loss[1,8] <- sqrt(mean((X_train %*% Parkhomenko_SCCA$a - Y_train %*%Parkhomenko_SCCA$b) ^2))
  l2loss_test[1,8] <- sqrt(mean((X_test %*% Parkhomenko_SCCA$a  - Y_test %*% Parkhomenko_SCCA$b)^2))
  results.x[['Parkhomenko-Author']] <- Parkhomenko_SCCA$uhat
  results.y[['Parkhomenko-Author']] <- Parkhomenko_SCCA$vhat
  recovered_corr [['Parkhomenko-Author']] <- cancors_Parkhomenko
  
  RCC_cv<-estim.regul_crossvalidation(X,Y,n.cv=5)
  RCCA<-rcc(X,Y,RCC_cv$lambda1.optim,RCC_cv$lambda2.optim)
  MSEa[1,9]<-principal_angles(trueA,RCCA$xcoef)$angles[1,1]
  MSEb[1,9]<-principal_angles(trueB,RCCA$ycoef)$angles[1,1]
  cancors_RCCA<- RCCA$cor
  TPRa[1,9]<-TPR(trueA[, order(true_corr,decreasing=T)], RCCA$xcoef[,order(cancors_RCCA,decreasing=T)[1:rank]])
  TPRb[1,9]<-TPR(trueB[, order(true_corr,decreasing=T)], RCCA$ycoef[,order(cancors_RCCA,decreasing=T)[1:rank]])
  TNRa[1,9]<-TNR(trueA[, order(true_corr,decreasing=T)], RCCA$xcoef[,order(cancors_RCCA,decreasing=T)[1:rank]])
  TNRb[1,9]<-TNR(trueB[, order(true_corr,decreasing=T)], RCCA$ycoef[,order(cancors_RCCA,decreasing=T)[1:rank]])
  l2loss[1,9] <- sqrt(mean((X_train %*% RCCA$xcoef - Y_train %*%RCCA$ycoef) ^2))
  l2loss_test[1,9] <- sqrt(mean((X_test %*%  RCCA$xcoef - Y_test %*% RCCA$ycoef)^2))
  results.x[['Canonical Ridge-Author']] <- RCCA$xcoef
  results.y[['Canonical Ridge-Author']] <- RCCA$ycoef
  recovered_corr [['Canonical Ridge-Author']] <- cancors_RCCA
  
  #### Maube start by concatenating all the results
  all_results = data.frame(t(rbind(l2loss, l2loss_test,
                      TPRa, TPRb, TNRa, TNRb,
                      MSEa, MSEb)))
  colnames(all_results) = c("l2loss", "l2loss_test",
                            "TPRa", "TPRb", "TNRa", "TNRb",
                            "MSEa", "MSEb")
  all_results["method"] = methods
  
  correl = do.call(rbind, myList <- lapply(seq_along(recovered_corr), 
                                           function(i) {
    temp = data.frame(t(sort(recovered_corr[[i]])[1:rank]))
    colnames(temp) = sapply(1:nrow(temp), function(i){paste0("Correlation", i)})
    temp["method"] = names(recovered_corr)[[i]]
    return(temp)
  }))
  
  all_results = merge(all_results, correl, by="method")
  
  all_coeffs = do.call(rbind, myList <- lapply(seq_along(results.x), function(i) {
    temp = data.frame(results.x[[i]][, 1:k])
    colnames(temp) = sapply(1:ncol(temp), function(i){paste0("C", i)})
    temp["coefficient"] = 1:nrow(temp)
    temp["matrix"]  =  "x"
    temp["method"] = names(results.x)[[i]]
    return(temp)
  }))
  all_coeffs.y = do.call(rbind, myList <- lapply(seq_along(results.y), function(i) {
    temp = data.frame(results.y[[i]][, 1:k])
    colnames(temp) = sapply(1:ncol(temp), function(i){paste0("C", i)})
    temp["coefficient"] = 1:nrow(temp)
    temp["matrix"]  =  "y"
    temp["method"] = names(results.y)[[i]]
    return(temp)
  }))
  all_coeffs = rbind(all_coeffs, all_coeffs.y)
  
  
  all_results["l1"] = NA
  all_results["l2"] = NA
  all_results["l3"] = NA
  all_results["max.iter"] = NA
  all_results["fold"] = NA
  all_results["eta"] = NA
  write_csv(all_results, paste0(mysavedir, name_experiment, 'other_methods_all_results.csv'))
  write_csv(all_coeffs, paste0(mysavedir, name_experiment, 'other_methods_coeffs.csv'))
  
}
if(type_experiment == "SparseCCA"){
  lambda1 = args[4]
  lambda2 = args[5]
  lambda3 = args[6]
  
  
  
}

if(type_experiment == "SparseCCA"){
  lambda1 = args[5]
  lambda2 = args[6]
  lambda3 = args[7]
  max.iter = args[8]
  n.cv = args[9]
  eta = args[10]
  
  id_exp = paste0(type_experiment, lambda1, '-',  lambda2, '-',
                  lambda3, '-')
  
  
  train = sample(1:nrow(X_train), nrow(X_train))
  cv.sets <- split(train, ceiling(seq_along(1:nrow(X_train))/(nrow(X_train)/n.cv)))
  
  
  CCAChao.resultsfull <-  sparseCCA_Chao(X_train, Y_train, rank, lambda1, lambda2,
                                     lambda3,
                                     eta=eta, max_k=30, 
                                     zero_threshold = 1e-6, 
                                     epsilon=0.1,
                                     max.iter = max.iter,
                                     verbose=TRUE)
  
  all_results = c()
  for (i in 1:n.cv.sample){
    validation.sample<- train[c(cv.sets[[i]])]
    training.sample<-train[-c(cv.sets[[i]])]
    Xcv = X_train[training.sample, ]
    Ycv = Y_train[training.sample, ]
    Xval=X_train[validation.sample,]
    Yval=Y_train[validation.sample,]
    CCAChao.results_temp <-  sparseCCA_Chao(Xcv, Ycv, rank, lambda1, lambda2,
                                       lambda3,
                                       eta=eta, max_k=30, 
                                       zero_threshold = 1e-6, 
                                       epsilon=0.1,
                                       max.iter = max.iter,
                                       verbose=TRUE)
    temp_results = data.frame("method" = "ChaoCCA",
                             "l1" = lambda1,
                             "l2" = lambda2,
                             "l3" = lambda3,
                             "eta" = eta,
                             "fold" = i,
                             "max.iter" = max.iter,
                             "l2loss" = sqrt(mean((Xcv %*% CCAChao.results_temp$xcoef - Ycv %*% CCAChao.results_temp$ycoef )^2)),
                             "l2loss_val" =sqrt(mean((Xval %*% CCAChao.results_temp$xcoef  - Yval %*% CCAChao.results_temp$ycoef ) ^2)),
                             "l2loss_test" =sqrt(mean((X_test %*% CCAChao.results$xcoef  - X_test %*% CCAChao.results$ycoef ) ^2)),
                             "TPRa" =TPR(trueA[, order(true_corr,decreasing=T)], CCAChao.results$xcoef[,order(CCAChao.results$cancors,decreasing=T)[1:rank]]), 
                             "TPRb"=TPR(trueB[, order(true_corr,decreasing=T)], CCAChao.results$xcoef[,order(CCAChao.results$cancors,decreasing=T)[1:rank]]), 
                             "TNRa"=TNR(trueA[, order(true_corr,decreasing=T)], CCAChao.results$xcoef[,order(CCAChao.results$cancors,decreasing=T)[1:rank]]), 
                             "TNRb"=TNR(trueB[, order(true_corr,decreasing=T)], CCAChao.results$xcoef[,order(CCAChao.results$cancors,decreasing=T)[1:rank]]),
                             "MSEa"=principal_angles(trueA,CCAChao.results$xcoef)$angles[1,1], 
                             "MSEb"=principal_angles(trueB,CCAChao.results$ycoeff)$angles[1,1]
                             )
    all_results = rbind(all_results, temp_results)
  }
  
  write_csv(all_results, paste0(mysavedir, name_experiment, id_exp,'_all_results.csv'))
  write_csv(all_coeffs, paste0(mysavedir, name_experiment, id_exp, '_coeffs.csv'))
  
}



if(type_experiment == "genCCA"){
  lambda1 = as.numeric(args[5])
  lambda2 = as.numeric(args[6])
  lambda3 = as.numeric(args[7])
  max.iter = as.numeric(args[8])
  
  id_exp = paste0(type_experiment, lambda1, '-',  lambda2, '-',
                  lambda3)
  print(id_exp)
  
  train = sample(1:nrow(X_train), nrow(X_train))
  cv.sets <- split(train, ceiling(seq_along(1:nrow(X_train))/(nrow(X_train)/n.cv)))
  
  genCCA.results <-  genCCA2(X_train, Y_train,
                              Da=D, Db=NULL,
                              lambdaA1=lambda1,
                              lambdaB1=lambda3,
                              lambdaA2=lambda2,
                              lambdaB2=0,
                              rank,
                              A.initial=NULL,B.initial=NULL,
                              max.iter=max.iter,
                              conv=10^-2,
                              solver =  "CGD",
                              standardize=TRUE,
                              verbose=FALSE)
  
  
  all_results = data.frame("method" = "genCCA",
                           "l1" = lambda1,
                           "l2" = lambda2,
                           "l3" = lambda3,
                           "max.iter" = max.iter,
                           "fold"=0,
                           "eta" = NA,
                           "l2loss" = sqrt(mean((X_train %*% genCCA.results$xcoef - Y_train %*% genCCA.results$ycoef ) ^2)),
                           "l2loss_val" = sqrt(mean((X_train %*% genCCA.results$xcoef - Y_train %*% genCCA.results$ycoef ) ^2)),
                           "l2loss_test" =sqrt(mean((X_test %*% genCCA.results$xcoef  - Y_test %*% genCCA.results$ycoef ) ^2)),
                           "TPRa" =TPR(trueA[, order(true_corr,decreasing=T)], genCCA.results$xcoef[,order(genCCA.results$cancors,decreasing=T)[1:rank]]), 
                           "TPRb"=TPR(trueB[, order(true_corr,decreasing=T)], genCCA.results$ycoef[,order(genCCA.results$cancors,decreasing=T)[1:rank]]), 
                           "TNRa"=TNR(trueA[, order(true_corr,decreasing=T)], genCCA.results$xcoef[,order(genCCA.results$cancors,decreasing=T)[1:rank]]), 
                           "TNRb"=TNR(trueB[, order(true_corr,decreasing=T)], genCCA.results$ycoef[,order(genCCA.results$cancors,decreasing=T)[1:rank]]),
                           "MSEa"=principal_angles(trueA, genCCA.results$xcoef)$angles[1,1], 
                           "MSEb"=principal_angles(trueB, genCCA.results$ycoef)$angles[1,1])
  for(k in 1:rank){
    all_results[paste0("Correlation", k)] = genCCA.results$cancors[k]
  }
  
  for (i in 1:n.cv){
    print(i)
    validation.sample<- train[c(cv.sets[[i]])]
    training.sample<-train[-c(cv.sets[[i]])]
    Xcv = X_train[training.sample, ]
    Ycv = Y_train[training.sample, ]
    Xval=X_train[validation.sample,]
    Yval=Y_train[validation.sample,]
    genCCA.results_temp <-  genCCA2(Xcv, Ycv,
                                     Da=D, Db=NULL,
                                     lambdaA1=lambda1,
                                     lambdaB1=lambda3,
                                     lambdaA2=lambda2,
                                     lambdaB2=0,
                                     rank,
                                     A.initial=NULL,B.initial=NULL,
                                     max.iter=max.iter,
                                     conv=10^-2,
                                     solver =  "CGD",
                                     standardize=TRUE,
                                     verbose=FALSE)
    temp_results =   data.frame("method" = "genCCA",
                                              "l1" = lambda1,
                                              "l2" = lambda2,
                                              "l3" = lambda3,
                                              "max.iter" = max.iter,
                                              "fold"=i,
                                              "eta" = NA, 
                                              "l2loss" = sqrt(mean((Xcv %*% genCCA.results_temp$xcoef - Ycv %*% genCCA.results_temp$ycoef ) ^2)),
                                              "l2loss_val" = sqrt(mean((Xval %*% genCCA.results_temp$xcoef - Yval %*% genCCA.results_temp$ycoef ) ^2)),
                                              "l2loss_test" =sqrt(mean((X_test %*% genCCA.results_temp$xcoef  - Y_test %*% genCCA.results_temp$ycoef ) ^2)),
                                              "TPRa" =TPR(trueA[, order(true_corr,decreasing=T)], genCCA.results_temp$xcoef[,order(genCCA.results_temp$cancors,decreasing=T)[1:rank]]), 
                                              "TPRb"=TPR(trueB[, order(true_corr,decreasing=T)], genCCA.results_temp$ycoef[,order(genCCA.results_temp$cancors,decreasing=T)[1:rank]]), 
                                              "TNRa"=TNR(trueA[, order(true_corr,decreasing=T)], genCCA.results_temp$xcoef[,order(genCCA.results_temp$cancors,decreasing=T)[1:rank]]), 
                                              "TNRb"=TNR(trueB[, order(true_corr,decreasing=T)], genCCA.results_temp$ycoef[,order(genCCA.results_temp$cancors,decreasing=T)[1:rank]]),
                                              "MSEa"=principal_angles(trueA, genCCA.results_temp$xcoef)$angles[1,1], 
                                              "MSEb"=principal_angles(trueB, genCCA.results_temp$ycoef)$angles[1,1])
    for(k in 1:rank){
      temp_results[paste0("Correlation", k)] = genCCA.results_temp$cancors[k]
    }
    all_results = rbind(all_results, temp_results)
  }
  
  

  all_coefs.x = data.frame(genCCA.results$ycoef)
  colnames(all_coefs.x) = sapply(1:ncol(all_coefs.x), function(i){paste0("C", i)})
  all_coefs.x["coefficient"] = 1:nrow(all_coefs.y)
  all_coefs.x["matrix"]  =  "y"
  all_coefs.x["method"] = "genCCA"
  all_coefs.x["l1"] = lambda1
  all_coefs.x["l2"] = lambda2
  all_coefs.x["l3"] = lambda3

  
  all_coefs.y = data.frame(genCCA.results$ycoef)
  colnames(all_coefs.y) = sapply(1:ncol(all_coefs.y), function(i){paste0("C", i)})
  all_coefs.y["coefficient"] = 1:nrow(all_coefs.y)
  all_coefs.y["matrix"]  =  "y"
  all_coefs.y["method"] = "genCCA"
  all_coefs.y["l1"] = lambda1
  all_coefs.y["l2"] = lambda2
  all_coefs.y["l3"] = lambda3
  all_coeffs = rbind(all_coefs.x, all_coefs.y)
  
  
  write_csv(all_results, paste0(mysavedir, name_experiment, id_exp,'_all_results.csv'))
  write_csv(all_coeffs, paste0(mysavedir, name_experiment, id_exp, '_coeffs.csv'))
  
}





