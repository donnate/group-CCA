#!/usr/bin/env Rscript
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source('experiments/alternative_methods/regular_CCA.R')
source("src/cca_bis.R")
source("src/CCA_alaChao.R")
source("src/GroupChao.R")

source('src/metrics.R')

library(tidyverse)
library(PMA)
args = commandArgs(trailingOnly=TRUE)

#### load(data)
name_experiment  = args[1]
type_experiment = args[2]
mysavedir = 'experiments/results/'
mydatadir = 'experiments/data/'
lambda1 = as.numeric(args[3])
lambda2 = as.numeric(args[4])
lambda3 = as.numeric(args[5])
max.iter = as.numeric(args[6])
penalty_type_chao = args[7]



load(paste0(mydatadir, name_experiment, '-environment.RData'))
print("yolo")
print(type_experiment)


Sigma_x = t(X_train)%*% X_train
Sigma_y = t(Y_train)%*% Y_train
Sigma_xy = t(X_train)%*% Y_train

if (type_experiment == "others"){
  
  reg.CCA<-regular_cca(X=X_train,Y=Y_train,rank=rank)
  result = evaluate_method(reg.CCA$xcoef,
                            reg.CCA$ycoef, 
                            trueA, trueB,
                            Sigma_x, Sigma_y,
                            X_train, Y_train,
                            X_test, Y_test,
                            D,
                            normalize=FALSE)
  result = data.frame(result)
  result["method"] = "regular-CCA"

  
  FIT_SAR_BIC<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                         lambdaAseq=seq(from=0.2,to=0.02,length=10),
                         lambdaBseq=seq(from=0.2,to=0.02,length=10),
                         max.iter=100,conv=10^-2,
                         selection.criterion=1,n.cv=5)
  result_temp = evaluate_method(FIT_SAR_BIC$ALPHA,
                            FIT_SAR_BIC$BETA,
                            trueA, trueB,
                            Sigma_x, Sigma_y,
                            X_train, Y_train,
                            X_test, Y_test,
                            D,
                            normalize=TRUE)
  result_temp = data.frame(result_temp)
  result_temp["method"] = "SAR-Author"
  result = rbind(result, result_temp)

  
  
  FIT_SAR_CV<-SparseCCA(X=X,Y=Y,rank=rank,
                        lambdaAseq=seq(from=0.2,to=0.02,length=10),
                        lambdaBseq=seq(from=0.2,to=0.02,length=10),
                        max.iter=100,conv=10^-2, selection.criterion=2, n.cv=5)
  result_temp = data.frame(evaluate_method(FIT_SAR_CV$ALPHA,
                            FIT_SAR_CV$BETA, 
                            trueA, trueB,
                            Sigma_x, Sigma_y,
                            X_train, Y_train,
                            X_test, Y_test,
                            D,
                            normalize=TRUE))
  result_temp["method"] = "SAR-CV"
  result = rbind(result, result_temp)
  
  
  Witten_Perm<-CCA.permute(x=X,z=Y,typex="standard",typez="standard",nperms=50)
  WittenSCCA_Perm<-CCA(x=X,z=Y,typex="standard",typez="standard",K=rank,penaltyx=Witten_Perm$bestpenaltyx,penaltyz=Witten_Perm$bestpenaltyz,trace=F)
  result_temp = data.frame(evaluate_method(WittenSCCA_Perm$u,
                            WittenSCCA_Perm$v, 
                            trueA, trueB,
                            Sigma_x, Sigma_y,
                            X_train, Y_train,
                            X_test, Y_test,
                            D,
                            normalize=TRUE))
  result_temp["method"] = "Witten-Author"
  result = rbind(result, result_temp)
  
  
  Witten_CV<-Witten.CV(X=X,Y=Y,n.cv=5,lambdax=matrix(seq(from=0,to=1,length=20),nrow=1),lambday=matrix(seq(from=0,to=1,length=20),nrow=1))
  WittenSCCA_CV<-CCA(x=X,z=Y,typex="standard",typez="standard",K=rank,penaltyx=Witten_CV$lambdax.opt,penaltyz=Witten_CV$lambday.opt,trace=F)
  result_temp = data.frame(evaluate_method(WittenSCCA_Perm$u,
                            WittenSCCA_Perm$v,
                            trueA, trueB,
                            Sigma_x, Sigma_y,
                            X_train, Y_train,
                            X_test, Y_test,
                            D,
                            normalize=TRUE))
  result_temp["method"] = "Witten-CV"
  result = rbind(result, result_temp)

  
  
  Waaijenborg_Delta<-Waaijenborg(X=X,Y=Y,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),rank=rank,selection.criterion=1)
  result_temp = data.frame(evaluate_method(Waaijenborg_Delta$vhat,
                            Waaijenborg_Delta$uhat, 
                            trueA, trueB,
                            Sigma_x, Sigma_y,
                            X_train, Y_train,
                            X_test, Y_test,
                            D,
                            normalize=TRUE))
  result_temp["method"] = "Waaijenborg-Author"
  result = rbind(result, result_temp)

  
  Waaijenborg_Test<-Waaijenborg(X=X,Y=Y,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                                rank=rank, selection.criterion=2)
  result_temp = data.frame(evaluate_method(Waaijenborg_Test$vhat,
                            Waaijenborg_Test$uhat, 
                            trueA, trueB,
                            Sigma_x, Sigma_y,
                            X_train, Y_train,
                            X_test, Y_test,
                            D,
                            normalize=TRUE))
  result_temp["method"] = "Waaijenborg-CV"
  result = rbind(result, result_temp)
 
  
  Parkhomenko_SCCA<-SCCA_Parkhomenko(x.data=X,y.data=Y, Krank=rank)
  result_temp = data.frame(evaluate_method(Parkhomenko_SCCA$a,
                            Parkhomenko_SCCA$b, 
                            trueA, trueB,
                            Sigma_x, Sigma_y,
                            X_train, Y_train,
                            X_test, Y_test,
                            D,
                            normalize=TRUE))
  result_temp["method"] = "Parkhomenko-Author"
  result = rbind(result, result_temp)

  
  RCC_cv<-estim.regul_crossvalidation(X,Y,n.cv=5)
  RCCA<-rcc(X,Y,RCC_cv$lambda1.optim,RCC_cv$lambda2.optim)
  result_temp = data.frame(evaluate_method(as.matrix(RCCA$xcoef[,1:rank], ncol=1),
                            as.matrix(RCCA$ycoef[,1:rank], ncol=1),
                            trueA, trueB,
                            Sigma_x, Sigma_y,
                            X_train, Y_train,
                            X_test, Y_test,
                            D,
                            normalize=TRUE))
  result_temp["method"] = "Canonical Ridge-Author"
  result = rbind(result, result_temp)
  
  result["l1"] = NA
  result["l2"] = NA
  result["l3"] = NA
  result["max.iter"] = NA
  result["fold"] = NA
  result["eta"] = NA
  write_csv(result, paste0(mysavedir, name_experiment, '-others-all-results.csv'))  
}

if(type_experiment == "original-ChaoCCA"){
  eta = 2
  id_exp = paste0(name_experiment, '-', type_experiment, '-', lambda1, '-',  lambda2, '-',
                  lambda3, '-')
  
  
  train = sample(1:nrow(X_train), nrow(X_train))
  cv.sets <- split(train, ceiling(seq_along(1:nrow(X_train))/(nrow(X_train)/n.cv)))
  
  
   
  CCAChao.resultsfull <- sparseCCA_Chao(X_train, Y_train, rank, lambda1, lambda2,
                                     lambda3,
                                     eta=eta, max_k=30, 
                                     zero_threshold = 1e-6, 
                                     epsilon=0.01,
                                     max.iter = max.iter,
                                     verbose=TRUE)
  
  result = data.frame(evaluate_method(CCAChao.resultsfull$xcoef,
                              CCAChao.resultsfull$ycoef, 
                              trueA, trueB,
                              Sigma_x, Sigma_y,
                              X_train, Y_train,
                              X_test, Y_test,
                              D,
                              normalize=FALSE))
  result["method"] =  "CCA-Chao-orginal"
  result["fold"] = 0

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
                                     epsilon=0.01,
                                     max.iter = max.iter,
                                     verbose=TRUE)
    temp_result = evaluate_method(CCAChao.results_temp$xcoef,
                              CCAChao.results_temp$ycoef, 
                              trueA, trueB,
                              Sigma_x, Sigma_y,
                              Xcv, Ycv,
                              Xval, Yval,
                              D,
                              normalize=TRUE)
    temp_result = data.frame(temp_result)
    temp_result["method"] = "CCA-Chao-orginal"
    temp_result["fold"] = i
    result = rbind(result, temp_result)
  }
  result["l1"] = lambda1
  result["l2"] = lambda2
  result["l3"] = lambda3
  result["max.iter"] = max.iter
  result["eta"] = eta
  write_csv(result, paste0(mysavedir, id_exp,'-all-results.csv'))
  
}

if(type_experiment == "genChaoCCA"){
  print("yolo2")
  eta = 2
  id_exp = paste0(name_experiment, '-', type_experiment, '-', penalty_type_chao, '-', lambda1, '-',  lambda2, '-',
                  lambda3, '-')
  
  
  train = sample(1:nrow(X_train), nrow(X_train))
  cv.sets <- split(train, ceiling(seq_along(1:nrow(X_train))/(nrow(X_train)/n.cv)))
  
  
  CCAChao.resultsfull <-  genCCA_Chao(D, X_train, Y_train, rank,
                                      lambda1, lambda2, lambda3,  
                                      penalty_type = penalty_type_chao,
                                      eta=eta, 
                                      max_k=30, zero_threshold = 1e-6, 
                                      epsilon=0.1,
                                      max.iter = 5000,
                                      verbose=TRUE)
                      #  sparseCCA_Chao(X_train, Y_train, rank, lambda1, lambda2,
                      #                lambda3,
                      #                eta=eta, max_k=30, 
                      #                zero_threshold = 1e-6, 
                      #                epsilon=0.1,
                      #                max.iter = max.iter,
                      #                verbose=TRUE)
  
  result = data.frame(evaluate_method(CCAChao.resultsfull$xcoef,
                              CCAChao.resultsfull$ycoef, 
                              trueA, trueB,
                              Sigma_x, Sigma_y,
                              X_train, Y_train,
                              X_test, Y_test,
                              D,
                              normalize=FALSE))
  result["method"] = paste0("CCA-Chao-", penalty_type_chao)
  result["fold"] = 0
  write_csv(result, paste0(mysavedir, id_exp,'-all-results.csv'))

  for (i in 1:n.cv.sample){
    validation.sample<- train[c(cv.sets[[i]])]
    training.sample<-train[-c(cv.sets[[i]])]
    Xcv = X_train[training.sample, ]
    Ycv = Y_train[training.sample, ]
    Xval=X_train[validation.sample,]
    Yval=Y_train[validation.sample,]
    CCAChao.results_temp <-  genCCA_Chao(D, Xcv, Ycv, rank,
                                      lambda1, lambda2, lambda3,  
                                      penalty_type = penalty_type_chao,
                                      eta=eta, 
                                      max_k=30, zero_threshold = 1e-5, 
                                      epsilon=0.01,
                                      max.iter = 5000,
                                      verbose=TRUE)
    temp_result = evaluate_method(CCAChao.results_temp$xcoef,
                              CCAChao.results_temp$ycoef, 
                              trueA, trueB,
                              Sigma_x, Sigma_y,
                              Xcv, Ycv,
                              Xval, Yval,
                              D,
                              normalize=TRUE)
    temp_result = data.frame(temp_result)
    temp_result["method"] = paste0("CCA-Chao-", penalty_type_chao)
    temp_result["fold"] = i
    result = rbind(result, temp_result)
    write_csv(result, paste0(mysavedir, id_exp,'-all-results.csv'))
  }
  result["l1"] = lambda1
  result["l2"] = lambda2
  result["l3"] = lambda3
  result["max.iter"] = max.iter
  result["eta"] = eta
  write_csv(result, paste0(mysavedir, id_exp,'-all-results.csv'))
  
}



if(type_experiment == "genCCA"){
  
  id_exp = paste0(name_experiment, '-', type_experiment, '-',  lambda2, '-',
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
                              conv=conv,
                              solver =  "CGD",
                              standardize=TRUE,
                              verbose=TRUE)
  
  print("Solved once")

  print(dim(genCCA.results$xcoef))
  result = data.frame(evaluate_method(genCCA.results$xcoef,
                              genCCA.results$ycoef,
                              trueA, trueB, 
                              Sigma_x, Sigma_y,
                              X_train, Y_train,
                              X_test, Y_test,
                              D,
                              normalize=TRUE))
    
  result["method"] = "genCCA"
  result["fold"] = 0
  write_csv(result, paste0(mysavedir, id_exp,'-all-results.csv'))
  
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
                                     conv=conv,
                                     solver =  "CGD",
                                     standardize=TRUE,
                                     verbose=FALSE)
    result_temp = data.frame(evaluate_method(genCCA.results_temp$xcoef,
                                  genCCA.results_temp$ycoef, 
                                  trueA, trueB,
                                  Sigma_x, Sigma_y,
                                  Xcv, Ycv,
                                  Xval, Yval,
                                  D,
                                  normalize=TRUE))
    result_temp["method"] = "genCCA"
    result_temp["fold"] = i
    result= rbind(result, result_temp)
    write_csv(result, paste0(mysavedir, id_exp,'-all-results.csv'))
  }

  result["l1"] = lambda1
  result["l2"] = lambda2
  result["l3"] = lambda3
  result["max.iter"] = max.iter
  result["eta"] = eta


  write_csv(result, paste0(mysavedir, id_exp,'-all-results.csv'))
  
}





