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
require(gtools)
numCores <- detectCores()

vfn <- function(x){
  ifelse(x=="x", 1, -1)
}


sort_cancors <- function(U, V, rank=3){
  all_perms = gtools::permutations(n = rank, r = rank, v = 1:rank)
  result = apply(all_perms, 1, function(x){
    #print(x)
    sqrt(mean((U - V[, x]) ^2))
    #print((U[1,] - V[1, x]) ^2)
  })
  return(list("loss"=result[which.min(result)], "perm"=all_perms[which.min(result),],
              "all_results" = result))
}
library(foreach)    # install.packages('foreach')
library(caret)      # install.packages('caret', dependencies = c("Depends", "Suggests"))
library(doParallel) # install.packages('doParallel')
registerDoParallel(makeCluster(4)) # Use 4 cores for parallel CV




simple_experiment <- function(n, p, q, sigma, k, effect_size = 2,
                              sigma_noise=0.1, power=1,
                              probs = list('11'= 0.8, '12'=0.05, '13'=0.05, 
                                           '22'= 0.7, '23' = 0.01,
                                           '33' = 0.2),
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
  if (type_graph == "pa"){
    G <- sample_pa(p, power = power)
  }else{
   if (type_graph == "SBM"){
    pref.matrix = rbind(c(probs$`11`, probs$`12`, probs$`13`),
                        c(probs$`12`, probs$`22`, probs$`23`),
                        c(probs$`13`, probs$`23`, probs$`33`))
    G <- sample_sbm(p, pref.matrix = pref.matrix,
                    block.sizes = c(floor(p/3), floor(p/3), p - 2*floor(p/3)))
    Z = c(rep(1, floor(p/3)), rep(2, floor(p/3)), rep(3,  p - 2*floor(p/3)))
    }
  }
  E = data.frame(as_edgelist(G))
  colnames(E)  = c("x", "y")
  E["e"] = 1:nrow(E)
  E = pivot_longer(E, cols=-c("e"))
  E["fill_value"] = sapply(E$name, vfn)
  D = pivot_wider(E, id_cols =c("e"), names_from = "value", values_from  =  fill_value)
  D[is.na(D)] <- 0
  D = as.matrix(D[, 2:(p+1)])
  
  X = matrix(rnorm(p * n, mean=0, sd=sigma_noise), n, p)
  Y = matrix(rnorm(q*n,  mean=0, sd=sigma_noise), n, q)
  source = matrix(0, k, p)
  colors = matrix(0, k, p)
  colors_l = rep(0, p)
  trueA = matrix(0, p, k)
  trueB = matrix(0, q, k)
  true_corr = rep(0, k)
  
  if ( type_graph  %in% c("pa")){
    indices <- sample(1:p, 1)
    #### Make sure the selected clusters are independent
    not_all_indices = TRUE
    while(not_all_indices){
      print("here")
      for (j in 2:k){
        found_vertex = FALSE
        iter = 1
        while(found_vertex==FALSE){
          iter  = iter + 1
          index <- sample(1:p, 1)
          d = sapply(indices, FUN=function(x){distances( G, v =x, to=index)})
          print(d)
          if(min(d)> 2 * egonet_size + 1){
            indices <- c(indices, index)
            found_vertex = TRUE
          }
          if(iter > 100){
            break;
          }
        }
      }
      not_all_indices = (length(indices) < k)
    }
    print(indices)
    for (i in 1:k){
      print(i)
      idx = indices[i]
      subg <- ego(G, order=egonet_size, nodes = idx, 
                  mode = "all", mindist = 0)[[1]]
      nodes_in_network <- as.numeric(subg)
      colors_l[nodes_in_network] = i
      colors[i, nodes_in_network]  = i
      source[i, nodes_in_network] = 1
      #mean_value <- sapply(1:n, function(i){rnorm(1, mean=effect_size, sd=sigma)})
      #X[, nodes_in_network] <- mean_value + X[, nodes_in_network] 
      Y[, i] =Y[, i] + X[, nodes_in_network] %*% matrix(rep(1, length(nodes_in_network)), nrow=length(nodes_in_network), ncol=1)
      true_corr[i] = cor(Y[, i], X[, nodes_in_network] %*% matrix(rep(1, length(nodes_in_network)), nrow=length(nodes_in_network), ncol=1))
      trueA[nodes_in_network, i] =  1/sqrt(sum(X[, source[i,]]^2))
      trueB[i,i] = 1/ sqrt(sum(Y[,i]^2))
    }
  }else{
    for (i in 1:k){
      print(i)
      nodes_in_network = which(Z==i)
      colors_l[nodes_in_network] = i
      colors[i, nodes_in_network]  = i
      source[i, nodes_in_network] = 1
      #mean_value <- sapply(1:n, function(i){rnorm(1, mean=effect_size, sd=sigma)})
      X[, nodes_in_network] <- delta + X[, nodes_in_network] 
      Y[, i] =Y[, i] + X[, nodes_in_network] %*% matrix(rep(2, length(nodes_in_network)), nrow=length(nodes_in_network), ncol=1)
      true_corr[i] = cor(Y[, i], X[, nodes_in_network] %*% matrix(rep(1, length(nodes_in_network)), nrow=length(nodes_in_network), ncol=1))
      trueA[nodes_in_network, i] =  1/sqrt(sum(X[, source[i,]]^2))
      trueB[i,i] = 1/ sqrt(sum(Y[,i]^2))
    }
    
  }
  
  
  
  
  
  
  
  if(plot){
    source_df = data.frame(source)
    source_df["component"] = 1:k
    ggplot(pivot_longer(source_df, cols=-c("component"))) +
      geom_tile(aes(x=name,y=component, fill=as.factor(value)))
    #### we can also simply plot the graph
    g <- simplify(G)
    V(g)$color= "black"
    V(g)$color[which(source[1,]>0)] =  "lightblue"
    V(g)$color[which(source[2,]>0)] =   "orange"
    V(g)$color[which(source[3,]>0)] =  "red"
    plot(g)
    
  }
  df.x <- data.frame(X) %>% mutate_all(~(scale(.) %>% as.vector))
  df.y <- data.frame(Y) %>% mutate_all(~(scale(.) %>% as.vector))
  X <- as.matrix(df.x)
  Y <- as.matrix(df.y)
  
  
  # Store results Estimation Accuracy
  Nsim = 1
  results.x <-list()
  results.y <-list()
  recovered_corr <- list()
  
  #### divide data into train and test
  train  = sample(1:n, size=as.integer(0.8 *n), replace = FALSE)
  test = setdiff(1:n, train)
  X_train = X[train, ]
  X_test = X[test, ]
  Y_train = Y[train, ]
  Y_test = Y[test, ]
  
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

  
  # if (p<n){
  Sigma_x = t(X_train)%*% X_train
  Sigma_y = t(Y_train)%*% Y_train
  Sigma_xy = t(X_train)%*% Y_train
  svd_x = svd(Sigma_x)
  inv_sqrt_Sigma_x = svd_x$u %*% diag(sapply(svd_x$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_x$v)
  svd_y = svd(Sigma_y)
  inv_sqrt_Sigma_y = svd_y$u %*% diag(sapply(svd_y$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_y$v)
  svd_for_cca = svd(inv_sqrt_Sigma_x %*% Sigma_xy %*% inv_sqrt_Sigma_y)
  xcoef = inv_sqrt_Sigma_x %*% svd_for_cca$u[, 1:rank]
  ycoef = inv_sqrt_Sigma_y %*% svd_for_cca$v[, 1:rank]
  cancor <- svd_for_cca$d
  
  MSEa[1,1]<-principal_angles(trueA, xcoef)$angles[1,1]
  MSEb[1,1]<-principal_angles(trueB, ycoef[, 1:rank])$angles[1,1]
  TPRa[1,1]<-TPR(trueA,xcoef)
  TPRb[1,1]<-TPR(trueB,ycoef)
  TNRa[1,1]<-TNR(trueA,xcoef)
  TNRb[1,1]<-TNR(trueB,ycoef)
  l2loss[1,1] <- sqrt(mean((X_train %*% xcoef- Y_train %*% ycoef)^2))
  l2loss_test[1, 1] <- sqrt(mean((X_test %*% xcoef- Y_test %*% ycoef)^2))
  results.x[['Regular-CCA']] <- xcoef 
  results.y[['Regular-CCA']] <- ycoef
  recovered_corr [['Regular-CCA']] <- cancor
  
  res_df = data.frame(xcoef)
  res_df["coef"] = 1:nrow(res_df)
  ggplot(pivot_longer(res_df, cols=-c("coef"))) +
    geom_tile(aes(x=name,y=coef, fill=value))
  

  # 
  save.image(file=paste0(mysavedir, namefile, 'environment.RData'))
  
  
  FIT_SAR_BIC<-SparseCCA(X=X_train,Y=Y_train,rank=k,
                         lambdaAseq=seq(from=20,to=0.02,length=20),
                         lambdaBseq=seq(from=20,to=0.02,length=20),
                         max.iter=100,conv=10^-2,
                         selection.criterion=1,n.cv=5)
  
  res = sort_cancors(X_test %*% FIT_SAR_BIC$ALPHA, Y_test %*% FIT_SAR_BIC$BETA, rank)
  ###### 
  svd_sol = svd(t(FIT_SAR_BIC$ALPHA) %*% Sigma_x %*%FIT_SAR_BIC$ALPHA  )
  inv_sqrt_sol = svd_sol$u %*% diag(sapply(svd_sol$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol$v)
  svd_sol2 = svd(t(FIT_SAR_BIC$BETA) %*% Sigma_y %*%FIT_SAR_BIC$BETA  )
  inv_sqrt_sol2 = svd_sol2$u %*% diag(sapply(svd_sol2$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol2$v)
  results.x[['SAR-Author']] <- FIT_SAR_BIC$ALPHA %*% inv_sqrt_sol 
  results.y[['SAR-Author']] <- FIT_SAR_BIC$BETA %*% inv_sqrt_sol2
  
  MSEa[1,2]<-principal_angles(trueA,results.x[['SAR-Author']])$angles[1,1]
  MSEb[1,2]<-principal_angles(trueB,results.y[['SAR-Author']])$angles[1,1]
  cancors_SAR_BIC <- FIT_SAR_BIC$cancors
  TPRa[1,2]<-TPR(trueA,FIT_SAR_BIC$ALPHA[,1:k])
  TPRb[1,2]<-TPR(trueB,FIT_SAR_BIC$BETA[,1:k])
  TNRa[1,2]<-TNR(trueA,FIT_SAR_BIC$ALPHA[,1:k])
  TNRb[1,2]<-TNR(trueB,FIT_SAR_BIC$BETA[,1:k])
  l2loss[1,2] <- sqrt(mean((X_train %*% results.x[['SAR-Author']] - Y_train %*% results.y[['SAR-Author']])^2))
  l2loss_test[1,2] <- sqrt(mean((X_test %*% results.x[['SAR-Author']] - Y_test %*% results.y[['SAR-Author']])^2))
  recovered_corr [['SAR-Author']] <- cancors_SAR_BIC
  
  
  FIT_SAR_CV<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                        lambdaAseq=seq(from=0.2,to=0.02,length=10),
                        lambdaBseq=seq(from=0.2,to=0.02,length=10),
                        max.iter=100,conv=10^-2, selection.criterion=2, n.cv=5)
  res = sort_cancors(X_test %*%   FIT_SAR_CV$ALPHA, Y_test %*% FIT_SAR_CV$BETA, rank)
  ###### 
  svd_sol = svd(t(FIT_SAR_CV$ALPHA) %*% Sigma_x %*%FIT_SAR_CV$ALPHA  )
  inv_sqrt_sol = svd_sol$u %*% diag(sapply(svd_sol$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol$v)
  svd_sol2 = svd(t(FIT_SAR_CV$BETA) %*% Sigma_y %*%FIT_SAR_CV$BETA  )
  inv_sqrt_sol2 = svd_sol2$u %*% diag(sapply(svd_sol2$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol2$v)
  results.x[['SAR-CV']] <- FIT_SAR_CV$ALPHA %*% inv_sqrt_sol 
  results.y[['SAR-CV']] <- FIT_SAR_CV$BETA %*% inv_sqrt_sol2
  recovered_corr [['SAR-CV']] <- cancors_SAR_CV
  
  MSEa[1,3]<-principal_angles(trueA, FIT_SAR_CV$ALPHA)$angles[1,1]
  MSEb[1,3]<-principal_angles(trueB, FIT_SAR_CV$BETA)$angles[1,1]
  cancors_SAR_CV <- FIT_SAR_CV$cancors
  TPRa[1,3]<-TPR(trueA,FIT_SAR_CV$ALPHA[,order(cancors_SAR_CV,decreasing=T)[1:k]])
  TPRb[1,3]<-TPR(trueB,FIT_SAR_CV$BETA[,order(cancors_SAR_CV,decreasing=T)[1:k]])
  TNRa[1,3]<-TNR(trueA,FIT_SAR_CV$ALPHA[,order(cancors_SAR_CV,decreasing=T)[1:k]])
  TNRb[1,3]<-TNR(trueB,FIT_SAR_CV$BETA[,order(cancors_SAR_CV,decreasing=T)[1:k]])
  l2loss[1,3] <- sqrt(mean((X_train %*% results.x[['SAR-CV']] - Y_train %*% results.y[['SAR-CV']])^2))
  l2loss_test[1, 3] <- sqrt(mean((X_test %*% results.x[['SAR-CV']] - Y_test %*% results.y[['SAR-CV']])^2))

  
  
  Witten_Perm<-CCA.permute(x=X_train,z=Y_train,typex="standard",typez="standard", nperms=50)
  WittenSCCA_Perm<-CCA(x=X_train, z=Y_train, typex="standard",typez="standard",K=rank,penaltyx=Witten_Perm$bestpenaltyx,penaltyz=Witten_Perm$bestpenaltyz,trace=F)
  svd_sol = svd(t(WittenSCCA_Perm$u) %*% Sigma_x %*%WittenSCCA_Perm$u  )
  inv_sqrt_sol = svd_sol$u %*% diag(sapply(svd_sol$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol$v)
  svd_sol2 = svd(t(WittenSCCA_Perm$v) %*% Sigma_y %*%WittenSCCA_Perm$v )
  inv_sqrt_sol2 = svd_sol2$u %*% diag(sapply(svd_sol2$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol2$v)
  results.x[['Witten-Author']] <- WittenSCCA_Perm$u%*% inv_sqrt_sol
  results.y[['Witten-Author']] <- WittenSCCA_Perm$v%*% inv_sqrt_sol2
  recovered_corr [['Witten-Author']] <- cancors_witten
  
  MSEa[1,4]<-principal_angles(trueA,WittenSCCA_Perm$u)$angles[1,1]
  MSEb[1,4]<-principal_angles(trueB,WittenSCCA_Perm$v)$angles[1,1]
  cancors_witten <- WittenSCCA_Perm$cors
  TPRa[1,4]<-TPR(trueA, WittenSCCA_Perm$u)
  TPRb[1,4]<-TPR(trueB, WittenSCCA_Perm$v)
  TNRa[1,4]<-TNR(trueA, WittenSCCA_Perm$u)
  TNRb[1,4]<-TNR(trueB, WittenSCCA_Perm$v)
  l2loss[1,4] <- sqrt(mean((X_train %*% results.x[['Witten-Author']] - Y_train %*% results.y[['Witten-Author']])^2))
  l2loss_test[1,4] <- sqrt(mean((X_test %*% results.x[['Witten-Author']] - Y_test %*% results.y[['Witten-Author']])^2))
  
  Witten_CV<-Witten.CV(X=X_train,Y=Y_train,n.cv=5,lambdax=matrix(seq(from=0,to=1,length=20),nrow=1),lambday=matrix(seq(from=0,to=1,length=20),nrow=1))
  WittenSCCA_CV<-CCA(x=X_train,z=Y_train,typex="standard",typez="standard",K=rank,penaltyx=Witten_CV$lambdax.opt,penaltyz=Witten_CV$lambday.opt,trace=F)
  svd_sol = svd(t(WittenSCCA_CV$u) %*% Sigma_x %*%WittenSCCA_CV$u  )
  inv_sqrt_sol = svd_sol$u %*% diag(sapply(svd_sol$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol$v)
  svd_sol2 = svd(t(WittenSCCA_CV$v) %*% Sigma_y %*%WittenSCCA_CV$v )
  inv_sqrt_sol2 = svd_sol2$u %*% diag(sapply(svd_sol2$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol2$v)
  results.x[['Witten-CV']] <-WittenSCCA_CV$u %*% inv_sqrt_sol
  results.y[['Witten-CV']] <- WittenSCCA_CV$v%*% inv_sqrt_sol2
  recovered_corr [['Witten-CV']] <- cancors_witten
  
  MSEa[1,5]<-principal_angles(trueA,WittenSCCA_CV$u)$angles[1,1]
  MSEb[1,5]<-principal_angles(trueB,WittenSCCA_CV$v)$angles[1,1]
  cancors_witten <- WittenSCCA_CV$cors
  TPRa[1,5]<-TPR(trueA,WittenSCCA_CV$u)
  TPRb[1,5]<-TPR(trueB,WittenSCCA_CV$v)
  TNRa[1,5]<-TNR(trueA,WittenSCCA_CV$u)
  TNRb[1,5]<-TNR(trueB,WittenSCCA_CV$v)
  l2loss[1,5] <- sqrt(mean((X_train %*% results.x[['Witten-CV']] - Y_train %*% results.y[['Witten-CV']])^2))
  l2loss_test[1,5] <- sqrt(mean((X_test %*% results.x[['Witten-CV']] - Y_test %*% results.y[['Witten-CV']])^2))
  
  Waaijenborg_Delta<-Waaijenborg(X=X_train,Y=Y_train,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                                 lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),rank=rank,selection.criterion=1)
  svd_sol = svd(t(Waaijenborg_Delta$vhat) %*% Sigma_x %*%Waaijenborg_Delta$vhat  )
  inv_sqrt_sol = svd_sol$u %*% diag(sapply(svd_sol$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol$v)
  svd_sol2 = svd(t(Waaijenborg_Delta$uhat) %*% Sigma_y %*%Waaijenborg_Delta$uhat )
  inv_sqrt_sol2 = svd_sol2$u %*% diag(sapply(svd_sol2$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol2$v)
  results.x[['Waaijenborg-Author']] <- Waaijenborg_Delta$vhat%*% inv_sqrt_sol
  results.y[['Waaijenborg-Author']] <- Waaijenborg_Delta$uhat%*% inv_sqrt_sol2
  recovered_corr [['Waaijenborg-Author']] <- cancors_Waaijenborg
  
  MSEa[1,6]<-principal_angles(trueA,Waaijenborg_Delta$vhat)$angles[1,1]
  MSEb[1,6]<-principal_angles(trueB,Waaijenborg_Delta$uhat)$angles[1,1]
  cancors_Waaijenborg <- Waaijenborg_Delta$cancors
  TPRa[1,6]<-TPR(trueA, Waaijenborg_Delta$vhat)
  TPRb[1,6]<-TPR(trueB, Waaijenborg_Delta$uhat)
  TNRa[1,6]<-TNR(trueA, Waaijenborg_Delta$vhat)
  TNRb[1,6]<-TNR(trueB, Waaijenborg_Delta$uhat)
  l2loss[1,6] <- sqrt(mean((X_train %*% results.x[['Waaijenborg-Author']]- Y_train %*% results.y[['Waaijenborg-Author']]) ^2))
  l2loss_test[1,6] <- sqrt(mean((X_test %*% results.x[['Waaijenborg-Author']]  - Y_test %*% results.y[['Waaijenborg-Author']] )^2))

  
  Waaijenborg_Test<-Waaijenborg(X=X_train,
                                Y=Y_train,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                                lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),
                                rank=k, selection.criterion=2)
  svd_sol = svd(t(Waaijenborg_Test$vhat) %*% Sigma_x %*%Waaijenborg_Test$vhat  )
  inv_sqrt_sol = svd_sol$u %*% diag(sapply(svd_sol$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol$v)
  svd_sol2 = svd(t(Waaijenborg_Test$uhat) %*% Sigma_y %*%Waaijenborg_Test$uhat )
  inv_sqrt_sol2 = svd_sol2$u %*% diag(sapply(svd_sol2$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol2$v)
  results.x[['Waaijenborg-CV']] <- Waaijenborg_Test$vhat %*% inv_sqrt_sol
  results.y[['Waaijenborg-CV']] <- Waaijenborg_Test$uhat %*% inv_sqrt_sol2
  recovered_corr [['Waaijenborg-CV']] <- cancors_Waaijenborg
  MSEa[1,7]<-principal_angles(trueA,Waaijenborg_Test$vhat)$angles[1,1]
  MSEb[1,7]<-principal_angles(trueB,Waaijenborg_Test$uhat)$angles[1,1]
  cancors_Waaijenborg <- Waaijenborg_Test$cancors
  TPRa[1,7]<-TPR(trueA, Waaijenborg_Test$vhat)
  TPRb[1,7]<-TPR(trueB, Waaijenborg_Test$uhat)
  TNRa[1,7]<-TNR(trueA, Waaijenborg_Test$vhat)
  TNRb[1,7]<-TNR(trueB, Waaijenborg_Test$uhat)
  l2loss[1,7] <- sqrt(mean((X_train %*% results.x[['Waaijenborg-CV']]- Y_train %*% results.y[['Waaijenborg-CV']]) ^2))
  l2loss_test[1,7] <- sqrt(mean((X_test %*% results.x[['Waaijenborg-CV']]  - Y_test %*% results.y[['Waaijenborg-CV']] )^2))

  
  Parkhomenko_SCCA<-SCCA_Parkhomenko(x.data=X_train,
                                     y.data=Y_train, Krank=k)
  svd_sol = svd(t(Parkhomenko_SCCA$a) %*% Sigma_x %*%Parkhomenko_SCCA$a )
  inv_sqrt_sol = svd_sol$u %*% diag(sapply(svd_sol$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol$v)
  svd_sol2 = svd(t(Parkhomenko_SCCA$b) %*% Sigma_y %*%Parkhomenko_SCCA$b )
  inv_sqrt_sol2 = svd_sol2$u %*% diag(sapply(svd_sol2$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol2$v)
  results.x[['Parkhomenko-Author']] <- Parkhomenko_SCCA$a %*% inv_sqrt_sol
  results.y[['Parkhomenko-Author']] <- Parkhomenko_SCCA$b %*% inv_sqrt_sol2
  recovered_corr [['Parkhomenko-Author']] <- cancors_Parkhomenko
  MSEa[1,8]<-principal_angles(trueA,Parkhomenko_SCCA$a)$angles[1,1]
  MSEb[1,8]<-principal_angles(trueB,Parkhomenko_SCCA$b)$angles[1,1]
  cancors_Parkhomenko <- Parkhomenko_SCCA$cancor
  TPRa[1,8]<-TPR(trueA, Parkhomenko_SCCA$a)
  TPRb[1,8]<-TPR(trueB, Parkhomenko_SCCA$b)
  TNRa[1,8]<-TNR(trueA, Parkhomenko_SCCA$a)
  TNRb[1,8]<-TNR(trueB, Parkhomenko_SCCA$b)
  l2loss[1,8] <- sqrt(mean((X_train %*% results.x[['Parkhomenko-Author']]- Y_train %*%results.y[['Parkhomenko-Author']]) ^2))
  l2loss_test[1,8] <- sqrt(mean((X_test %*% results.x[['Parkhomenko-Author']] - Y_test %*%results.y[['Parkhomenko-Author']])^2))

  
  RCC_cv<-estim.regul_crossvalidation(X_train,Y_train,n.cv=5)
  RCCA<-rcc(X,Y,RCC_cv$lambda1.optim,RCC_cv$lambda2.optim)
  svd_sol = svd(t(RCCA$xcoef[, 1:k]) %*% Sigma_x %*%RCCA$xcoef[, 1:k] )
  inv_sqrt_sol = svd_sol$u %*% diag(sapply(svd_sol$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol$v)
  svd_sol2 = svd(t(RCCA$ycoef[, 1:k]) %*% Sigma_y %*%RCCA$ycoef[, 1:k] )
  inv_sqrt_sol2 = svd_sol2$u %*% diag(sapply(svd_sol2$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol2$v)
  results.x[['Canonical Ridge-Author']] <- RCCA$xcoef[, 1:k] %*% inv_sqrt_sol
  results.y[['Canonical Ridge-Author']] <- RCCA$ycoef[, 1:k] %*% inv_sqrt_sol2
  recovered_corr [['Canonical Ridge-Author']] <- cancors_RCCA
  MSEa[1,9]<-principal_angles(trueA,RCCA$xcoef[, 1:k])$angles[1,1]
  MSEb[1,9]<-principal_angles(trueB,RCCA$ycoef[, 1:k])$angles[1,1]
  cancors_RCCA<- RCCA$cor
  TPRa[1,9]<-TPR(trueA, RCCA$xcoef[, 1:k])
  TPRb[1,9]<-TPR(trueB, RCCA$ycoef[, 1:k])
  TNRa[1,9]<-TNR(trueA, RCCA$xcoef[, 1:k])
  TNRb[1,9]<-TNR(trueB, RCCA$ycoef[, 1:k])
  l2loss[1,9] <- sqrt(mean((X_train %*% results.x[['Canonical Ridge-Author']] - Y_train %*%results.y[['Canonical Ridge-Author']]) ^2))
  l2loss_test[1,9] <- sqrt(mean((X_test %*%  results.x[['Canonical Ridge-Author']]- Y_test %*%results.y[['Canonical Ridge-Author']])^2))

  
  ### Should probably time this
  Fit.genCCA <- genCCA2(X=X_train,Y=Y_train, Da = D, Db=NULL, 
                        A.initial =RCCA$xcoef[, 1:k],
                        B.initial =RCCA$ycoef[, 1:k],
                        rank=rank, 
                        lambdaA1=0.1, 
                        lambdaA2=0.1,
                        lambdaB1 =0,
                        lambdaB2 =0,
                        max.iter=max.iter, 
                        conv=conv, 
                        solver= solver)
  
  genCCA.results <-  genCCA.CV(X_train, Y_train, D, k, n.cv=2,
                            lambda1seq=seq(from=0.1, 10, by=10),
                            lambda2seq=seq(from=0.1, 10, by=10),
                            lambdaB1seq=c(0),
                            lambdaB2seq=c(0))
  svd_sol = svd(t(Fit.genCCA$xcoef[, 1:k]) %*% Sigma_x %*%Fit.genCCA$xcoef[, 1:k] )
  inv_sqrt_sol = svd_sol$u %*% diag(sapply(svd_sol$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol$v)
  svd_sol2 = svd(t(Fit.genCCA$ycoef[, 1:k]) %*% Sigma_y %*%Fit.genCCA$ycoef[, 1:k] )
  inv_sqrt_sol2 = svd_sol2$u %*% diag(sapply(svd_sol2$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol2$v)
  results.x[['genCCA']] <- Fit.genCCA$xcoef %*% inv_sqrt_sol
  results.y[['genCCA']] <- Fit.genCCA$ycoef %*% inv_sqrt_sol2
  #recovered_corr [['genCCA']] <- genCCA.results$cancors
  MSEa[1,10]<-principal_angles(trueA,Fit.genCCA$xcoef)$angles[1,1]
  MSEb[1,10]<-principal_angles(trueB,Fit.genCCA$ycoef)$angles[1,1]
  #cancors = diag(t(X %*% genCCA.results$fit.gcc$xcoef) %*% (Y %*% genCCA.results$fit.gcc$ycoef))
  TPRa[1,10]<-TPR(trueA, Fit.genCCA$xcoef)
  TPRb[1,10]<-TPR(trueB, Fit.genCCA$ycoef)
  TNRa[1,10]<-TNR(trueA, Fit.genCCA$xcoef)
  TNRb[1,10]<-TNR(trueB, Fit.genCCA$ycoef)
  l2loss[1,10] <- sqrt(mean((X_train %*%  results.x[['genCCA']] - Y_train %*% results.y[['genCCA']]) ^2))
  l2loss_test[1,10] <- sqrt(mean((X_test %*%   results.x[['genCCA']] - Y_test %*%  results.y[['genCCA']])^2))
  
  
  ###
  if(plot){
    source_df = data.frame(source)
    source_df["component"] = 1:k
    ggplot(pivot_longer(source_df, cols=-c("component"))) +
      geom_tile(aes(x=name,y=component, fill=as.factor(value)))
    #### we can also simply plot the graph
    g <- simplify(G)
    V(g)$color <- scales::dscale(as.numeric(genCCA.results$fit.gcc$xcoef[,1]) %>% cut(10), diverging_pal)
    plot(g)
    V(g)$color <- scales::dscale(as.numeric(genCCA.results$fit.gcc$xcoef[,2]) %>% cut(10), diverging_pal)
    plot(g)
    V(g)$color <- scales::dscale(as.numeric(genCCA.results$fit.gcc$xcoef[,3]) %>% cut(10), diverging_pal)
    plot(g)
    
    
    res_df = data.frame(genCCA.results$fit.gcc$xcoef)
    res_df["coef"] = 1:nrow(res_df)
    ggplot(pivot_longer(res_df, cols=-c("coef"))) +
      geom_tile(aes(x=name,y=coef, fill=value))
    
    plot(genCCA.results$fit.gcc$xcoef[,3], trueA[,1])
    plot(genCCA.results$fit.gcc$xcoef[,1], trueA[,2])
    plot(genCCA.results$fit.gcc$xcoef[,2], trueA[,3])
    
    res_df = data.frame(cor(trueA, genCCA.results$fit.gcc$xcoef))
    res_df["component"] = 1:k
    ggplot(pivot_longer(res_df, cols=-c("component"))) +
      geom_tile(aes(x=name,y=component, fill=abs(value)))
  }
  
  
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
  
  
  

  

  
  ### Select the best out of all the different options
  ### selected lambdas = which.min(cv.results)
  #### look at the predicted correlations
  
  #### Assess the quality of the reconstruction
  
  
  
  
  return(list(X=X, Y=Y, Z=Z, 
              gencca_results = gencca_results.final,
              cc_results=cc_results, cv.results= cv.results ))
}


plot_results_on_graph <- function(xcoef, i, G){
  g <- simplify(G)
  #V(g)$color <- scales::dscale(as.numeric(xcoef[,1]) %>% cut(10), diverging_pal)
  V(g)$color <- scales::dscale(as.numeric(xcoef[,i]) %>% cut(10), diverging_pal)
  plot(g)
}

plot_results<- function(xcoef, trueA){
  res_df = data.frame(cor(trueA, xcoef))
  res_df["component"] = 1:k
  ggplot(pivot_longer(res_df, cols=-c("component"))) +
    geom_tile(aes(x=name,y=component, fill=abs(value)))
}


plot_results<- function(xcoef, trueA){
  res_df = data.frame(xcoef)
  res_df["type"] = "estimated"
  res_df["component"] = 1:ncol(xcoef)
  res_df2 = data.frame(trueA)
  res_df2["type"] = "GT"
  res_df2["component"] = 1:ncol(xcoef)
  res_df = rbind( res_def, 
                  res_def2)
  ggplot(pivot_longer(res_df, cols=-c("component"))) +
    geom_tile(aes(x=name,y=component, fill=abs(value)))
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
