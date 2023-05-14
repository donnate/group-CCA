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
library(igraph)
#library(corpcor)
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

PROBA = list('11'= 0.08, '12'=0.005, '13'=0.005, 
            '22'= 0.07, '23' = 0.001,
            '33' = 0.02)



generate_localized_example <- function(n, p1,  nnzeros = 5,
                             theta = diag( c(0.9,  0.8)),
                             a = 0, r=2, signal_strength="normal",
                             type_graph="SBM", power=1.8, probs=PROBA,
                             threshold_limit = 0.8){
  #n  <- 200;
  #p1 <- 100;
  #p2 <- 100;
  # r <- 2
  p2  <- p1
  p <- p1 + p2;
  pp <- c(p1,p2);
  print('--------------------------------------');
  print('Generating data ...');
  #a <- 0.3;
  generated_graph = FALSE
  trials = 1

  while(generated_graph == FALSE){
      if (type_graph == "PA"){
        G <- sample_pa(p1, power = power)
        #G <- make_lattice(length = 10, dim = 3)
      }else{
        if (type_graph == "SBM"){
          pref.matrix = rbind(c(probs$`11`, probs$`12`, probs$`13`),
                              c(probs$`12`, probs$`22`, probs$`23`),
                              c(probs$`13`, probs$`23`, probs$`33`))
          G.sbm <- sample_sbm(p1, pref.matrix = pref.matrix,
                              block.sizes = c(floor(p1/3), floor(p1/3), p1- 2*floor(p1/3)))
          comp = components(G.sbm)
          G <- delete_vertices(G.sbm, which(comp$membership!=1))
          Z.u = c(rep(1, floor(p1/3)), rep(2, floor(p1/3)), rep(3,  p1 - 2*floor(p1/3)))
          Z.u = Z.u[which(comp$membership==1)]
        }else{
          if(as.integer(sqrt(p1))^2 == p1){
            G <- make_lattice(length = as.integer(sqrt(p1)), 
                              dim = as.integer(sqrt(p1)) )
          }else{
            print("p1 must have an integer square root")
          }
        }
      }
      
      components.G = components(G)
      E = data.frame(as_edgelist(G))
      colnames(E)  = c("x", "y")
      E["e"] = 1:nrow(E)
      E = pivot_longer(E, cols=-c("e"))
      E["fill_value"] = sapply(E$name,  vfn <- function(x){
        ifelse(x=="x", 1, -1)
      })
      D = tidyr::pivot_wider(E, id_cols =c("e"), names_from = "value", values_from  =  "fill_value")
      D[is.na(D)] <- 0
      D = as.matrix(D %>% dplyr::select(-e))
      daggerDx = pinv(D)
      daggerD = bdiag(daggerDx,daggerDx)
      
      A = as_adjacency_matrix(G,sparse=FALSE)
      cov = diag(rep(1, nrow(A))) 
      #inv_deg = diag(apply(A,1, function(x){ ifelse(sum(x) >0, sqrt(1/sum(x)), 0)}))
      L = laplacian_matrix(G, sparse=FALSE, normalized = TRUE)  #### Use this as the covariance matrix
      sole = which(diag(L) == 0)
      alpha=2.1
      L = L[-c(sole), -c(sole)] + alpha * diag(rep(1, (length(V(G)) -length(sole))))
      #inv_deg = diag(apply(A,1, function(x){  sqrt(1/(alpha +sum(x)))}))
      #L = inv_deg %*% (L + alpha * diag(rep(1, length(V(G))))) %*% inv_deg
      L = (L + t(L))/2
      #mean(eigen(L)$values >0)
      #### invert
      cov[-c(sole), -c(sole)] =  solve(L) 
      cov[which(abs(cov)<1e-3)] = 0
      cov = diag(1/sqrt(diag(cov))) %*% cov %*% diag(1/sqrt(diag(cov))) 
      print("cov done")
      generated_graph = min(eigen(cov)$values) >0
      if(trials > 10){
        print("Error: unable to generate appropriate covariance")
        return()
      }
      trials = trials + 1
  }

  # generate covariance matrix for X and Y
  #T1 = toeplitz(a^(0:(pp[1]-1)));
  Sigma <- diag(p1+p2)
  print(c(pp[1], pp[2],r))
  
  # generate covariance matrix for X and Y
  u = matrix(0, pp[1], r)
  v = matrix(0, pp[2], r)
  print(c(p1, nnzeros))
  s  = sample(1:p1, nnzeros)

  T1 = cov;
  T2 = cov;
  Tss = T1;
  Sigma[1:p1, 1:p1] = T1;
  Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)] = T1;

  u[s,1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  v[s,1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  for (i in 1:r){
    #old = u[,i]
    u[,i] = mvrnorm(n=1, mu=u[,i], Sigma = Tss)
    u[which(abs(u[,i])<threshold_limit),i]=0
    v[,i] = mvrnorm(n=1, mu=v[,i], Sigma = Tss)
    v[which(abs(v[,i])<threshold_limit),i]=0
    #plot(u[,i], old)
  }
  s_selected.x = which(apply(u, 1, norm)>0)
  s_selected.y = which(apply(v, 1, norm)>0)
  u[s_selected.x,] <- u[s_selected.x,] %*%(sqrtm(t(u[s_selected.x,1:r]) %*% cov[s_selected.x, s_selected.x] %*% u[s_selected.x,1:r])$Binv)
  v[s_selected.y,] <- v[s_selected.y,] %*%(sqrtm(t(v[s_selected.y,1:r]) %*% cov[s_selected.y, s_selected.y] %*% v[s_selected.y,1:r])$Binv)
  u[which(abs(u) <1e-5)] = 0
  v[which(abs(v) <1e-5)] = 0

  
  Sigma[(p1+1):(p1+p2), 1:p1] = T2 %*%  v  %*% theta %*% t(u) %*% T1;
  Sigma[1:p1, (p1+1):(p1+p2)] = t(Sigma[(p1+1):(p1+p2), 1:p1])
  Sigma[which(abs(Sigma)<1e-3)] = 0
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
  Mask_dual = bdiag(matrix(1, length(E(G)),length(E(G))),
                     matrix(1, length(E(G)),length(E(G))))
  Sigma0 = Sigma * Mask;
  
  S <- cov(Data)
  sigma0hat <- S * Mask
  
  # Estimate the subspace spanned by the largest eigenvector using convex relaxation and TGD
  # First calculate ground truth
  result = geigen::geigen(Sigma, Sigma0)
  evalues <- result$ values
  evectors <-result$vectors
  evectors <- evectors[,p:1]
  a <- evectors[,1:r]
  #scale <- a %*% sqrtm(diag(r)+t(a) %*% Sigma %*% a/lambda)$B;
  return(list(Sigma=Sigma, Sigma0=Sigma0, daggerD = daggerD,
              daggerDx =daggerDx,
              daggerDy = daggerDx,
              components.G = components.G,
         S = S, sigma0hat =  sigma0hat, Mask= Mask,
         Mask_dual=Mask_dual,
         X=X, Y = Y, Data=Data,u=u, v=v, 
         Sigmax=Sigmax, Sigmay=Sigmay,
         a=a
        ))
}


pipeline_group_lasso <- function(Data, Mask, sigma0hat, Gamma, r, nu=1, Sigmax, 
                                    Sigmay, maxiter=30, lambdax=NULL, lambday=NULL,
                                    adaptive=TRUE, kfolds=5, param1=10^(seq(-4, 2, by = 0.25)),
                                    create_folds=TRUE, init ="Fantope", normalize=FALSE){
  
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
  }else{
    S1 <- cov(Data)
    S3 <- cov(Data)
    X = Data[,1:p1]
    Y = Data[,(p2+1):p]
    X1 = Data[,1:p1]  %*% daggerD1
    Y1 = Data[,(p2+1):p]  %*% daggerD2
  }
  S1 = daggerD %*% S1 %*%  t(daggerD) 
  sigma0hat1 = daggerD %*% (S1 * Mask) %*%  t(daggerD) 
  if (init == "Fantope"){
      ag <- sgca_init(A=S1, B=sigma0hat1, rho=0.5 * sqrt(log(dim(sigma0hat1)[1])/dim(X)[1]),
                  K=r ,nu=nu,trace=FALSE, maxiter = maxiter) ###needs to be changed to be a little more scalable
      ainit <- init_process(ag$Pi, r) 
  
  }else{
      RCC_cv<-estim.regul_crossvalidation(X1, 
                                          Y1,
                                          n.cv=5)
      method<-rcc(X1, Y1, RCC_cv$lambda1.optim, RCC_cv$lambda2.optim)
      ainit= rbind(method$xcoef[,1:r], method$ycoef[,1:r])
  }

  init <- gca_to_cca(ainit, S3, pp)
  print("Init done")
  initu<- init$u
  initv <- init$v


  if (is.null(lambdax)){
    ### do CV
    resultsx <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_gamma_sparse_function(X, Y, kfolds, initu, initv,
                                                  lambdax = .x, adaptive=adaptive,
                                                   normalize=normalize)))

    
    # print best hyperparameters and corresponding RMSE
    best_hyperparams <- resultsx[which.min(resultsx$rmse), ]
    which_lambdax = which(abs(resultsx$rmse-min(resultsx$rmse))/(1e-6 + min(resultsx$rmse)) <0.05)
    lambdax = max(resultsx$param1[which_lambdax])
  }
  if (is.null(lambday)){
    ### do CV
    results <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_gamma_sparse_function(Y, X, kfolds, initv, initu,
                                                  lambdax = .x, adaptive=adaptive)))
    best_hyperparams <- results[which.min(results$rmse), ]
    which_lambday = which(abs(results$rmse-min(results$rmse))/(1e-6 + min(results$rmse)) <0.05)
    lambday = max(results$param1[which_lambday])
  }
  
  ufinal = gamma_sparse_lasso(X, Y %*% initv, Gamma,
                              adaptive=adaptive, lambdax, 
                              max.iter=5000, 
                              max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
  
  vfinal = gamma_sparse_lasso(Y, X %*% initu, Gamma,
                       adaptive=adaptive, lambday, max.iter=5000, 
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