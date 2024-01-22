library(dplyr)
library(tidyr)
library(Matrix)
library(glmnet)
library(gglasso)
library(rrpack)
library(foreach)
library(doParallel)
library(CVXR)





CCA_graph_rrr.folds<- function(X, Y, 
                               Gamma, 
                               Sx=NULL, Sy=NULL, kfolds=5,
                               lambda=0.01,
                               r=2, Kx = NULL,
                               do.scale = FALSE, 
                               lambda_Kx=0,
                               LW_Sy = FALSE,
                               rho=10,
                               niter=1e4,
                               thresh=1e-4,
                               Gamma_dagger = NULL) {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- rep(1e8, kfolds)
  p1 <- dim(X)[2]
  p2 <- dim(Y)[2]
  p <- p1 + p2;
  n <- nrow(X)
  
  #init <- gca_to_cca(ainit, S0, pp)
  # loop over folds
  no_cores <- detectCores() - 2  # Save three cores for system processes
  registerDoParallel(cores=no_cores) 
  
  rmse <- foreach(i=seq_along(folds), .combine=c, .packages=c('CVXR', 'Matrix')) %dopar% {
    # ... same code inside your for loop ...
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]
    tryCatch(
      {
        final = CCA_graph_rrr(X_train, Y_train, 
                              Gamma, Sx=NULL,
                              Sy =NULL,
                              lambda = lambda, Kx=Kx, r=r, 
                              lambda_Kx=lambda_Kx,
                              do.scale = do.scale,
                              LW_Sy = LW_Sy,
                              rho=rho,
                              niter=niter,
                              thresh = thresh,
                              Gamma_dagger = Gamma_dagger)
        # make predictions on validation data
        # compute RMSE on validation data
        return(mean((X_val %*% final$U - Y_val%*% final$V)^2))
        #print(rmse)
      },
      error = function(e) {
        #    # If an error occurs, assign NA to the result
        print("An error has occured")
        return(NA)
      })
  }
  
  print(c(lambda, rmse))
  # return mean RMSE across folds
  if (mean(is.na(rmse)) == 1){
    return(1e8)
  }else{
    return(mean(rmse, na.rm=TRUE))
  }
}


CCA_graph_rrr.CV<- function(X, Y, 
                            Gamma,
                            r=2, Kx = NULL, lambda_Kx = 0,
                            param_lambda=10^seq(-3, 1.5, length.out = 10),
                            kfolds=5, 
                            parallelize = FALSE,
                            do.scale=FALSE,
                            LW_Sy = FALSE,
                            rho=10,
                            niter=1e4,
                            thresh=1e-4,
                            Gamma_dagger = NULL
){
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  # if(q >n){
  #   X_temp <- X
  #   X <- Y
  #   Y <- X_temp
  # }
  # if (scale){
  #   X <- scale(X)
  #   Y <- scale(Y)
  # }
  # Sx = t(X) %*% X /n
  # Sy = t(Y) %*% Y /n
  # if (LW_Sy){
  #   lw_cov <- corpcor::cov.shrink(Y)
  #   Sy <- as.matrix(lw_cov)
  # }
  
  if ( n <  min(q,p)){
    print("Warning!!!! Both X and Y are high dimensional, method may fail")
  }
  if (parallelize){
    no_cores <- detectCores() - 5  # Save five cores for system processes
    registerDoParallel(cores=no_cores)
    resultsx <- foreach(lambda=param_lambda, .combine=rbind, .packages=c('CVXR',
                                                                         'Matrix')) %dopar% {
                                                                           rmse <- CCA_graph_rrr.folds(X, Y, Gamma,
                                                                                                       Sx=NULL, Sy=NULL, 
                                                                                                       kfolds=kfolds,
                                                                                                       lambda=lambda, r=r, 
                                                                                                       Kx = Kx, 
                                                                                                       lambda_Kx=lambda_Kx,
                                                                                                       do.scale = do.scale,
                                                                                                       LW_Sy = LW_Sy,
                                                                                                       rho=rho,
                                                                                                       niter=niter,
                                                                                                       thresh=thresh,
                                                                                                       Gamma_dagger= Gamma_dagger)
                                                                           data.frame(lambda=lambda, rmse=rmse)
                                                                         }
  }else{
    print('here')
    resultsx <- expand.grid(lambda = param_lambda) %>%
      mutate(rmse = map_dbl(lambda, ~CCA_graph_rrr.folds(X, Y, Gamma, 
                                                         Sx=NULL, Sy=NULL,
                                                         kfolds=kfolds,
                                                         lambda=.x,
                                                         r=r,
                                                         Kx = Kx, lambda_Kx=lambda_Kx,
                                                         do.scale = do.scale,
                                                         LW_Sy = LW_Sy,
                                                         rho=rho,
                                                         niter=niter,
                                                         thresh=thresh,
                                                         Gamma_dagger = Gamma_dagger)))
  }
  
  
  resultsx$rmse[which(is.na(resultsx$rmse))] = 1e8
  resultsx$rmse[which((resultsx$rmse) ==0)] = 1e8
  resultsx = resultsx %>% filter(rmse > 1e-5)
  #plot(log(resultsx$lambda), resultsx$rmse)
  opt_lambda <- resultsx$lambda[which.min(resultsx$rmse)]
  print(c("selected", opt_lambda))
  rmse = resultsx$rmse
  final <-CCA_graph_rrr(X, Y, Gamma, 
                        Sx = NULL, Sy = NULL, 
                  lambda =opt_lambda, Kx = Kx,
                  r = r,
                  lambda_Kx=lambda_Kx,
                  do.scale = do.scale,
                  LW_Sy = LW_Sy,
                  rho=rho,
                  niter=niter,
                  thresh=thresh,
                  Gamma_dagger = Gamma_dagger)
  
  print(resultsx)
  return(list( ufinal = final$U, 
               vfinal = final$V,
               lambda=opt_lambda,
               resultsx=resultsx,
               rmse = rmse,
               cor = sapply(1:r, function(i){cov(X %*% final$U[,i], Y %*% final$V[,i])})
  ))
  
}


CCA_graph_rrr = function(X, Y, 
                         Gamma, 
                         Sx=NULL, Sy=NULL,
                         Sxy = NULL,
                         lambda =0, Kx, r,
                         do.scale = FALSE, lambda_Kx=0,
                         LW_Sy = FALSE,
                         rho=10,
                         niter=1e4,
                         thresh=1e-4,
                         Gamma_dagger = NULL){
  # solve RRR: ||Y-XB|| + tr(Bt K B)
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  if(q >n){
    X_temp <- X
    X <- Y
    Y <- X_temp
  }
  if (do.scale){
    X <- scale(X)
    Y <- scale(Y)
  }else{
    X <- scale(X, scale=FALSE)
    Y <- scale(Y, scale=FALSE)
  }
  
  if ( n <  min(q,p)){
    print("Warning!!!! Both X and Y are high dimensional, method may fail")
  }
  if (is.null(Sx)){
    Sx = t(X) %*% X /n
  }
  if (is.null(Sy)){
    ###
    Sy = t(Y) %*% Y /n
    if (LW_Sy){
      lw_cov <- corpcor::cov.shrink(Y)
      Sy <- as.matrix(lw_cov)
    }
  }
  
  svd_Sy = svd(Sy)
  sqrt_inv_Sy = svd_Sy$u %*% diag(sapply(svd_Sy$d, function(x){ifelse(x > 1e-4, 
                                                                      1/sqrt(x), 0)}))  %*% t(svd_Sy$v)
  tilde_Y = Y %*% sqrt_inv_Sy
  
  if(is.null(Gamma_dagger)){
    Gamma_dagger = pinv(Gamma)
  }
  Pi = diag(rep(1, p)) - Gamma_dagger %*% Gamma
  XPi = X %*% Pi
  XGamma_dagger = X %*% Gamma_dagger
  #### Start by computing the intercept
  Projection =  pinv(t(XPi) %*% XPi) %*% t(XPi) %*% tilde_Y
  new_Ytilde = tilde_Y - XPi %*% Projection 
  
  
  
  
  if(!is.null(Kx)){
    Sx_tot = Sx + lambda_Kx * as.matrix(Kx) 
  }else{
    Sx_tot = Sx
  }
  

  
  print("Using CVXR")
  print("lambda is")
  print(lambda)
  ### Use CVXR
  new_p = ncol(XGamma_dagger)
  
  #### Let's use ADMM
  Sx_tot_newp = t(XGamma_dagger) %*% XGamma_dagger/n
  U = matrix(0, new_p, q)
  Z = matrix(0, new_p, q)
  prod_xy = t(XGamma_dagger) %*% new_Ytilde / n
  invSx = solve(Sx_tot_newp + rho *diag(rep(1, new_p)))
  for (i in 1:niter){
    Uold = U
    Zold = Z
    B = invSx %*% (prod_xy  + rho* (Z - U))
    Bold = B
    Z = B + U
    norm_col = sapply(1:nrow(Z), function(i){sqrt(sum(Z[i,]^2))})
    index_0 = which(norm_col < lambda/rho)
    index_pos = which(norm_col  > lambda/rho)
    if(length(index_0)>0){
      Z[index_0,] = 0
    }
    print(length(index_pos))
    if(length(index_pos)>1){
      Z[index_pos,] = diag(sapply(index_pos, function(x){ 1- (lambda/rho)/norm_col[x]})) %*% Z[index_pos,] 
    }else{
      if(length(index_pos) == 1){
        Z[index_pos,] = ( 1- (lambda/rho)/norm_col[index_pos]) * Z[index_pos,] 
      }
    }
    U = U + B - Z
    print(c("ADMM iter", i, norm(Z - B), norm(Zold - Z), norm(Uold - U)))
    if (max(c(norm(Z - B)/sqrt(new_p), norm(Zold - Z)/sqrt(new_p))) <thresh){
      break
    }
  }
  
  B_opt = Gamma_dagger %*% B + Pi %*% Projection
  B_opt[which(abs(B_opt)<1e-5)] = 0
  print(B_opt)
    
  svd_Sx = svd(Sx_tot)
  sqrt_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, sqrt(x), 0)}))  %*% t(svd_Sx$u)
  sqrt_inv_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sx$u)
  sol = svd(sqrt_Sx %*% B_opt)
  print(sqrt_Sx %*% B_opt)
  #sol = svd(t(B_OLS[I, ]), nu = r, nv=r)
  V = sqrt_inv_Sy %*% sol$v[, 1:r]
  #U = matrix(0, p, r)
  # B = U \tilde{V} 
  inv_D = diag(sapply(1:r, FUN=function(x){ifelse(sol$d[x]<1e-4, 0, 1/sol$d[x])}))
  U = B_opt %*% sol$v[, 1:r] %*% inv_D ### = U\lambda
  print(t(U) %*% Sx %*% U)
  print(t(V) %*% Sy %*% V)
  
  
  
  loss = mean((Y %*% V - X %*% U)^2)
  return(list(U = U, V = V, loss = loss,
              cor = sapply(1:r, function(i){cov(X %*% U[,i], Y %*% V[,i])})))
}
