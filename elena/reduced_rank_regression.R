library(dplyr)
library(tidyr)
library(Matrix)

CCA_rrr = function(X, Y, Sx=NULL, Sy=NULL,
                  lambda =0, Kx, r, highdim=FALSE, 
                  penalty = "l21", lambda_Kx=0, solver="rrr",
                  LW_Sy = FALSE){
  # solve RRR: ||Y-XB|| + tr(Bt K B)
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  if(q >n){
    X_temp <- X
    X <- Y
    Y <- X_temp
  }
  X <- scale(X, scale = FALSE)
  Y <- scale(Y, scale = FALSE)
  
  if ( n <  min(q,p)){
    print("Warning!!!! Both X and Y are high dimensional, method may fail")
  }
  if (is.null(Sx)){
    Sx = t(X) %*% X /n
  }
  if (is.null(Sy)){
    Sy = t(Y) %*% Y /n
    if (LW_Sy){
      lw_cov <- corpcor::cov.shrink(Y)
      Sy <- as.matrix(lw_cov)
    }
  }
  
  #Sy = t(Y) %*% Y /n
  svd_Sy = svd(Sy)
  sqrt_inv_Sy = svd_Sy$u %*% diag(sapply(svd_Sy$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sy$u)
  tilde_Y = Y %*% sqrt_inv_Sy
  
  if(!is.null(Kx)){
    Sx_tot = Sx + lambda_Kx * as.matrix(Kx) 
  }else{
    Sx_tot = Sx
  }
  Sxy = t(X) %*% tilde_Y/ n 
  if(!highdim){
    B_OLS = solve(Sx_tot) %*% Sxy
    svd_Sx = svd(Sx)
    sqrt_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, sqrt(x), 0)}))  %*% t(svd_Sx$u)
    sqrt_inv_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sx$u)
    sol = svd(sqrt_Sx %*% B_OLS)
    V = sqrt_inv_Sy %*% sol$v[, 1:r]
    U = sqrt_inv_Sx %*% sol$u[, 1:r]
    #sol = svd(t(B_OLS), nu = r, nv=r)
    #V = sqrt_inv_Sy %*% sol$u[, 1:r]
    #U_temp = B_OLS %*%  sol$u[, 1:r]
    #U = sqrt_inv_Sx %*% svd(sqrt_Sx %*% U_temp)$u
    #U = sqrt_inv_Sx %*% U_temp
    print(t(U) %*% Sx %*% U)
    print(t(V) %*% Sy %*% V)
  } else {
    if (penalty == "l21"){
      if (solver =="CVXR"){
        print("Using CVSR")
        ### Use CVXR
        B <- Variable(p, q)
        objective <- Minimize(1/n * sum_squares(tilde_Y - X %*% B) + lambda * sum(norm2(B, axis=1)))
        problem <- Problem(objective)
        result <- solve(problem)
        B_opt <- result$getValue(B)
      }else{
        print("Using Solver")
        test <- cv.srrr( tilde_Y, X, nrank = r,
                         method ="glasso",
                         nfold = 2, norder = NULL,
                        A0 = NULL,   V0 = NULL,
                                          modstr = list("lamA" = rep(lambda, 10),
                                                        "nlam" = 10))
        B_opt <- test$coef # test$U  %*% test$D %*% t(test$V)
      }
      B_opt[which(abs(B_opt)<1e-5)] = 0
      I = which(apply(B_opt^2, 1, sum) > 0)
      if (length(I) >0){
        svd_Sx = svd(Sx[I, I])
        sqrt_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, sqrt(x), 0)}))  %*% t(svd_Sx$u)
        sqrt_inv_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sx$u)
        sol = svd(sqrt_Sx %*% B_opt[I,])
        #sol = svd(t(B_OLS[I, ]), nu = r, nv=r)
        V = sqrt_inv_Sy %*% sol$v[, 1:r]
        U = matrix(0, p, r)
        U[I,] = sqrt_inv_Sx %*% sol$u[, 1:r]
      }else{
        U = matrix(0, p, r)
        V = sqrt_inv_Sy %*% test$V[, 1:r]
      }
      
      #svd_Sx = svd(Sx[I, I])
      #sqrt_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, sqrt(x), 0)}))  %*% t(svd_Sx$u)
      #sqrt_inv_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sx$u)
      #sol = svd(sqrt_Sx %*% B_opt[I,])
      #sol = svd(t(B_OLS[I, ]), nu = r, nv=r)
      #V = sqrt_inv_Sy %*% test$V[, 1:r]
      #U = matrix(0, p, r)
      #U[I,] = sqrt_inv_Sx %*% test$U[1:r]
      #U_temp = B_OLS[I, ] %*%  sol$u[, 1:r]
      
      #U[I, ] = sqrt_inv_Sx %*% svd(sqrt_Sx %*% U_temp)$u
      print(t(U) %*% Sx %*% U)
      print(t(V) %*% Sy %*% V)
    }

  }
  loss = mean((Y %*% V - X %*% U)^2)
  return(list(U = U, V = V, loss = loss))
}




CCA_rrr.CV<- function(X, Y, 
                     r=2, Kx = NULL, lambda_Kx = 0,
                     param_lambda=10^seq(-3, 1.5, length.out = 100),
                     kfolds=10,penalty="l21", solver="rrr",
                     LW_Sy = FALSE
                     ){
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  if(q >n){
    X_temp <- X
    X <- Y
    Y <- X_temp
  }
  X <- scale(X, scale = FALSE)
  Y <- scale(Y, scale = FALSE)
  Sx = t(X) %*% X /n
  Sy = t(Y) %*% Y /n
  if (LW_Sy){
    lw_cov <- corpcor::cov.shrink(Y)
    Sy <- as.matrix(lw_cov)
  }
  
  if ( n <  min(q,p)){
    print("Warning!!!! Both X and Y are high dimensional, method may fail")
  }
  if (solver=="CVRX"){
    resultsx <- expand.grid(lambda = param_lambda) %>%
      mutate(rmse = map_dbl(lambda, ~CCA_rrr.folds(X, Y, Sx=Sx, Sy=Sy,
                                                   kfolds=kfolds, init,
                                                   lambda=.x,
                                                   r=r, penalty=penalty,
                                                   Kx = Kx, lambda_Kx =lambda_Kx)))
    resultsx = resultsx %>% filter(rmse > 1e-5) 
    opt_lambda <- resultsx$lambda[which.min(resultsx$rmse)]
    print(c("selected", opt_lambda))
    
    final <-CCA_rrr(X, Y, lambda =opt_lambda, Kx,
                    r, highdim=TRUE, 
                    penalty = penalty, lambda_Kx=lambda_Kx)
    
    print(resultsx)
  }else{
    svd_Sy = svd(Sy)
    sqrt_inv_Sy = svd_Sy$u %*% diag(sapply(svd_Sy$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sy$u)
    tilde_Y = Y %*% sqrt_inv_Sy
    print("Using rrr solver")
    test <- cv.srrr(tilde_Y, X, nrank = r,
                     method ="glasso",
                     nfold = kfolds,
                     norder = NULL,
                     A0 = NULL,
                     V0 = NULL,
                     modstr = list("lamA" = param_lambda, 
                                   "nlam" = length(param_lambda)))
    opt_lambda <- test$lambda[test$minid]
    resultsx = data.frame("lambda" = test$lambda,
                          "rmse" = test$cv.path[1:length(param_lambda)])
    final <- list(U = test$U,
                  V = sqrt_inv_Sy %*% test$V)
  }
  return(list( ufinal = final$U, 
               vfinal = final$V,
               lambda=opt_lambda,
               resultsx=resultsx
  ))
  
}


CCA_rrr.folds<- function(X, Y, Sx, Sy, kfolds=5, init,
                  lambda=0.01,
                  r=2, penalty="l21", Kx = NULL,
                  lambda_Kx = 0) {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- rep(1e8, kfolds)
  p1 <- dim(X)[2]
  p2 <- dim(Y)[2]
  p <- p1 + p2;
  n <- nrow(X)
  
  #init <- gca_to_cca(ainit, S0, pp)
  # loop over folds
  for (i in seq_along(folds)) {
    #print("here")
    # split data into training and validation sets
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]
    
    # fit model on training data with hyperparameters
    tryCatch(
      {
        final = CCA_rrr(X_train, Y_train, Sx=Sx,
                        Sy =Sy,
                        lambda = lambda, Kx=Kx, r, 
                    highdim=TRUE, penalty = "l21", lambda_Kx=lambda_Kx)
        # make predictions on validation data
        # compute RMSE on validation data
        rmse[i] <- mean((X_val %*% final$U - Y_val%*% final$V)^2)
        #print(rmse)
      },
      error = function(e) {
        #    # If an error occurs, assign NA to the result
        print("An error has occured")
        rmse[i] <- NA
      })
    #print(rmse)
    #}
  }
  
  # return mean RMSE across folds
  if (mean(is.na(rmse)) == 1){
    return(1e8)
  }else{
    return(mean(rmse, na.rm=TRUE))
  }
}



CCA_step = function(X, Y, Kx, Ky, k, highdim = FALSE){
  tr = transform(Y, Ky)
  Y = scale(Y, center = tr$mu, scale = F)
  Ytilde = Y %*% tr$tilde
  rrr = RRR(X, Ytilde, Kx, k, highdim)
  #for imputation
  Yhat = rrr$Yhat %*% tr$notilde 
  Yhat = scale(Yhat, center = -tr$mu, scale = F)
  #cca solution
  U = rrr$U
  V = tr$tilde %*% rrr$V
  return(list(U = U, V = V, Yhat = Yhat, rrrloss = rrr$loss))
}

CCAimpute = function(X, Y, Kx = NULL, Ky = NULL, k = min(nrow(X), min(ncol(X), ncol(Y))), init = NULL, eps = 1e-6, maxiter = 1000, verbose = F){
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  if(p > n & is.null(Kx)) stop("Regularization of X coefficients is requred when p > n")
  if(q > n & is.null(Ky)) stop("Regularization of Y coefficients is requred when q > n")
  
  missY = is.na(Y)
  missX = is.na(X)
  if(sum(missY) > 0){
    if(sum(missX) > 0) rrrside = "both"   
    else rrrside = "Y" 
  } else {
    if(p > q) rrrside = "Y"
    else rrrside = "X"
  }
  
  #initialize
  iter = 0
  delta = Inf
  loss = Inf
  info = c()
  if(is.null(init)){
    X = data.frame(X) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))) %>% as.matrix()
    Y = data.frame(Y) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))) %>% as.matrix()
  } else {
    X = init$X
    Y = init$Y
  }
  
  while(delta > eps & iter < maxiter){
    iter = iter + 1
    loss0 = loss
    
    #CCA step (impute Y)
    if(rrrside %in% c("Y", "both")){
      if(verbose) cat("RRR with response Y\n")
      if(p > n){
        step = CCA_step(X, Y, Kx, Ky, k, highdim = T)
        if(verbose) cat("High dim on X\n")
      } 
      else step = CCA_step(X, Y, Kx, Ky, k)
      Y[missY] = step$Yhat[missY]
      imploss = mean((Y - step$Yhat)^2)
      U = step$U
      V = step$V
    }
    
    #CCA step (impute X)
    if(rrrside %in% c("X", "both")){
      if(verbose) cat("RRR with response X\n")
      if(q > n){
        step = CCA_step(Y, X, Ky, Kx, k, highdim = T)
        if(verbose) cat("High dim on Y\n")
      } 
      else step = CCA_step(Y, X, Ky, Kx, k)
      X[missX] = step$Yhat[missX]
      imploss = mean((X - step$Yhat)^2)
      U = step$V
      V = step$U
    }
    
    XU = X %*% U
    YV = Y %*% V
    loss = diag(cor(XU, YV))[1]
    info = rbind(info, c(iter, loss, step$rrrloss, imploss, delta))
    if(iter > 1) delta = abs((loss0 - loss)/loss0)
    if(verbose) cat("iter ", iter, "cor", loss, "imp.loss", imploss, "delta", delta, "\n\n")
  }
  colnames(info) = c("iter", "cor", "rrrloss", "imploss", "delta")
  sx = diag(t(XU) %*% XU)/(n-1)
  if(!is.null(Kx)) sx = sx + diag(t(U) %*% Kx %*% U)
  sy = diag(t(YV) %*% YV)/(n-1)
  if(!is.null(Ky)) sy = sy + diag(t(V) %*% Ky %*% V)
  
  return(list(U = scale(U, center = F, scale = sqrt(sx)), V = scale(V, center = F, scale = sqrt(sy)), X = X, Y = Y, info = data.frame(info)))
}

CCAimpute_l1= function(X, Y, lambdax = NULL, lambday = NULL, k = k, eps = list(outer = 1e-6, inner = 1e-6), maxiter = list(outer = 1000, inner = 1000), verbose = list(outer = F, inner = F)){
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  
  #initialize
  epoch = 0
  delta = Inf
  deltaU = Inf
  deltaV = Inf
  loss = Inf
  info = list(outer = c(), inner = c())
  if(!is.null(lambdax)) Kx = sparseMatrix(i = 1:p, j = 1:p, x = lambdax, dims = c(p, p))
  else Kx = NULL
  if(!is.null(lambday)) Ky = sparseMatrix(i = 1:q, j = 1:q, x = lambday, dims = c(q, q))
  else Ky = NULL
  init = NULL
  u0 = rep(1, p)
  v0 = rep(1, q)
  
  while(delta > eps$outer & epoch < maxiter$outer){
    epoch = epoch + 1
    loss0 = loss
    
    impute = CCAimpute(X, Y, Kx = Kx, Ky = Ky, init = init, k = k, eps = eps$inner, maxiter = maxiter$inner, verbose = verbose$inner)
    XU = impute$X %*% impute$U
    YV = impute$Y %*% impute$V
    loss = diag(cor(XU, YV))[1]
    
    info$inner = rbind(info$inner, data.frame(epoch = epoch, impute$info)) %>% data.frame()
    info$outer = rbind(info$outer, c(epoch, loss, delta)) %>% data.frame()
    
    #update kernels
    kx = lambdax/pmax(abs(impute$U[,1]), 1e-6)
    #kx[abs(impute$U[,1]) < 1e-6] = 0
    ky = lambday/pmax(abs(impute$V[,1]), 1e-6)
    #ky[abs(impute$V[,1]) < 1e-6] = 0
    Kx = diag(kx)
    Ky = diag(ky)
    init = list(X = impute$X, Y = impute$Y)
    u = impute$U[,1]
    v = impute$V[,1]
    
    if(epoch > 1){ 
      delta = abs((loss0 - loss)/loss0)
      deltaU = sum((u0 - u)^2)/sum(u0^2)
      deltaV = sum((v0 - v)^2)/sum(v0^2)
    }
    u0 = u
    v0 = v
    if(verbose$outer) cat("===================================\nepoch", epoch, "cor", loss, "delta", delta, "deltaU", deltaU, "deltaV", deltaV, "\n===================================\n\n")
    
  }
  colnames(info$outer) = c("epoch", "cor", "delta")
  return(list(U = impute$U, V = impute$V, X = impute$X, Y = impute$Y, info = info))
}

