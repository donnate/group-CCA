library(dplyr)
library(tidyr)
library(Matrix)

transform = function(Y, Ky){
  mu = colMeans(Y)
  Sy = cov(Y)
  if(!is.null(Ky)) Sy = Sy + Ky 
  ED = eigen(Sy, symmetric = T)
  r = sum(ED$values > 1e-6)
  Sigma_sqrt_inv = ED$vectors[,1:r] %*% diag(1/sqrt(ED$values[1:r])) %*% t(ED$vectors[,1:r])
  Sigma_sqrt = ED$vectors[,1:r] %*% diag(sqrt(ED$values[1:r])) %*% t(ED$vectors[,1:r])
  Sigma_inv = ED$vectors[,1:r] %*% diag(1/(ED$values[1:r])) %*% t(ED$vectors[,1:r])
  return(list(Sigma_sqrt_inv = Sigma_sqrt_inv, Sigma_sqrt = Sigma_sqrt, 
              Sigma_inv=Sigma_inv, mu = mu))
}

RRR = function(X, Y, Kx, k, highdim){
  # solve RRR: ||Y-XB|| + tr(Bt K B)
  n = nrow(X)
  if(!highdim){
    Sx = t(X) %*% X
    if(!is.null(Kx)) Sx = Sx + (n-1) * as.matrix(Kx) 
    Sxy = t(X) %*% Y
    Bhat = solve(Sx) %*% Sxy
    ED = eigen(t(Sxy) %*% Bhat, symmetric = T)
  } else {
    if(sum(Kx) == sum(diag(Kx))){
      A = t(scale(X, center = F, scale = diag(Kx)))
      A[diag(Kx) == 0,] = 0
    } else A = as.matrix(solve(Kx, t(X)))
    Sx = X %*% A
    Bhat = A %*% solve(Sx + diag(n-1, n)) %*% Y
    Sxy = t(X) %*% Y
    ED = eigen(t(Sxy) %*% Bhat, symmetric = T)
  }
  V = ED$vectors[, 1:k, drop = F]
  U = Bhat %*% V
  Yhat = (X %*% U) %*% t(V)
  loss = mean((Y - Yhat)^2)
  return(list(U = U, V = V, Yhat = Yhat, loss = loss))
}

CCA_step = function(X, Y, Kx, Ky, k, highdim = FALSE){
  tr = transform(Y, Ky)
  #Y = scale(Y, center = tr$mu, scale = F) # scaling is taken care of in the next step
  Ytilde = Y %*% tr$Sigma_sqrt_inv
  rrr = RRR(X, Ytilde, Kx, k, highdim)
  #for imputation
  Yhat = rrr$Yhat %*% tr$Sigma_sqrt 
  #Yhat = scale(Yhat, center = -tr$mu, scale = F)
  #cca solution
  U = rrr$U
  U = U %*% sqrtm(t(U) %*% cov(X) %*% U)$Binv
  V = tr$Sigma_inv %*% rrr$V #Sigma_Y^{-1/2}V
  V = V %*% sqrtm(t(V) %*% cov(Y) %*% V)$Binv 
  return(list(U = U, V = V, Yhat = Yhat, rrrloss = rrr$loss))
}

CCAimpute = function(X, Y, Kx = NULL, Ky = NULL, k = min(nrow(X), min(ncol(X), ncol(Y))), 
                     init = NULL, eps = 1e-6, maxiter = 1000, verbose = F){
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  X = scale(X, scale = FALSE)
  Y = scale(Y, scale = FALSE)
  if(p > n & is.null(Kx)) stop("Regularization of X coefficients is required when p > n")
  if(q > n & is.null(Ky)) stop("Regularization of Y coefficients is required when q > n")
  
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
    X = data.frame(X) %>% mutate_all(~replace_na(., median(., na.rm = TRUE))) %>% as.matrix()
    Y = data.frame(Y) %>% mutate_all(~replace_na(., median(., na.rm = TRUE))) %>% as.matrix()
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
  sx = t(XU) %*% XU /(n-1)
  if(!is.null(Kx)) sx = sx + diag(t(U) %*% Kx %*% U)
  sy = t(YV) %*% YV /(n-1)
  if(!is.null(Ky)) sy = sy + diag(t(V) %*% Ky %*% V)
    
  return(list(U = U %*% sqrtm(sx)$Binv, 
              V = V%*% sqrtm(sy)$Binv, 
              X = X, Y = Y, info = data.frame(info)))
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




