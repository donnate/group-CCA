CCA_RRR = function(X, Y, Kx = NULL, Ky = NULL, 
                   lambdax = 0, lambday = 0, reg = "none",
                   k = min(nrow(X), min(ncol(X), ncol(Y))), 
                   init = NULL, eps = 1e-6, maxiter = 1000, verbose = F){
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  X = scale(X, scale = FALSE)
  Y = scale(Y, scale = FALSE)
  if(p > n & is.null(Kx)) stop("Regularization of X coefficients is required when p > n")
  if(q > n & is.null(Ky)) stop("Regularization of Y coefficients is required when q > n")

  
  #initialize
  iter = 0
  delta = Inf
  loss = Inf
  info = c()
  
  while(delta > eps & iter < maxiter){
    iter = iter + 1
    loss0 = loss
    
    #CCA step (impute Y)
    if(p > n){
        step = CCA_step(X, Y, Kx, Ky, k, highdim = T)
        if(verbose) cat("High dim on X\n")
      } 
      else{
        step = CCA_step(X, Y, Kx, Ky, k)
        U = step$U
        V = step$V
    }
    
    #update kernels
    if (reg %in% c("X", "both")){
      kx = lambdax/(apply(U^2, 1, sum))
      #kx[abs(impute$U[,1]) < 1e-6] = 0
      kx = lambday/sqrt(kx) 
      Kx = diag(kx)
    }else{
      Kx = diag(1, p)
    }
    if (reg %in% c("Y", "both")){
      ky = lambdax/(apply(V^2, 1, sum))
      #x[abs(impute$U[,1]) < 1e-6] = 0
      ky = lambday/sqrt(ky) 
      Ky = diag(ky)
    }else{
      Ky = diag(1, q)
    }
    
    # 
    # #CCA step (impute X)
    # if(rrrside %in% c("X", "both")){
    #   if(verbose) cat("RRR with response X\n")
    #   if(q > n){
    #     step = CCA_step(Y, X, Ky, Kx, k, highdim = T)
    #     if(verbose) cat("High dim on Y\n")
    #   } 
    #   else step = CCA_step(Y, X, Ky, Kx, k)
    #   X[missX] = step$Yhat[missX]
    #   imploss = mean((X - step$Yhat)^2)
    #   U = step$V
    #   V = step$U
    # }
    
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