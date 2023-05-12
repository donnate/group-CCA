library(fda)


############## RRR #################

CCA_step = function(X, Y, k){
  #whiten
  ED = eigen(cov(Y))
  tr = ED$vectors %*% diag(1/sqrt(ED$values)) %*%t(ED$vectors)
  invtr = ED$vectors %*% diag(sqrt(ED$values)) %*%t(ED$vectors)
  Ytilde = Y %*% tr
  
  #solve RRR
  Bhat = solve(t(X) %*% X) %*% t(X) %*% Ytilde
  Yhat = X %*% Bhat 
  SVD = svd(Yhat)
  V = SVD$v[, 1:k, drop = F]
  U = Bhat %*% V
  return(list(U = U, V = invtr %*% V, Ucca = U, Vcca = tr %*% V))
}

CCA_RRR = function(X, Y, k = min(ncol(X), ncol(Y)), eps = 1e-4, maxiter = 1000){
  p = ncol(X)
  q = ncol(Y)
  
  maskY = !is.na(Y)
  maskX = !is.na(X)
  
  #initialize
  iter = 0
  delta = Inf
  loss = Inf
  X = data.frame(X) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))) %>% as.matrix()
  Y = data.frame(Y) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))) %>% as.matrix()
  
  while(delta > eps & iter < maxiter){
    iter = iter + 1
    
    loss0 = loss
    upd = CCA_step(Y, X, k)
    #impute X
    if(sum(1 - maskX) > 0) X[!maskX] = (Y %*% upd$U %*% t(upd$V))[!maskX] 
    
    upd = CCA_step(X, Y, k)
    #impute Y
    if(sum(1 - maskY) > 0) Y[!maskY] = (X %*% upd$U %*% t(upd$V))[!maskY]
    
    Ucca = upd$Ucca 
    Vcca = upd$Vcca 
    loss = sum(diag(cor(X %*% Ucca, Y %*% Vcca)))
    
    if(iter > 1) delta = abs((loss0 - loss)/loss0)
    cat("iter ", iter, "loss", loss, "delta", delta, "\n")
  }
  return(list(U = Ucca, V = Vcca, cors = diag(cor(X %*% Ucca, Y %*% Vcca))))
}  


############## evaluation #################

principal_angles = function(A, B){
  angles = rep(0, ncol(A))
  QRa = qr(A)
  QRb = qr(B)
  Qa = qr.Q(QRa)
  Qb = qr.Q(QRb)
  
  d = svd(t(Qa) %*% Qb)$d
  if(QRa$rank <= QRb$rank) Q = Qb - Qa %*% (t(Qa) %*% Qb)
  else Q = Qa - Qb %*% (t(Qb) %*% Qa)
  s = sort(svd(Q)$d)
  
  for(i in 1:min(QRa$rank, QRb$rank)){
    if(d[i]^2 < 0.5) angles[i] = acos(d[i])
    else if(s[i]^2 <= 0.5) angles[i] = asin(s[i])
  }
  angles
}

angles = function(A, B){
  p = ncol(A)
  angs = rep(0, p)
  for(i in 1:p){
    angs[i] = principal_angles(A[,1:i,drop = F], B[,1:i,drop = F])[1]
  } 
  angs
}

evaluate = function(X, Y, U, V, U0, V0){
  data.frame(comp = 1:ncol(V), cors = diag(cor(X %*% U, Y %*% V)),
             mses = colMeans((X %*% U - Y %*% V)^2),
             Uangs = angles(U, U0), 
             Vangs = angles(V, V0),
             XUangs = angles(X %*% U, X %*% U0), 
             YVangs = angles(Y %*% V, Y %*% V0))
}

############## Iteratively re-weighted CCCA #################

D = function(alpha){
  as.matrix(dist(alpha, diag = T, upper = T))
}

kernel = function(alpha, lambda, gamma, A){
  W = lambda * A * (gamma/pmax(D(alpha), 1e-6) + (1 - gamma))
  w = rowSums(W)
  return(2 * (diag(w) - W))
}

rho = function(X, Y, alpha, beta, KX, KY){
  alpha = as.matrix(alpha)
  beta = as.matrix(beta)
  Cxx = var(X, use = "pairwise") + KX
  Cyy = var(Y, use = "pairwise") + KY
  Cxy = cov(X, Y, use = "pairwise")
  return(c(t(alpha) %*% Cxy %*% beta/sqrt(t(alpha) %*% Cxx %*% alpha) / sqrt(t(beta) %*% Cyy %*% beta)))
}

CCA = function(X, Y, KX, KY){
  Cxx = var(X, use = "pairwise") + KX
  Cyy = var(Y, use = "pairwise") + KY
  Cxy = cov(X, Y, use = "pairwise")
  cc = geigen(Cxy, Cxx, Cyy)
  names(cc) = c("cor", "xcoef", "ycoef")
  list(alpha = cc$xcoef[,1], beta = cc$ycoef[,1], rho = cc$cor[1])
}


CCCA = function(X, Y, lambda, gamma, A, mu, maxiter = 100, eps = 1e-6, verbose = FALSE){
  #lambda: penalty magnitude
  #gamma: if 1 then l1, if 0 then l2
  #A: weights in the penalty
  p = ncol(X)
  q = ncol(Y)
  iter = 0
  alpha = rnorm(p)
  beta = rnorm(q)
  cc = list(alpha = alpha, beta = beta)
  cor0 = rho(X, Y, cc$alpha, cc$beta, kernel(alpha, lambda, gamma, A), diag(mu, q))
  delta = Inf
  if(verbose) cat("iter ", iter, "cor", cor0, "\n")
  while(delta > eps & iter < maxiter){
    iter = iter + 1
    cc = CCA(X, Y, kernel(cc$alpha, lambda, gamma, A), diag(mu, q))
    cor =  rho(X, Y, cc$alpha, cc$beta, kernel(alpha, lambda, gamma, A), diag(mu, q))
    delta = abs((cor - cor0)/cor0)
    cor0 = cor
    if(verbose) cat("iter ", iter, "cor", cor, "delta", delta, "\n")
  }
  return(cc)
}


CCCA_gridsearch = function(Xtrain, Ytrain, Xtest, Ytest, method, lambdas, gammas, A, mus, eps = 1e-6, verbose = F){
  p = ncol(Xtrain)
  q = ncol(Ytrain)
  result = c()
  if(method == "RCCA") gammas = 0
  for(lambda in lambdas){
    for(mu in mus){
      for(gamma in gammas){
        cat("\nmethod:", method, "lambda:", lambda, "gamma:", gamma, "mu:", mu, "\n")
        if(method == "RCCA") cc = CCA(Xtrain, Ytrain, KX = diag(lambda, p), KY = diag(mu, q))
        if(method == "CCCA") cc = CCCA(Xtrain, Ytrain, lambda, gamma, A, mu, eps = eps, verbose = verbose)
        result = rbind(result, data.frame(cors = c(cor(Xtrain %*% cc$alpha, Ytrain %*% cc$beta), cor(Xtest %*% cc$alpha, Ytest %*% cc$beta)),
                                          mses = c(mean((Xtrain %*% cc$alpha - Ytrain %*% cc$beta)^2), mean((Xtest %*% cc$alpha - Ytest %*% cc$beta)^2)),
                                          type = c("train", "test"), 
                                          method = ifelse(mean(A) == 1, paste(method, "full"), method), 
                                          lambda, gamma, mu))
      }
    }
  }
  return(result)
}
