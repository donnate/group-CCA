############## Iteratively re-weighted CCCA #################

D = function(alpha){
  as.matrix(dist(alpha, diag = T, upper = T))
}

kernel = function(alpha, lambda, gamma, A){
  W = lambda * A * (gamma / pmax(D(alpha), 1e-6) + (1 - gamma))
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

my_CCA = function(X, Y, KX, KY){
  index = which()
  Cxx = cov(X) + KX
  Cyy = cov(Y) + KY
  Cxy = cov(X, Y, use = "pairwise")
  cc = fda::geigen(Cxy, Cxx, Cyy)
  names(cc) = c("cor", "xcoef", "ycoef")
  list(alpha = cc$xcoef[,1:r], beta = cc$ycoef[,1:r], rho = cc$cor[1:r])
}


CCCA = function(X, Y, lambda, gamma, A, mu, maxiter = 100, eps = 1e-6, verbose = FALSE){
  #lambda: penalty magnitude
  #gamma: if 1 then l1, if 0 then l2
  #A: weights in the penalty
  p = ncol(X)
  q = ncol(Y)
  iter = 0
  alpha = CCA::cc(X[,1:5], Y[,1:5])$xcoef[, 1:r]
  alpha = rbind(alpha, matrix(0, p-5, r))
  beta = CCA::cc(X[,1:5], Y[,1:5])$ycoef[, 1:r]
  beta = rbind(beta, matrix(0, p-5, r))
  cc = list(alpha = alpha, beta = beta)
  #cor0 = rho(X, Y, cc$alpha, cc$beta, kernel(alpha, lambda, gamma, A), diag(mu, q))
  cor0 = sum((X %*% cc$alpha - Y %*% cc$beta)^2)
  delta = Inf
  if(verbose) cat("iter ", iter, "cor", cor0, "\n")
  
  while(delta > eps & iter < maxiter){
    iter = iter + 1
    if (gamma ==  "21" ){
      all_norms = sqrt(apply(alpha^2, 1, sum))
      all_norms_inv = sapply(all_norms, function(x){ifelse(x<1e-5, 0, 1/x)})
      Kx = lambdax * diag(all_norms_inv) 
      all_norms = sqrt(apply(beta^2, 1, sum))
      all_norms_inv = sapply(all_norms, function(x){ifelse(x<1e-5, 0, 1/x)})
      Ky = lambday * diag(all_norms_inv) 
    }else{
      if (gamma ==  "22" ){
        Kx = diag(1, p)
        Ky = diag(1, q) 
      }else{
        Kx = NULL
        Ky = NULL
      }
    }

    cc = my_CCA(X, Y, Kx, Ky)
    alpha = cc$alpha
    beta = cc$beta
    #cor =  rho(X, Y, cc$alpha, cc$beta, kernel(alpha, lambda, gamma, A), diag(mu, q))
    cor = mean((X %*% cc$alpha - Y %*% cc$beta)^2)
    delta = abs((cor - cor0)/cor0)
    cor0 = cor
    if(verbose) cat("iter ", iter, "cor", cor, "delta", delta, "\n")
  }
  return(cc)
}


CCCA_gridsearch = function(Xtrain, Ytrain, Xtest, Ytest, method, subtype = "", lambdas, gammas, A, mus, eps = 1e-6, verbose = F){
  p = ncol(Xtrain)
  q = ncol(Ytrain)
  result = c()
  if(method == "RCCA") gammas = 0
  cat("\n\nmethod:", method, subtype)
  for(lambda in lambdas){
    cat("\nlambda:", lambda)
    for(mu in mus){
      cat(" mu:", mu)
      for(gamma in gammas){
        if(method == "RCCA") cc = CCA(Xtrain, Ytrain, KX = diag(lambda, p), KY = diag(mu, q))
        if(method == "CCCA") cc = CCCA(Xtrain, Ytrain, lambda, gamma, A, mu, eps = eps, verbose = verbose)
        result = rbind(result, data.frame(cors = c(cor(Xtrain %*% cc$alpha, Ytrain %*% cc$beta), cor(Xtest %*% cc$alpha, Ytest %*% cc$beta)),
                                          mses = c(mean((Xtrain %*% cc$alpha - Ytrain %*% cc$beta)^2), mean((Xtest %*% cc$alpha - Ytest %*% cc$beta)^2)),
                                          type = c("train", "test"), 
                                          method = paste(method, subtype), 
                                          lambda, gamma, mu))
      }
    }
  }
  return(result)
}

gen = generate_example_non_trivial_pca(n, p, q,
                                       r_pca = r_pca,
                                       nnzeros = nnzeros,
                                       noise = s,
                                       theta = diag(seq(0.9, 0.7, length.out = r)),
                                       lambda_pca = 1,
                                       r = r, 
                                       overlapping_amount = overlapping_amount,
                                       normalize_diagonal = FALSE,
                                       prop_missing = prop_missing) 
