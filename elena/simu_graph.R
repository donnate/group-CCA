


get_edge_incidence <- function(g, weight = 1){
  n_nodes = vcount(g)
  d_max = max(degree(g))
  #d_max = 1
  edges = data.frame(as_edgelist(g)) %>%
    arrange(X1, X2)
  W = matrix(0, nrow = n_nodes, ncol=n_nodes)
  Gamma = matrix(0, nrow(edges), n_nodes)
  
  weights = rep(1, nrow(edges))
  edges["weights"]= weights
  # Make beta_v into a matrix
  for (e in 1:nrow(edges)){
    W[edges$X1[e], edges$X2[e]] = edges$weights[e] #* (beta_v[edges$X2[e]] + beta_v[edges$X1[e]])/2
    W[edges$X2[e], edges$X1[e]] = edges$weights[e] #* (beta_v[edges$X2[e]] + beta_v[edges$X1[e]])/2
    Gamma[e,edges$X1[e]] = edges$weights[e] 
    Gamma[e,edges$X2[e]] = - edges$weights[e]
  }
  return(list(Gamma=Gamma,
              Gamma_dagger = pinv(Gamma),
              W = W))
}

### gamma sparse CCA
CCA_rrr_graph = function(X, Y, A, Sx=NULL, Sy=NULL,
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
  ### turn adjacency matrix into edge incidence matrix
  
  
  ###
  X_transf = X %*% Gamma
  
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
  Sxy = t(X_transf) %*% tilde_Y/ n 
  test <- cv.srrr( tilde_Y, X_transf, nrank = r,
                   method ="glasso",
                   nfold = nfolds, norder = NULL,
                   A0 = NULL,   V0 = NULL,
                   modstr = list("lamA" = rep(lambda, 10),
                                 "nlam" = 10))
  B_opt <- test$coef # test$U  %*% test$D %*% t(test$V)
  B_opt[which(abs(B_opt)<1e-5)] = 0
  
  B_opt = Gamma_dagger %*% B_opt #### transform it back
  #### Make it sparse
  B_opt_intercept = apply(B_opt, 2, mean)
  B_opt = B_opt - diag(B_opt_intercept) * matrix(1, p, q)
  I = which(apply(B_opt^2, 1, sum) > 0)
  #sol = svd(t(B_OLS[I, ]), nu = r, nv=r)
  
  if (length(I) >0 & length(I) < n){
    lw_covX <- corpcor::cov.shrink(X[,I])
    Sx_loc<- as.matrix(lw_covX)
    svd_Sx = svd(Sx_loc)
    sqrt_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, sqrt(x), 0)}))  %*% t(svd_Sx$u)
    sqrt_inv_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sx$u)
    sol = svd(sqrt_Sx %*% B_opt[I,])
    #sol = svd(t(B_OLS[I, ]), nu = r, nv=r)
    V = sqrt_inv_Sy %*% sol$v[, 1:r]
    U = matrix(0, p, r)
    U[I,] = sqrt_inv_Sx %*% sol$u[, 1:r]
  }else{
    if(length(I) == 0){
      U = matrix(0, p, r)
      V = sqrt_inv_Sy %*% test$V[, 1:r] 
    }else{
      #### Solve another regression problem
      sol = svd(B_opt[I,])
      V = sqrt_inv_Sy %*% sol$v[, 1:r]
      B_optV = B_opt %*% t(sol$v[, 1:r])
      ###### Solve another group lasso problem
      # Create a block diagonal matrix for X
      Y_long <- as.vector(B_optV)
      X_block <- kronecker(diag(1, r), X %*% Gamma)
      groups_long <- rep(1:p, each = r)
      #### run the CV procedure)
      test = cv.gglasso(X_block, Y_long, group=groups_long, lambda = lambdax,
                        intercept = FALSE)
      ind_lambda = which(test$gglasso.fit$df >0)
      ind = which.min(test$cvm[ind_lambda])
      beta = test$gglasso.fit$beta[, ind_lambda[ind]]
      U <- Gamma_dagger %*% matrix(beta, nrow = p, ncol = q, byrow = TRUE)
      U <- U %*% sqrtm(t(U) %*% Sx %*% U)$Binv
    }
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
  
  loss = mean((Y %*% V - X %*% U)^2)
  return(list(U = U, V = V, loss = loss))
}




CCA_rrr_graph.CV<- function(X, Y, 
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









