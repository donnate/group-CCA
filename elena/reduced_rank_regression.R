library(dplyr)
library(tidyr)
library(Matrix)
library(glmnet)
library(gglasso)

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




RRR_CCA_group <- function(X, Y,  groups, Sy=NULL, 
                          Sx=NULL, r=2, 
                          lambdax =0, LW_Sy = TRUE){
  
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
  
  if(is.null(Sx)){
    Sx = t(X) %*% X /n
  }
  if (is.null(Sy)){
    Sy = t(Y) %*% Y /n
    if (LW_Sy){
      lw_cov <- corpcor::cov.shrink(Y)
      Sy <- as.matrix(lw_cov)
    }
  }
  
  
  
  # Create a block diagonal matrix for X
  X_block <- kronecker(diag(1, q), X)
  
  groups_long <- rep(groups, each = q)
  # The dimension of B is now pq x q
  svd_y = svd(Sy)
  Sigma_y_sqrt_inv = svd_y$u %*% diag(sapply(svd_y$d, function(x){ifelse(x >1e-4, 1/sqrt(x), 0)})) %*% t(svd_y$u)
  Y_tilde = Y %*% Sigma_y_sqrt_inv
  if(is.null(group_intercept) ==FALSE){
    Y_tilde = Y_tilde
  }
  Y_long <- as.vector(Y_tilde)
  #### run the CV procedure)
  test = gglasso(X_block, Y_long, group=groups_long, lambda = c(lambdax),
                 intercept = FALSE)
  
  beta = test$gglasso.fit$beta[, ind_lambda[ind]]
  B_opt <- matrix(beta, nrow = p, ncol = q, byrow = TRUE)
  B_opt[which(abs(B_opt)<1e-5)] = 0
  I = which(apply(B_opt^2, 1, sum) > 0)
  if (length(I) >0){
    svd_Sx = svd(Sx[I, I])
    sqrt_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, sqrt(x), 0)}))  %*% t(svd_Sx$u)
    sqrt_inv_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sx$u)
    sol = svd(sqrt_Sx %*% B_opt[I,])
    #sol = svd(t(B_OLS[I, ]), nu = r, nv=r)
    V = Sigma_y_sqrt_inv %*% sol$v[, 1:r]
    U = matrix(0, p, r)
    U[I,] = sqrt_inv_Sx %*% sol$u[, 1:r]
  }else{
    U = matrix(0, p, r)
    V = Sigma_y_sqrt_inv %*% test$V[, 1:r]
  }
  loss = mean((Y %*% V - X %*% U)^2)
  return(list(U = U, V = V, loss = loss))
}

cv.RRR_CCA_group <- function(X, Y, groups, Sy = NULL,
                             Sx = NULL,
                             r=2, lambdax =10^seq(-5, 1,length.out=50),
                             group_intercept = NULL,
                             LW_Sy = TRUE){
  
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
  
  if(is.null(Sx)){
    Sx = t(X) %*% X /n 
  }
  if(is.null(Sy)){
    Sy = t(Y) %*% Y /n
    if (LW_Sy){
      lw_cov <- corpcor::cov.shrink(Y)
      Sy <- as.matrix(lw_cov)
    }
  }
  
  
  
  # Create a block diagonal matrix for X
  X_block <- kronecker(diag(1, q), X)
  
  groups_long <- rep(groups, each = q)
  # The dimension of B is now pq x q
  svd_y = svd(Sy)
  Sigma_y_sqrt_inv = svd_y$u %*% diag(sapply(svd_y$d, function(x){ifelse(x >1e-4, 1/sqrt(x), 0)})) %*% t(svd_y$u)
  Y_tilde = Y %*% Sigma_y_sqrt_inv
  Y_long <- as.vector(Y_tilde)
  #### run the CV procedure)
  test = cv.gglasso(X_block, Y_long, group=groups_long, lambda = lambdax,
                    intercept = FALSE)
  ind_lambda = which(test$gglasso.fit$df >0)
  ind = which.min(test$cvm[ind_lambda])
  beta = test$gglasso.fit$beta[, ind_lambda[ind]]
  B_opt <- matrix(beta, nrow = p, ncol = q, byrow = TRUE)
  B_opt[which(abs(B_opt)<1e-5)] = 0
  I = which(apply(B_opt^2, 1, sum) > 0)
  if (length(I) >0){
    svd_Sx = svd(Sx[I, I])
    sqrt_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, sqrt(x), 0)}))  %*% t(svd_Sx$u)
    sqrt_inv_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sx$u)
    sol = svd(sqrt_Sx %*% B_opt[I,])
    #sol = svd(t(B_OLS[I, ]), nu = r, nv=r)
    V = Sigma_y_sqrt_inv %*% sol$v[, 1:r]
    U = matrix(0, p, r)
    U[I,] = sqrt_inv_Sx %*% sol$u[, 1:r]
  }else{
    U = matrix(0, p, r)
    V = Sigma_y_sqrt_inv %*% test$V[, 1:r]
  }
  loss = mean((Y %*% V - X %*% U)^2)
  return(list(U = U, V = V, loss = loss))
}

