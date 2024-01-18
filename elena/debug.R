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



CCA_rrr_impute = function(X, Y, Sx=NULL, Sy=NULL,
                          lambda =0, Kx, r, highdim=FALSE, 
                          penalty = "l21", lambda_Kx=0, solver="rrr",
                          LW_Sy = FALSE){
  # solve RRR: ||Y-XB|| + tr(Bt K B)
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  ### start filling in missing values
  if(q >n){
    X_temp <- X
    X <- Y
    Y <- X_temp
  }
  missY = is.na(Y)
  missX = is.na(X)
  
  X = data.frame(X) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))) %>% as.matrix()
  Y = data.frame(Y) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))) %>% as.matrix()
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
  while( (it <max_it) & (delta > 1e-3)){
    if(!highdim){
      ### iterate until converge
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
        if (solver =="CVX"){
          print("Using CVXR")
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
        print(t(U) %*% Sx %*% U)
        print(t(V) %*% Sy %*% V)
      }
      X[missX] = rep(apply(X_orig, 2, mean), nrow(X[missX])) +  (U %*% Lambda %*% t(V) )[missX]
      Y[missY] = rep(apply(Y_orig, 2, mean), nrow(Y[missY])) +  (V %*% Lambda %*% t(U) )[missY]
    }
    
  }
  loss = mean((Y %*% V - X %*% U)^2)
  return(list(U = U, V = V, loss = loss))
}

