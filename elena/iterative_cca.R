library(Matrix)
library(CVXR)
source("elena/init_cca.R")
get.subset <-function (x, method = c("theory", "method"), alpha.method = 0.05,
            alpha.theory = 1.5, sigma = NA, df = Inf) {
    if (method == "theory") {
      ans <- which(x/sigma^2/df > (1 + alpha.theory * sqrt(log(df)/df)))
    }
    if (method == "method") {
      ans <- holm.robust(x = x, alpha = alpha.method)
    }
    ans
}

# init_cca <- function(X, Y, r, lambdax, lambday, Gamma_u, Gamma_v){
#   
#   n <- dim(X)[1]
#   p <- dim(X)[2]
#   q <- dim(Y)[2]
#   if (lambdax >0  & is.null(Gamma_u)){
#     rownorm2 <- apply(atanh(cov(X, Y))^2, 1, sum) *n
#     all_I_p = p.adjust(1-pchisq(rownorm2, q))
#     I0 = which(all_I_p < 0.1)
#   }else{
#     I0  = 1:p
#     if (is.null(Gamma_u) == FALSE & lambdax >0){
#       L = t(Gamma_u) %*% Gamma_u
#       L <- L -diag(diag(L))
#       A = -L 
#       X = A %*% X
#     }
#   }
#   
#   if (lambday >0  & is.null(Gamma_v)){
#     colnorm2 <- apply(atanh(cov(X, Y))^2, 2, sum) * n
#     all_J_p = p.adjust(1-pchisq(colnorm2, p))
#     J0 = which(all_J_p < 0.1)
#   }else{
#     J0  = 1:q
#     if (is.null(Gamma_v) == FALSE & lambday >0){
#       L = t(Gamma_v) %*% Gamma_v
#       L <- L -diag(diag(L))
#       A = -L 
#       Y = A %*% Y
#     }
#   }
#   
#   if(length(I0)<2){
#     I0 = 1:p
#   }
#   if(length(J0)<2){
#     J0 = 1:p
#   }
#   
#   init = CCA::cc(X[, I0], Y[, J0])
#   U = matrix(0, p, r)
#   U[I0, ] <- init$xcoef[,1:r]
#   V = matrix(0, q, r)
#   V[J0, ] <- init$ycoef[,1:r]
#   return(list(U=U, V=V))
# 
# }

alternating_cca <- function(X, Y, r, init_coef = NULL, lambdax = 0,
                            lambday=0, Gamma_u=NULL, Gamma_v=NULL,
                            thres=1e-5, max_iter=100){
  # Initializations
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  # Check for missingness
  missY = is.na(Y)
  missX = is.na(X)
  
  if(max(missX) >0){
    X = data.frame(X) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))) %>% as.matrix()
  } 
  if(max(missY) >0){
    Y = data.frame(Y) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))) %>% as.matrix()
  } 
  
  if (is.null(init_coef)){
    init_coef <- init_cca(X,Y, r, penalize_X=TRUE,
                          penalize_Y=TRUE)
  }
  U <- init_coef$u_init
  V <- init_coef$v_init
  ZV <- Y %*% V
  ZU <- X %*% U
  U.old = U
  V.old = V
  U_init = U
  V_init = V
  it = 1
  eps = 1e10
  while(it<max_iter & eps > thres) {
    if (lambdax >0){
      U <- Variable(p, r)
      if (is.null(Gamma_u)){
        objective <- Minimize(1/n * sum_squares(X %*% U - ZV) + lambdax * sum(norm2(U, axis=1)))
        problem <- Problem(objective)
        result <- solve(problem)
        U <- result$getValue(U)
        U_init <- U
      }else{
        objective <- Minimize(1/n * sum_squares(X %*% U - ZV) + lambdax * sum(norm2(Gamma_u %*% U, axis=1)))
        problem <- Problem(objective)
        result <- solve(problem)
        U <- result$getValue(U)
        U_init <- U
      }

    }else{
      U = solve(t(X) %*%X) %*% t(X) %*% ZV
    }
    #U = U %*% diag(1/sqrt(apply(U^2, 2, sum)))
    if (max(U) < 1e-3){
      ### Solution is essentially 0
      return(list(U=matrix(0, p, r), V = V))
    }else{
      norms = apply(U, 1, norm)
      ind = which(norms > 1e-3)
      if (length(ind)>1){
        U[ind,] = U[ind,] %*% sqrtm(t(U)[,ind] %*% cov(X[, ind], X[, ind]) %*% U[ind, ])$Binv
      }else{
        U[ind,] = U[ind,] %*% sqrtm( cov(X[, ind], X[, ind]) * (t(U)[,ind] %*% U[ind, ]))$Binv
      }
      
    }
    ZU <- X %*% U

    if (lambday >0){
      V <- Variable(q, r)
      
      if (is.null(Gamma_v)){
        objective <- Minimize(1/n * sum_squares(ZU - Y %*% V) + lambday * sum(norm2(V, axis=1)))
        problem <- Problem(objective)
        result <- solve(problem)
        V <- result$getValue(V)
        V_init <- V
      }else{
        objective <- Minimize(1/n * sum_squares(ZU - Y %*% V) + lambday * sum(norm2(Gamma_v %*% V, axis=1)))
        problem <- Problem(objective)
        result <- solve(problem)
        V <- result$getValue(V)
        V_init <- V
      }

    }else{
      V = solve(t(Y) %*%Y) %*% t(Y) %*% ZU
    }
    #V = V %*% diag(1/sqrt(apply(V^2, 2, sum)))
    if (max(V) < 1e-3){
      ### Solution is essentially 0
      return(list(U=U, V = matrix(0, q, r)))
    }else{
      norms = apply(V, 1, norm)
      ind = which(norms > 1e-3)
      if (length(ind)>1){
        V[ind,] = V[ind,] %*% sqrtm(t(V)[,ind] %*% cov(Y[, ind], Y[, ind]) %*% V[ind, ])$Binv
      }else{
        V[ind,] = V[ind,] %*% sqrtm( cov(Y[, ind], Y[, ind]) * (t(V)[,ind] %*% V[ind, ]))$Binv
      }
      
     # V = V %*% sqrtm(t(V) %*% cov(Y) %*% V)$Binv
    }
    
    #### Solve for Y and X
    if(sum(missY)>0){
     Y[missY] <- ((X %*% U) %*%  t(V) %*% pinv((V) %*% t(V)))[missY]
    }
    if(sum(missX)>0){
     X[missX] <- ((Y %*% V) %*%  t(U)  %*% pinv((U) %*% t(U)))[missX]
    }
    ZV <- Y %*% V
    eps = subdistance( U, U.old) + subdistance(V, V.old)
    U.old = U
    V.old = V
    it = it + 1
    print(c(it, eps))
  }
  print(it)
  return(list(U=U, V = V))
}


pipeline_alternating_CCA<- function(X, Y, 
                             r=2,
                             param_lambda=c(1, 0.1, 0.01, 0.001),
                             kfolds=10,
                             maxiter=5000, convergence=1e-4, 
                             eta=1e-3, expand2D = FALSE){
  p1 <- dim(X)[2]
  p2 <- dim(Y)[2]
  p <- p1 + p2;
  n <- nrow(X)
  pp <- c(p1,p2);
  Sigma_X = cov(X)
  Sigma_Y = cov(Y)
  
  init <- init_cca(X,Y, r, penalize_X=TRUE,
                   penalize_Y=TRUE)
  
  if (expand2D){
    resultsx <- expand.grid(lambdax = param_lambda, lambday = param_lambda) %>%
      mutate(rmse = map2_dbl(lambdax, lambday,  ~cv_function_alternating_cca(X, Y, 
                                                                             kfolds=kfolds, init,
                                                                             lambdax= .x,
                                                                             lambday = .y,
                                                                             r=r, 
                                                                             maxiter=maxiter, 
                                                                             convergence=convergence)))
    
    resultsx = resultsx %>% filter(rmse >0) 
    opt_lambdax <- resultsx$lambdax[which.min(resultsx$rmse)]
    opt_lambday <- resultsx$lambday[which.min(resultsx$rmse)]
  }else{
    resultsx <- expand.grid(lambdax = param_lambda) %>%
      mutate(rmse = map_dbl(lambdax, ~cv_function_alternating_cca(X, Y, 
                                                                  kfolds=kfolds, init,
                                                                  lambdax= .x,
                                                                  lambday = .x,
                                                                  r=r, 
                                                                  maxiter=maxiter, 
                                                                  convergence=convergence)))
    
    resultsx = resultsx %>% filter(rmse >0) 
    opt_lambdax <- resultsx$lambdax[which.min(resultsx$rmse)]
    opt_lambday <- resultsx$lambdax[which.min(resultsx$rmse)]
  }

  



  
  
  #X, Y, Mask, kfolds=5, ainit,lambda, r=2, k=20,  
  #maxiter=1000, eta=0.001, convergence=1e-
  
  #print(resultsx)
  ###### (X, Y, Mask, kfolds=5, ainit, lambda, k=20)
  
  # print best hyperparameters and corresponding RMSE
  

  #which_lambdax = which(abs(resultsx$rmse-min(resultsx$rmse))/(1e-6  + min(resultsx$rmse)) <0.01)
  #lambda = max(resultsx$lambda[which_lambdax])
  #k = max(resultsx$k[which_lambdax])
  print(c("selected", opt_lambdax, opt_lambday))
  
  final <- alternating_cca(X, Y, r, 
                           init_coef = init, lambdax = opt_lambdax, 
                           lambday = opt_lambday,
                           Gamma_u=NULL, Gamma_v=NULL,
                           thres=convergence, 
                           max_iter=maxiter)
  print(final)
  if (max(abs(final$U)) > 1e-5){
    Uhat = final$U %*% sqrtm(t(final$U) %*% Sigma_X %*% final$U)$Binv 
  }else{
    Uhat = final$U
  }
  if (max(abs(final$V)) > 1e-5){
    Vhat = final$V %*% sqrtm(t(final$V) %*% Sigma_Y %*% final$V)$Binv 
  }else{
    Vhat = final$V
  }
  print(resultsx)
  return(list( ufinal = Uhat, 
               vfinal = Vhat,
               initu=init$u_init, 
               initv=init$v_init,
               lambdax=opt_lambdax,
               lambday=opt_lambday,
               resultsx=resultsx
  ))
  
}


cv_function_alternating_cca<- function(X, Y, kfolds=5, init,
                           lambdax=0.01,
                           lambday =0.01,
                           r=2, 
                           maxiter=20, 
                           convergence=1e-4) {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- rep(1e10, kfolds)
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
        final = alternating_cca(X_train, Y_train, r, 
                                init_coef = init, lambdax = lambdax, 
                                lambday = lambday,
                                Gamma_u=NULL, Gamma_v=NULL,
                                thres=convergence, 
                                max_iter=maxiter)
        # make predictions on validation data
        # compute RMSE on validation data
        rmse[i] <- mean((X_val %*% final$U - Y_val%*% final$V)^2)
        print(rmse)
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
