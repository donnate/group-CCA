# Function for Thredholded Gradient Descent
# Inputs:
# =======
# A, B:         Pair of matrix to calculate generalized eigenspace (sigma, sigma0)
# r:            Latent dimension
# init:         Initialize estimator of generlized eigenspace, obtained from sgca_int. Need to be k-sparse
# lambda:       penalty parameter of Lagrangian function f(L), default to be 0.01
# convergence:  tolerance level for convergence in TGD
# maxiter:      maximum number of iterations in TGD
# plot:         if set to True will plot intermediate iterations, need to specify scale variable V (as shown in paper)
#               default set to False

# Outputs:
# ========
# final_estimate:   final estimation of leading r sparse generalized eigenspace
library(pracma)
source("elena/init_cca.R")

compute_grad <- function(ut, Gamma = NULL, Z = NULL){
  
  if (is.null(Gamma) == FALSE){
    Gu = Gamma %*% ut
    norm_gamma_u <- sqrt(apply(Gu^2, 1, sum))
    norm_gamma_u_inv <- sapply(norm_gamma_u, function(x){ifelse(x > 1e-5, 1/x, 0)})
    grad_fn <-  diag(norm_u_inv ) %*% t(Gamma) %*% Gu
    
  }else{
    if (is.null(Z) == FALSE){
      list_Z <- unique(Z)
      norm_u <-  apply(ut^2, 1, sum)
      norm_G <- sapply(1:length(list_Z), function(i){
        z = list_Z[i]
        ind = which(Z == z)
        return(sqrt(sum(norm_u[ind])))
      })
      grad_fn <- ut
      for (i in 1:length(list_Z)){
            ind = which(Z == list_Z[i])
            grad_fn[ind, ] <- diag(rep(ifelse(norm_G[i] > 1e-5), 
                                              1/norm_G[i] ,
                                              0 ), 
                                      length(ind)) %*% norm_u[ind, ]
      }
    }else{
      norm_u <- sqrt(apply(ut^2, 1, sum))
      norm_u_inv <- sapply(norm_u, function(x){ifelse(x > 1e-5, 1/x, 0)})
      grad_fn <-  diag(norm_u_inv ) %*% ut
    }
  }
  return(grad_fn)
}
  


cca_gd <- function(X, Y, r, init, k,  lambda_x = 0.01, 
                   lambda_y = 0.01, eta=0.01, convergence=1e-3,
                   penalize_X = TRUE, penalize_Y = TRUE, 
                   maxiter=10000, plot = FALSE,
                   Zu=NULL, A_u=NULL,
                   Zv=NULL, A_v=NULL,
                   Gamma_u = NULL,
                   Gamma_v = NULL){
    init <- init_cca(X,Y, r, penalize_X=penalize_X,
                     penalize_Y=penalize_Y,
                     Zu=Zu, A_u=A_u,
                     Zv=Zv, A_v=A_v)
    
    Sigma_X <- cov(X)
    Sigma_Y <- cov(Y)
    Sigma_XY <- cov(X, Y)
    ut <- init$u_init
    vt <- init$v_init
    criteria <- 1e10
    iter <- 0
    error <- rep(0, maxiter)
    # renormalization 
    while(iter <= maxiter){
      #perform gradient descent on X
      grad_fn <- compute_grad(ut, Gamma = Gamma_u, Z = Zu)
      grad <- 1/n * (Sigma_X %*% ut - Sigma_XY %*% vt)  + lambda_x * grad_fn;
      new_u <- ut - eta * grad
      # Apply Regularisation:
      new_u <- new_u %*% sqrtm(t(new_u) %*% Sigma_X %*% new_u)$Binv

      grad_fn <- compute_grad(vt, Gamma = Gamma_v, Z = Zv)
      grad <- 1/n * (Sigma_Y %*% vt - t(Sigma_XY) %*% new_u)  + lambda_y * grad_fn;
      new_v <- vt - eta * grad
      # Apply Regularisation:
      new_v <- new_v %*% sqrtm(t(new_v) %*% Sigma_Y %*% new_v)$Binv

      criteria <- mean(c( subdistance(ut, new_u),  subdistance(vt, new_v)))
      ut <- new_u
      vt <- new_v
      print(criteria)
      iter <- iter+1
    }
    if (plot){
      plot(error[1:iter], type='l',  main="Matrix distance and iterations", 
           xlab="Number of iterations", ylab="Matrix distance",)
    }
    final_estimate <- ut %*% sqrtm(t(ut) %*% B %*% ut)$Binv
    return(final_estimate)
}

cca_gd2 <- function(X, Y, r, init, kx=20, ky=20,
                    lambda_x = 0.01, 
                    lambda_y = 0.01, 
                    eta=0.01, 
                    convergence=1e-3,
                    penalize_X = TRUE, penalize_Y = TRUE,
                    lambda_gamma_x = 0.001,
                    lambda_gamma_y = 0.001,
                    maxiter=10000, plot = FALSE,
                    Zu=NULL, #A_u=NULL,
                   Zv=NULL, #A_v=NULL,
                   Gamma_u = NULL,
                   Gamma_v = NULL){
  #init <- init_cca(X,Y, r, penalize_X=penalize_X,
  #                 penalize_Y=penalize_Y,
  #                 Zu=Zu, A_u=A_u,
  #                 Zv=Zv, A_v=A_v)
  
  Sigma_X <- cov(X)
  Sigma_Y <- cov(Y)
  Sigma_XY <- cov(X, Y)
  ut <- init$u_init
  vt <- init$v_init
  criteria <- 1e10
  iter <- 0
  error <- rep(0, maxiter)
  # renormalization 
  while(iter <= maxiter & criteria > convergence){
    #perform gradient descent on X
    #grad_fn <- compute_grad(ut, Gamma = Gamma_u, Z = Zu)
    grad <- 1/n * (Sigma_X %*% ut - Sigma_XY %*% vt)  + lambda_x * Sigma_X %*% ut %*% (t(ut) %*% Sigma_X %*% ut- diag(r));
    new_u <- ut - eta * grad
    # Apply Regularisation:
    if (penalize_X){
      if(is.null(Gamma_u) == FALSE){
        new_u <- get_gamma_sparse_U(new_u, lambda_gamma_x)
      }else{
        if(is.null(Zu) == FALSE){
          new_u <- group_sparse(new_u, kx, r, Zu)
        }else{
          new_u <- hard(new_u, kx, r)
        }
      }
    }
    #new_u <- new_u %*% sqrtm(t(new_u) %*% Sigma_X %*% new_u)$Binv
    
    grad <- 1/n * (Sigma_Y %*% vt - t(Sigma_XY) %*% new_u)  + lambda_y  * Sigma_Y %*% vt %*% (t(vt) %*% Sigma_Y %*% vt- diag(r));
    new_v <- vt - eta * grad
    # Apply Regularisation:
    if (penalize_Y){
      if(is.null(Gamma_v) == FALSE){
        new_v <- get_gamma_sparse_U(new_v, lambda_gamma_y)
      }else{
        if(is.null(Zu) == FALSE){
          new_v <- group_sparse(new_v, ky, r, Zv)
        }else{
          new_v <- hard(new_v, ky, r)
        }
      }
    }
    #new_v <- new_v %*% sqrtm(t(new_v) %*% Sigma_Y %*% new_v)$Binv
    
    criteria <- mean(c( subdistance(ut, new_u),  subdistance(vt, new_v)))
    ut <- new_u
    vt <- new_v
    #print(criteria)
    iter <- iter+1
  }
  if (plot){
    plot(error[1:iter], type='l',  main="Matrix distance and iterations", 
         xlab="Number of iterations", ylab="Matrix distance",)
  }
  final_estimate_u <- ut %*% sqrtm(t(ut) %*% Sigma_X %*% ut)$Binv
  final_estimate_v <- vt %*% sqrtm(t(vt) %*% Sigma_Y %*% vt)$Binv
  return(list(Uhat =final_estimate_u,
              Vhat = final_estimate_v))
}



hard <-
  function(U, k, r){
    if(r>1){
      truncate.value <- sort(apply(U, 1, FUN = function(x) sum(x^2)),decreasing=TRUE)[k]
      U[which(apply(U, 1, FUN = function(x) sum(x^2))<truncate.value),] <- 0
    }else{
      truncate.value <- sort(abs(U),decreasing=TRUE)[k]
      U[which(abs(U)<truncate.value)] <- 0
    }
    return(U)
  }


compute_group_norm <- function(ut, Z){
  list_Z <- unique(Z)
  norm_u <-  apply(ut^2, 1, sum)
  norm_G <- sapply(1:length(list_Z), function(i){
    z = list_Z[i]
    ind = which(Z == z)
    return(sqrt(sum(norm_u[ind])))
  })
  return(norm_G)
}

group_hard <-
  function(U, k, r, Zu){
    list_Z <- unique(Zu)
      norms_u <- compute_group_norm(U, Zu)
      truncate.value <- sort(norms_u,decreasing=TRUE)[k]
      for (i in 1:length(list_Z)){
        ind = which(Zu == list_Z[i])
        if (norms_u[i] < truncate.value){
          U[ind,] <- 0
        }
      }
    return(U)
  }

get_gamma_sparse_U <- function(U, Gamma_u){
  U_opt <- Variable(nrow(U), ncol(U))
  constraints <- list(p >= 0, p <= 1, p[subject_0] == 1)
  # Define the quadratic loss
  loss <- sum((U_opt - U)^2)/ nrow(U) 
  # Define the L-1 norm term
  l21_norm <- cvxr_norm(Gamma %*% U_opt, p = 2, axis = 1)
  # Define the objective
  objective <- Minimize(loss + lambda * l21_norm)
  # Formulate the problem
  problem <- Problem(objective, constraints)
  # Solve the problem
  result_problem <- solve(problem)
  # Get the optimal value of p
  U_sol <- result_problem$getValue(U_opt)
  return(U_sol)
}



pipeline_CCA_gd2 <- function(X, Y, 
                             r=2,
                             param_lambda=c(1, 0.1, 0.01, 0.001),
                             param_k=c(5, 10, 20, 50, 1000), 
                             kfolds=10,
                             maxiter=5000, convergence=1e-4, 
                             eta=1e-3){
  p1 <- dim(X)[2]
  p2 <- dim(Y)[2]
  p <- p1 + p2;
  n <- nrow(X)
  pp <- c(p1,p2);
  Sigma_X = cov(X)
  Sigma_Y = cov(Y)

  init <- init_cca(X,Y, r, penalize_X=TRUE,
                   penalize_Y=TRUE,
                   Zu=NULL, A_u=NULL,
                   Zv=NULL, A_v=NULL)
  
  resultsx =c()
  for (lambda in param_lambda){
    
    resultsx_temp <- expand.grid(kx = param_k, ky = param_k) %>%
      mutate(rmse_str = map2_chr(kx, ky, ~ cv_function_gd(X, Y,  
                                                      kfolds=kfolds, init,
                                                      lambda_x = lambda,
                                                      lambda_y = lambda,
                                                      kx = .x,
                                                      ky = .y, 
                                                      r=r,
                                                      maxiter=maxiter, 
                                                      eta=eta, convergence=convergence)))%>%
      separate(rmse_str, 
               into = c("rmse", "sd"), sep = "_", convert = TRUE)
    resultsx_temp["lambda"] = lambda
    resultsx = rbind(resultsx, resultsx_temp)
  }

    #X, Y, Mask, kfolds=5, ainit,lambda, r=2, k=20,  
    #maxiter=1000, eta=0.001, convergence=1e-
    
    #print(resultsx)
    ###### (X, Y, Mask, kfolds=5, ainit, lambda, k=20)
    
    # print best hyperparameters and corresponding RMSE
  ind = which.min(resultsx$rmse)
  candidates = which( resultsx$rmse < resultsx$rmse[ind] + 0.5 * resultsx$sd[ind]) ### within 
  opt = which.min(apply(resultsx[candidates,c("kx", "ky")], 1, sum))
  opt_kx <- resultsx$kx[candidates[opt]]
  opt_ky <- resultsx$ky[candidates[opt]]
  opt_lambda <- resultsx$lambda[candidates[opt]]
  #which_lambdax = which(abs(resultsx$rmse-min(resultsx$rmse))/(1e-6  + min(resultsx$rmse)) <0.01)
  #lambda = max(resultsx$lambda[which_lambdax])
  #k = max(resultsx$k[which_lambdax])
  print(c("selected", opt_kx, opt_ky, opt_lambda))
  
  final <- cca_gd2(X, Y, r, init, kx=opt_kx, ky=opt_ky,
                   lambda_x = opt_lambda, lambda_y = opt_lambda, 
                   eta=eta, 
                   convergence=convergence,
                   penalize_X = TRUE, penalize_Y = TRUE,
                   lambda_gamma_x = 0,
                   lambda_gamma_y = 0,
                   maxiter=maxiter, plot = FALSE)
  Uhat = final$Uhat %*% sqrtm(t(final$Uhat) %*% Sigma_X %*% final$Uhat)$Binv 
  Vhat = final$Vhat %*% sqrtm(t(final$Vhat) %*% Sigma_Y %*% final$Vhat)$Binv 
  return(list( ufinal = Uhat, vfinal = Vhat,
               initu=init$u_init, 
               initv=init$v_init,
               lambda=opt_lambda, 
               kx=opt_kx,
               ky=opt_ky,
               resultsx=resultsx
  ))
  
}


cv_function_gd <- function(X, Y, kfolds=5, init,
                           lambda_x=0.01,
                           lambda_y =0.01,
                           r=2, kx=20, ky=20,
                           maxiter=2000, eta=0.001,
                           convergence=1e-4) {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- numeric(length = kfolds)
  p1 <- dim(X)[2]
  p2 <- dim(Y)[2]
  p <- p1 + p2;
  n <- nrow(X)
  
  #init <- gca_to_cca(ainit, S0, pp)
  # loop over folds
  for (i in seq_along(folds)) {
    # split data into training and validation sets
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]
    
    # fit model on training data with hyperparameters
    tryCatch(
      {
        final = cca_gd2(X_train, Y_train, r, 
                        init, kx=kx, ky=ky,
                        lambda_x = lambda_x, lambda_y = lambda_y, 
                        eta=eta, 
                        convergence=convergence,
                        penalize_X = TRUE, penalize_Y = TRUE,
                        lambda_gamma_x = 0,
                        lambda_gamma_y = 0,
                        maxiter=maxiter, plot = FALSE)
        # make predictions on validation data
        # compute RMSE on validation data
        #rmse[i] <- subdistance(X_val %*% final$Uhat, Y_val%*% final$Vhat)
        rmse[i] <- mean((X_val %*% final$Uhat -  Y_val%*% final$Vhat)^2)
        print(rmse)
      },
      error = function(e) {
        # If an error occurs, assign NA to the result
        rmse[i] <- NA
      })
  }
  
  # return mean RMSE across folds
  if (mean(is.na(rmse)) == 1){
    return(1e8)
  }else{
    return(paste0(mean(rmse, na.rm=TRUE), "_", sd(rmse, na.rm=TRUE)))
  }
}

