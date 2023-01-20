
library(RSpectra)
#### We have to turn this into a pipeline.

#### Question ---- do we expect this to be sparse?
#### Not necessarily --- there could be many brain regions that are activated. 

SVCST <- function(W, k, max_singval=5){
  #### g_i min of (1, (w_i - gamma)_+)
  #### minimze gamma such that sum g_i \leq r
  if(max_singval > min(dim(W))){
    svd_W = svd(W, max_singval)
  }else{
    svd_W = svds(W, max_singval)
  }
  lambdas = svd_W$d
  res <- sapply(lambdas[which(lambdas>0)], FUN=function(x){
    if(length(which(lambdas>x))>0){
      sum(sapply( lambdas[which(lambdas >x)] -x, FUN=function(y){min(1, y)}))}
    else{
      0
    }
  })
  indices = which(res<=k)
  index = max(indices)
  if(index > min(dim(W))){
    svd_W = svd(W, index)
  }else{
    svd_W = svds(W, index)
  }
  
  return(svd_W$u %*% diag(svd_W$d) %*% t(svd_W$v))
}

convex_formulation <- function(D = NULL,Sigma_x, Sigma_y, Sigma_xy, H,
                              G, lambda, eta = 0.5,
                              max_k = 10, max.iter=10,
                              ZERO_THRESHOLD= 1e-6, verbose = FALSE){
  
  svd_x = svds(Sigma_x, k = min(max_k, min(dim(Sigma_x))-1))
  sqrt_lambas_x = sqrt(svd_x$d)
  inv_sqrt_lambas_x = sapply(sqrt_lambas_x, function(x){ifelse(x ==0, 0, 1/x)})
  sqrt_Sigma_x = svd_x$u %*% diag(sqrt_lambas_x) %*% t(svd_x$v)
  inv_sqrt_Sigma_x = svd_x$u %*% diag(inv_sqrt_lambas_x) %*% t(svd_x$v)
  
  svd_y = svds(Sigma_y, k = min(max_k, min(dim(Sigma_y))-1))
  sqrt_lambas_y = sqrt(svd_y$d)
  inv_sqrt_lambas_y = sapply(sqrt_lambas_y, function(x){ifelse(x ==0, 0, 1/x)})
  sqrt_Sigma_y = svd_y$u %*% diag(sqrt_lambas_y) %*% t(svd_y$v)
  inv_sqrt_Sigma_y = svd_y$u %*% diag(inv_sqrt_lambas_y) %*% t(svd_y$v)
  
  
  #Gamma = kronecker(sqrt_Sigma_y, sqrt_Sigma_x, FUN = "*")
  #b <- as.vector(G - 1/eta *H + 1/eta * inv_sqrt_Sigma_x %*% Sigma_xy %*% inv_sqrt_Sigma_y)
  #N <- length(b)
  #p <- ncol(Gamma)
  p <- ncol(sqrt_Sigma_x)
  m <- nrow(sqrt_Sigma_y)
  b = G - 1/eta *H + 1/eta * inv_sqrt_Sigma_x %*% Sigma_xy %*% inv_sqrt_Sigma_y
  ## Define_Parameters
  Fhat <- Variable(p, m)
  ## Build_Penalty_Terms
  #penalty_term1 <- sum(cvxr_norm(hstack(beta, theta), 2, axis = 1))
  #penalty_term2 <- sum(cvxr_norm(theta, 2, axis = 1))
  penalty_term1 <- sum(cvxr_norm(D %*% Fhat, p=1, axis=2))
  
  ## Compute_Fitted_Value
  y_hat <- sqrt_Sigma_x %*% Fhat %*% sqrt_Sigma_y
  ## Build_Objective
  objective <- eta / 2 * sum_squares(b - y_hat) + lambda  * penalty_term1
  ## Define_and_Solve_Problem
  prob <- Problem(Minimize(objective))
  result <- solve(prob, verbose = TRUE, num_iters = max.iter)
  ## Return_Values
  Fhat <- result$getValue(Fhat)
  
  ## Zero out stuff before returning
  Fhat[abs(Fhat) < ZERO_THRESHOLD] <- 0.0
  y_hat <- sqrt_Sigma_x %*% Fhat %*% sqrt_Sigma_y
  list(
    Fhat = Fhat,
    yhat = y_hat,
    criterion = result$value)
}



gen_lasso <- function(D, Sigma_x, Sigma_xy, V,
                        lambda, lambda2, max.iter=5000,
                        max_k = 10, verbose = FALSE,
                        ZERO_THRESHOLD=1e-5){
  svd_x = svds(Sigma_x, k = min(max_k, min(dim(Sigma_x))-1))
  sqrt_lambas_x = sqrt(svd_x$d)
  inv_sqrt_lambas_x = sapply(sqrt_lambas_x, function(x){ifelse(x ==0, 0, 1/x)})
  sqrt_Sigma_x = svd_x$u %*% diag(sqrt_lambas_x) %*% t(svd_x$v)
  inv_sqrt_Sigma_x = svd_x$u %*% diag(inv_sqrt_lambas_x) %*% t(svd_x$v)
  
  p <- nrow(sqrt_Sigma_x)
  ## Define_Parameters
  Uhat <- Variable(p, k)
  ## Build_Penalty_Terms
  penalty_term1 <- sum(cvxr_norm(D %*% Uhat, 1, axis = 2))
  penalty_term2 <- sum((D %*% Uhat)^2)
  
  ## Compute_Fitted_Value
  y_hat <- sqrt_Sigma_x %*% Uhat 
  b <- inv_sqrt_Sigma_x %*% Sigma_xy %*% V
  ## Build_Objective
  objective <- 1 / 2 * sum_squares(y_hat - b) + lambda  * penalty_term1 + lambda2  * penalty_term2
  ## Define_and_Solve_Problem
  prob <- Problem(Minimize(objective))
  result <- solve(prob, verbose = TRUE, max.iter=max.iter)
  ## Return_Values
  Uhat <- result$getValue(Uhat)
  
  ## Zero out stuff before returning
  Uhat[abs(Uhat) < ZERO_THRESHOLD] <- 0.0
  list(
    Uhat = Uhat,
    criterion = result$value)
}


smooth_lasso <- function(D, Sigma_x, Sigma_xy, V,
                      lambda, lambda2, max.iter=5000,
                      max_k = 10, verbose = FALSE,
                      ZERO_THRESHOLD=1e-5){
  svd_x = svds(Sigma_x, k = min(max_k, min(dim(Sigma_x))-1))
  sqrt_lambas_x = sqrt(svd_x$d)
  inv_sqrt_lambas_x = sapply(sqrt_lambas_x, function(x){ifelse(x ==0, 0, 1/x)})
  sqrt_Sigma_x = svd_x$u %*% diag(sqrt_lambas_x) %*% t(svd_x$v)
  inv_sqrt_Sigma_x = svd_x$u %*% diag(inv_sqrt_lambas_x) %*% t(svd_x$v)
  
  p <- nrow(sqrt_Sigma_x)
  ## Define_Parameters
  Uhat <- Variable(p, k)
  ## Build_Penalty_Terms
  penalty_term1 <- sum(cvxr_norm(Uhat, 1, axis = 2))
  if (is.null(D) ==FALSE){
    penalty_term2 <- sum(cvxr_norm(D %*% Uhat, p=1, axis = 2))
  }

  
  ## Compute_Fitted_Value
  y_hat <- sqrt_Sigma_x %*% Uhat 
  b <- inv_sqrt_Sigma_x %*% Sigma_xy %*% V
  ## Build_Objective
  if (is.null(D)==FALSE){
  objective <- 1 / 2 * sum_squares(y_hat - b) + lambda  * penalty_term1 + lambda2  * penalty_term2
  }else{
    objective <- 1 / 2 * sum_squares(y_hat - b) + lambda  * penalty_term1 
  } 
  ## Define_and_Solve_Problem
  prob <- Problem(Minimize(objective))
  result <- solve(prob, verbose = TRUE, max.iter=max.iter)
  ## Return_Values
  Uhat <- result$getValue(Uhat)
  
  ## Zero out stuff before returning
  Uhat[abs(Uhat) < ZERO_THRESHOLD] <- 0.0
  list(
    Uhat = Uhat,
    criterion = result$value)
}


sparse_lasso <- function(D, Sigma_x, Sigma_xy, V,
                         lambda, lambda2, max.iter=5000,
                         max_k = 10, verbose = FALSE,
                         ZERO_THRESHOLD=1e-5){
  svd_x = svds(Sigma_x, k = min(max_k, min(dim(Sigma_x))-1))
  sqrt_lambas_x = sqrt(svd_x$d)
  inv_sqrt_lambas_x = sapply(sqrt_lambas_x, function(x){ifelse(x ==0, 0, 1/x)})
  sqrt_Sigma_x = svd_x$u %*% diag(sqrt_lambas_x) %*% t(svd_x$v)
  inv_sqrt_Sigma_x = svd_x$u %*% diag(inv_sqrt_lambas_x) %*% t(svd_x$v)
  
  p <- nrow(sqrt_Sigma_x)
  ## Define_Parameters
  Uhat <- Variable(p, k)
  ## Build_Penalty_Terms
  penalty_term1 <- sum(cvxr_norm(Uhat, 1, axis = 1))
  if (is.null(D) ==FALSE){
    penalty_term2 <- sum(cvxr_norm(D %*% Uhat,1, axis = 1))
  }
  
  
  ## Compute_Fitted_Value
  y_hat <- sqrt_Sigma_x %*% Uhat 
  b <- inv_sqrt_Sigma_x %*% Sigma_xy %*% V
  ## Build_Objective
  if (is.null(D)==FALSE){
    objective <- 1 / 2 * sum_squares(y_hat - b) + lambda  * penalty_term1 + lambda2  * penalty_term2
  }else{
    objective <- 1 / 2 * sum_squares(y_hat - b) + lambda  * penalty_term1 
  } 
  ## Define_and_Solve_Problem
  prob <- Problem(Minimize(objective))
  result <- solve(prob, verbose = TRUE, max.iter=max.iter)
  ## Return_Values
  Uhat <- result$getValue(Uhat)
  
  ## Zero out stuff before returning
  Uhat[abs(Uhat) < ZERO_THRESHOLD] <- 0.0
  list(
    Uhat = Uhat,
    criterion = result$value)
}

genCCA_Chao<- function(D, X, Y, k, lambda1, lambda2, lambda3,  
                       penalty_type = "GEN",
                       eta=1, 
                       max_k=30, zero_threshold = 1e-6, 
                       epsilon=0.01,
                       max.iter = 5000,
                       verbose=TRUE){
  Sigma_x = t(X) %*% X
  
  Sigma_y = t(Y) %*% Y
  Sigma_xy = t(X) %*% Y
  #### Initialize 
  old.F = SVCST(Sigma_xy, k=max_k)
  old.G = matrix(0, nrow=nrow(Sigma_xy), ncol = ncol(Sigma_xy))
  old.H = matrix(0, nrow=nrow(Sigma_xy), ncol = ncol(Sigma_xy))
  converged = FALSE
  it <- 0
  while(converged == FALSE){
    lasso_sol <- convex_formulation(D, Sigma_x, Sigma_y, Sigma_xy, old.H,
                                   old.G, lambda1, eta = eta,
                                   max_k = max_k, max.iter=max.iter,
                                   ZERO_THRESHOLD= zero_threshold, 
                                   verbose = verbose)
    new.F = lasso_sol$Fhat
    y_hat = lasso_sol$yhat
    
    new.G <-SVCST(1/eta * old.H + y_hat, k=k)
    new.H <- old.H + eta * (y_hat - new.G )
    it <- it + 1
    converged = ((max(CVXR::norm(new.F - old.F, 'F')/CVXR::norm(old.F, 'F'), CVXR::norm(new.G - old.G, 'F')/CVXR::norm(old.G, 'F')) <= epsilon) | it > max.iter)
    print(c(it, CVXR::norm(new.F - old.F, 'F')/CVXR::norm(old.F), CVXR::norm(new.G - old.G, 'F')/CVXR::norm(old.G, 'F')))
    old.F = new.F
    old.G = new.G
    old.H = new.H
  }
  
  
  A = new.F
  if(k >= min(dim(A))){
    svd_A  = svd(A, k)
  }else{
    svd_A  = svds(A, k)
  }
  
  if (penalty_type =="GEN"){
    res =  gen_lasso(D, Sigma_x, Sigma_xy, svd_A$v,
                     lambdaA1, lambdaA2, max.iter=5000,
                     max_k = 10, verbose = FALSE,
                     ZERO_THRESHOLD=1e-5)
  }else{
    res =  smooth_lasso(D, Sigma_x, Sigma_xy, svd_A$v,
                     lambdaA1, lambdaA2, max.iter=5000,
                     max_k = 10, verbose = FALSE,
                     ZERO_THRESHOLD=1e-5)
  }
  
  res_y = smooth_lasso(Sigma_y, t(Sigma_xy), svd_A$u,
                    lambdaB1, 0, max.iter=5000,
                    max_k = 10, verbose = FALSE,
                    ZERO_THRESHOLD=1e-5)
  
  
  return(list(
    xcoef = res$Uhat,
    ycoef = res$Vhat))
}


sparseCCA.CV<-function(X, Y, D, rank, n.cv=5,
                       lambda1seq=matrix(seq(from=0,to=10,by=0.1),nrow=1),
                       lambda2seq=matrix(seq(from=0,to=10,by=0.1),nrow=1),
                       lambda3seq=matrix(seq(from=0,to=10,by=0.1),nrow=1),
                       eta=1, 
                       max_k=30, zero_threshold = 1e-6, 
                       epsilon=0.1,
                       max.iter = 15,
                       verbose=TRUE){ 
  # Code to  maximizing test sample correlation with genCCA approach
  
  ### INPUT
  # X                   : (nxp) data matrix
  # Y                   : (nxq) data matrix
  # n.cv                : n.cv-fold cross-validation
  # lambdax             : grid of sparsity parameters for lambdax
  # lambday             : grid of sparsity parameters for lambday
  
  ### OUTPUT
  # lambdax.opt         : value of the sparsity parameter lambdax
  # lambday.opt         : value of the sparsity parameter lambday
  
  ### START CODE
  
  # Dimensions
  n = nrow(X)
  p = ncol(X)
  n.cv.sample<-trunc(n/n.cv)
  whole.sample<-seq(1,n)
  
  train = sample(1:n, n)
  cv.sets <- split(train, ceiling(seq_along(1:n)/(n.cv.sample)))
  print("here")
  
  cv.results <- do.call(rbind, lapply(1:n.cv, function(i){
    testing.sample<- train[c(cv.sets[[i]])]
    training.sample<-train[-c(cv.sets[[i]])]
    Xcv = X[training.sample, ]
    Ycv = Y[training.sample, ]
    Xtest=X[testing.sample,]
    Ytest=Y[testing.sample,]
    return(as.data.frame(do.call(rbind, 
                                 #mclapply(lambda1seq, genCCA.cv.lambda1_2, 
                                 lapply(lambda1seq, sparseCCA.cv.lambda1_2, 
                                        Xtrain=Xcv,Ytrain=Ycv, D=D, Xtest=Xtest,
                                        Ytest=Ytest,
                                        lambday=lambda2seq,
                                        eta=eta, max_k=max_k, zero_threshold = zero_threshold, 
                                        epsilon=epsilon,
                                        max.iter=max.iter,
                                        verbose=verbose)))) 
  }))
  
  cv.results_m = cv.results %>% group_by(lambda1, lambda2) %>% summarize(m=mean(score))
  indx= which.min(cv.results$m)
  Fit.sparseCCA <- sparseCCA_Chao(X=Xtrain, Y=Ytrain, rank, cv.results$lambda1[indx], 
                             cv.results$lambda2[indx],
                             eta=eta, max_k=max_k, zero_threshold = zero_threshold, 
                             epsilon=epsilon,
                             max.iter=max.iter,
                             verbose=verbose)
  # for (i in 1:n.cv){
  #   testing.sample<- train[c(cv.sets[[i]])]
  #   training.sample<-train[-c(cv.sets[[i]])]
  #   Xcv = X[training.sample, ]
  #   Ycv = Y[training.sample, ]
  #   Xtest=X[testing.sample,]
  #   Ytest=Y[testing.sample,]
  #   #profvis({
  #   cvscore[,,i]<-matrix(mclapply(lambda1seq, genCCA.cv.lambda1_2, 
  #                        Xtrain=Xcv,Ytrain=Ycv, D=D, Xtest=Xtest,
  #                        Ytest=Ytest,
  #                        lambday=lambda2seq), ncol=3, byrow=length(lambda1seq))
  # }
  
  # toc()
  
  lambdax.opt<-res_m$lambda1[indx]
  lambday.opt<-res_m$lambda2[indx]
  
  ### OUTPUT
  out<-list(fit.gcc = Fit.sparseCCA,
            lambdax.opt=lambdax.opt,lambday.opt=lambday.opt)
}



sparseCCA.cv.lambda1_2<-function(U,Xtrain,Ytrain, Xtest,Ytest,lambday,
                                 rank=3,
                                 eta=1,
                                 max_k=30,
                                 zero_threshold = 1e-6,
                                 epsilon=0.1,
                                 max.iter = 15){ #AUXILIARY FUNCTION
  testcorrelations<-c()
  A.init = A.initial
  B.init = B.initial
  for (V in lambday){
    print(c(V, "lambda1_2"))
    res <- sparseCCA.cv.lambda2_2(V, Xtrain,Ytrain, D, 
                                  Xtest,Ytest, U,
                                  rank=rank,
                                  eta=eta,
                                  max_k=max_k,
                                  zero_threshold = zero_threshold,
                                  epsilon=epsilon,
                                  max.iter = max.iter)
    A.init <- res$xcoef
    B.init <- res$ycoef
    testcorrelations<-c(testcorrelations,
                        res$score)
  }
  return(data.frame("score" = testcorrelations,
                    "lambda1" = U,
                    "lambda2" =lambday ))
}


sparseCCA.cv.lambda2_2<-function(V,Xtrain,Ytrain, 
                                 Xtest,Ytest,lambdaxfixed,
                                 rank=3, eta=1,
                                 max_k=30,
                                 zero_threshold = 1e-6,
                                 epsilon=0.1,
                                 max.iter = 15,
                                 verbose=TRUE,
                                 A.initial=NULL){ #AUXILIARY FUNCTION
  #X, Y, Da, Db, lambdaA1=NULL, lambdaB1=NULL,lambdaA2=1.,lambdaB2=1.,rank,A.initial=NULL,B.initial=NULL,max.iter=20,conv=10^-2,mode = c("glmnet", "ECOS", "CD")
  print(c(V,lambdaxfixed ))
  Fit.genCCA <- sparseCCA_Chao(Xtrain, Ytrain, rank, V, lambdaxfixed, 
                          eta=eta,max_k=max_k, 
                          zero_threshold = zero_threshold, 
                          epsilon=epsilon,
                          max.iter = max.iter,
                          verbose=verbose)
  return(list("score" =  mean((Xtest%*%Fit.genCCA$xcoef - Ytest%*%Fit.genCCA$ycoef)^2),#trabs(cor(Xtest%*%Fit.genCCA$xcoef,Ytest%*%Fit.genCCA$ycoef)),
              "prediction error" = mean((Xtest%*%Fit.genCCA$xcoef - Ytest%*%Fit.genCCA$ycoef)^2),
              "xcoef"= Fit.genCCA$xcoef,
              "ycoef"= Fit.genCCA$ycoef))
}  


