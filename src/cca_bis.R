# LIBRARIES
library(glmnet)
library(MASS)
library(pracma)
library(CVXR)

source("src/cd_solver.R")
source("src/penalties.R")
source("src/losses.R")
source("experiments/alternative_methods/SAR.R")
source("src/genglm.R")

genCCA2<-function(X, Y,
                 Da, Db,
                 lambdaA1=NULL,
                 lambdaB1=NULL,
                 lambdaA2=1.,
                 lambdaB2=1.,
                 rank,
                 A.initial=NULL,B.initial=NULL,
                 max.iter=20,conv=10^-2,
                 solver = c("glmnet", "ECOS", "CGD"),
                 standardize=TRUE,
                 verbose=FALSE,
                 reexpress=FALSE){
  ### Function to perform Sparse Canonical Correlation Analysis using alternating regressions

  ### INPUT
  # X                   : (nxp) data matrix
  # Y                   : (nxq) data matrix
  # lambdaA1          : grid of sparsity parameters for lambdaA
  # lambdaB1          : grid of sparsity parameters for lambdaB
  # lambdaA2          : grid of sparsity parameters for lambda2 for A
  # lambdaB2          : grid of sparsity parameters for lambda2 for B
  # rank                : number of canonical vector pairs to extract
  # selection.criterion : 1 for BIC, 2 for Cross-validation to select sparsity parameters
  # n.cv                : n.cv-fold cross-validation
  # A.initial           : starting value for the canonical vector A
  # B.initial           : starting value for the canonical vector B
  # max.iter            : maximum number of iterations
  # conv                : tolerance value for convergence

  ### OUTPUT
  # ALPHA               : (pxr) estimated canonical vectors correpsonding to the first data set
  # BETA                : (qxr) estimated canonical vectors correpsonding to the second data set
  # cancors             : r estimated canonical correlations
  # U_ALL               : (nxr) estimated canonical variates corresponding to the first data set
  # V_ALL               : (nxr) estimated canonical variates corresponding to the second data set
  # lamdbaA             : value of the sparsity parameter lambdaA
  # lamdbaB             : value of the sparsity parameter lambdaB
  # it                  : number of iterations


  ### STORE RESULTS
  ##### We'll add CV later
  df.x <- data.frame(X) %>% mutate_all(~(scale(.) %>% as.vector))

  df.y <- data.frame(Y) %>% mutate_all(~(scale(.) %>% as.vector))
  print("here")
  print(rank)

  X = as.matrix(df.x)
  Y = as.matrix(df.y)
  print(c(ncol(X), ncol(Y), rank, max.iter))
  ALPHA_ALL<-matrix(NA,ncol=rank,nrow=ncol(X))
  BETA_ALL<-matrix(NA,ncol=rank,nrow=ncol(Y))
  U_ALL<-matrix(NA,ncol=rank,nrow=nrow(X))
  V_ALL<-matrix(NA,ncol=rank,nrow=nrow(Y))
  cancors<-matrix(NA,ncol=rank,nrow=1)

  ALPHA_trans<-matrix(NA,ncol=rank,nrow=ncol(X))
  BETA_trans<-matrix(NA,ncol=rank,nrow=ncol(Y))
  U_trans<-matrix(NA,ncol=rank,nrow=nrow(X))
  V_trans<-matrix(NA,ncol=rank,nrow=nrow(Y))
  cancors_trans<-matrix(NA,ncol=rank,nrow=1)

  lambdaA_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  lambdaB_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  iterations<-matrix(NA,ncol=rank,nrow=1)
  obj.it<-matrix(NA,ncol=rank,nrow=max.iter+1)

  p = dim(X)[2]
  n = dim(Y)[1]
  q = dim(Y)[2]


  ### START CODE
  # Starting Values: Canonical Ridge Solution
  if(is.null(A.initial)){
    #cancor_regpar<-estim.regul_crossvalidation(X,Y,n.cv=n.cv)
    #cancor_regul<-rcc(X,Y,cancor_regpar$lambda1.optim,cancor_regpar$lambda2.optim)
    A.initial_ridge<- matrix(rnorm(rank * ncol(X)),nrow=ncol(X),ncol=rank)
    #cancor_regul <- cancor(X,Y)
    #A.initial_ridge<-matrix(cancor_regul$xcoef[,1:rank],ncol=rank,nrow=ncol(X))
    #B.initial_ridge<-matrix(cancor_regul$ycoef[,1:rank],ncol=rank,nrow=ncol(Y))
    A.initial<-apply(A.initial_ridge,2, NORMALIZATION_UNIT)
  }
  
  B.initial_ridge<- matrix(rnorm(rank *ncol(Y)),nrow=ncol(Y),ncol=rank)
  B.initial<-apply(B.initial_ridge,2,NORMALIZATION_UNIT)

  # Sequentially extract the r canonical vector pairs
  for (i.r in 1:rank){
    
    if (verbose) {print(paste0("Dimension is " ,i.r))}
    if (i.r==1){ # for r=1: start from original data sets
      X_data<-X
      Y_data<-Y
    }

    # STARTING VALUES: canonical ridge
    A.STARTING<-matrix(A.initial[,i.r],ncol=1)
    B.STARTING<-matrix(B.initial[,i.r],ncol=1)
  
    # CONVERGENCE CRITERION
    obj.initial<-mean((X_data%*%matrix(A.STARTING,ncol=1)-Y_data%*%matrix(B.STARTING,ncol=1))^2)
    obj.it[1,i.r]<-obj.initial
    # INITIALIZE CONVERGENCE PARAMETERS
    it<-1
    diff.obj<-conv * 10
  

    # FROM i until convergence: canonical vectors
    print(paste0("max.iter is", max.iter))
    while( (it<max.iter) & (diff.obj>conv)  ){

      # Estimating A conditional on B
      if(verbose){print(it)}
      print(it)
      
      FIT.A <- alternating.optimization(X=X_data, y=Y_data%*%B.STARTING,
                                        D=Da,
                                        lambda1=lambdaA1,
                                        lambda2=lambdaA2,
                                        solver=solver,
                                        eps = conv,
                                        max_it = 30)

      AHAT_FINAL <- FIT.A/ (1e-8 +  norm(X_data %*% FIT.A, type="F"))
      A.STARTING <-  AHAT_FINAL
      # Estimating B conditional on A

      # Estimating A conditional on B
      FIT.B <- alternating.optimization(X=Y_data, y=X_data%*%A.STARTING,
                                        D = Db,
                                        lambda1=lambdaB1,
                                        lambda2=lambdaB2,
                                        solver=solver,
                                        eps = conv,
                                        max_it = 30)

      BHAT_FINAL<- FIT.B/ (1e-8 + norm(Y_data %*% FIT.B, type="F"))
      B.STARTING <-  BHAT_FINAL
      if(norm(BHAT_FINAL, type='F') ==0){
        print("Regularization too strong!!!")
        break();
      }
      print("here")
      #print(AHAT_FINAL)

      #lambdaB_ALL[it,i.r]<-FIT.B$lambda2

      # Check convergence
      obj.new<-mean((X_data%*%matrix(AHAT_FINAL,ncol=1)-Y_data%*%matrix(BHAT_FINAL,ncol=1))^2)
      obj.it[it+1,i.r]<-obj.new
      diff.obj<-abs(obj.new-obj.initial)/obj.initial
      if (verbose) { 
        print(c(it, diff.obj))}

      # Updated starting values
      #B.STARTING<-BHAT_FINAL
      #A.STARTING<-AHAT_FINAL
      obj.initial<-obj.new
      it<-it+1
      # if((mean(BHAT_FINAL) == BHAT_FINAL[1]) || (mean(AHAT_FINAL) == AHAT_FINAL[1])){
      #   it = max.iter + 1
      # }
    } # end while-loop

    # Number of ITERATIONS
    iterations[1,i.r]<-it

    # CANONICAL VARIATES after convergence
    Uhat<-X_data%*%AHAT_FINAL
    Vhat<-Y_data%*%BHAT_FINAL
    
     # Express canonical vectors in terms of ORIGINAL DATA MATRICES
    if (i.r==1){# FIRST DIMENSION
      # Final estimates of canonical vectors,  variates and canonical correlation
      ALPHA_ALL[, i.r] <- AHAT_FINAL
      BETA_ALL[, i.r] <- BHAT_FINAL
      U_ALL[, i.r] <- Uhat
      V_ALL[, i.r] <- Vhat
      cancors[1, i.r] <- cor(Uhat, Vhat)

      # Deflated data matrices
      X_data <-  round(X_data  - Uhat%*%solve(t(Uhat)%*%Uhat)%*%t(Uhat)%*%X_data, 10)
      Y_data <-  round(Y_data - Vhat%*%solve(t(Vhat)%*%Vhat)%*%t(Vhat)%*%Y_data, 10)
      print(i.r)



    } else {# HIGHER ORDER DIMENSIONS
       print("higher d")

      # # A expressed in terms of original data set X
      if (reexpress){
        FIT.Aorig <- alternating.optimization(X=X, y=Uhat,
                                          D=Da,
                                          lambda1=lambdaA1,
                                          lambda2=lambdaA2,
                                          solver=solver,
                                          eps = conv,
                                          max_it = max.iter)
        ALPHAhat <- FIT.Aorig
        # ALPHAhat <- AHAT_FINAL

        # B expressed in terms of original data set Y
        FIT.Borig <- alternating.optimization(X=Y, y=Vhat,
                                          D=Db,
                                          lambda1=lambdaB1,
                                          lambda2=lambdaB2,
                                          solver=solver,
                                          eps = conv,
                                          max_it = max.iter)
        BETAhat <- FIT.Borig
        #BETAhat <- BHAT_FINAL
      }else{
        ALPHAhat <- AHAT_FINAL
        BETAhat <- BHAT_FINAL
      }

      # Final estimates of canonical vectors,  variates and canonical correlation
      ALPHA_ALL[, i.r] <- ALPHAhat
      BETA_ALL[, i.r] <- BETAhat
      Uhat <- X%*%ALPHAhat
      Vhat <- Y%*%BETAhat
      U_ALL[, i.r] <- Uhat
      V_ALL[, i.r] <- Vhat
      cancors[1, i.r] <- cor(Uhat, Vhat)

      # Deflated data matrices: regress original data sets on all previously found canonical variates
      X_data <-  X  - U_ALL[, 1:i.r]%*%solve(t(U_ALL[, 1:i.r])%*%U_ALL[, 1:i.r])%*%t(U_ALL[, 1:i.r])%*%X
      Y_data <-  Y -  V_ALL[, 1:i.r]%*%solve(t(V_ALL[, 1:i.r])%*%V_ALL[, 1:i.r])%*%t(V_ALL[, 1:i.r])%*%Y
      
    }
    # ALPHA_ALL[,i.r]<-AHAT_FINAL
    # BETA_ALL[,i.r]<-BHAT_FINAL
    # U_ALL[,i.r]<-Uhat
    # V_ALL[,i.r]<-Vhat
    # if (norm(Vhat, 'F')== 0 |norm(Uhat, 'F')== 0 ){
    #   cancors[1,i.r]<- 0
    # }else{
    #   cancors[1,i.r]<-cor(Uhat,Vhat)
    # }
    # if (norm(Uhat, 'F')>0){
    # X_data <-  X_data  - Uhat%*%solve(t(Uhat)%*%Uhat)%*%t(Uhat)%*%X_data
    # }else{
    #   X_data <-  X_data 
    # }
    # if (norm(Vhat, 'F')>0){
    #   Y_data <-  Y_data - Vhat%*%solve(t(Vhat)%*%Vhat)%*%t(Vhat)%*%Y_data 
    # }else{
    #   X_data <-  X_data 
    # }
    
  }
  out<-list(xcoef=ALPHA_ALL,ycoef=BETA_ALL, 
            cancors=cancors,
            U_ALL=U_ALL,V_ALL=V_ALL,
            it=iterations)
}



alternating.optimization <- function(X, y, D, lambda1, lambda2,
                                     solver="CGD", 
                                     eps = 1e-2,
                                     max_it = 200){
  p = dim(X)[2]
  n = dim(X)[1]
  if (is.null(lambda1)){
     lambda1 = 0
  }
  if (is.null(lambda2)){
     lambda2 = 0
  }

  if(is.null(D)){
      glm.model = glmnet(X, y, lambda = lambda1,
                         intercept = FALSE)
      FIT.A <- as.matrix(glm.model$beta)
  }
  else{
      res = genglm(X, y, D, family = "gaussian",
                  lambda1,
                  lambda2=lambda2,
                  standardize = FALSE,
                  solver = solver,
                  intercept = FALSE,
                  eps=eps,
                  max.it = max_it)
      FIT.A <- res$xcoef
  }

  return(FIT.A)
}

genCCA.CV<-function(X, Y, D, rank, n.cv=5,
                    lambda1seq=matrix(seq(from=0,to=10,by=0.1),nrow=1),
                    lambda2seq=matrix(seq(from=0,to=10,by=0.1),nrow=1),
                    lambdaB1seq=matrix(seq(from=0,to=10,by=0.1),nrow=1),
                    lambdaB2seq=NULL,
                    A.initial=NULL,B.initial=NULL,
                    max.iter=20,conv=10^-2,
                    solver ="CGD",
                    standardize=TRUE,
                    verbose=FALSE
                    ){ 
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
  m = ncol(Y)
  n.cv.sample<-trunc(n/n.cv)
  whole.sample<-seq(1,n)
  
  lambda1seq <- matrix(lambda1seq * 32 * sqrt(log(p)/n), ncol=1)
  lambda2seq <- matrix(lambda2seq * 32 * sqrt(log(p)/n), nrow=1)
  
  train = sample(1:n, n)
  cv.sets <- split(train, ceiling(seq_along(1:n)/(n.cv.sample)))
  print("here")
  
  cv.results <- do.call(rbind, mclapply(1:n.cv, function(i){
    testing.sample<- train[c(cv.sets[[i]])]
    training.sample<-train[-c(cv.sets[[i]])]
    Xcv = X[training.sample, ]
    Ycv = Y[training.sample, ]
    Xval=X[testing.sample,]
    Yval=Y[testing.sample,]
    return(as.data.frame(do.call(rbind, 
                                 mclapply(lambda1seq, genCCA.cv.lambda1_2, 
                                 #lapply(lambda1seq, genCCA.cv.lambda1_2, 
                  Xtrain=Xcv,Ytrain=Ycv, D=D, Xtest=Xval,
                  Ytest=Yval,
                  lambday=lambda2seq,
                  lambdaB1seq = lambdaB1seq,
                  lambdaB2seq = lambdaB2seq,
                  A.initial=A.initial, B.initial=B.initial,
                  max.iter=max.iter, 
                  conv=conv, 
                  solver= solver)))) 
    }))
  
  #print(cv.results)
  cv.results_m = cv.results %>% 
    group_by(lambda1, lambda2, lambdaB1, lambdaB2) %>% 
    summarize(m=mean(score)) %>%
    ungroup()
  indx= which.min(cv.results_m$m)
  print(cv.results_m)
  print(paste0("Selected index: ", indx  , " - ",  cv.results_m$lambda1[indx], " - ", cv.results_m$lambda2[indx]))
  Fit.genCCA <- genCCA2(X=X,Y=Y, Da = D, Db=NULL, 
                        A.initial = NULL,
                        B.initial = NULL,
                        rank=rank, 
                        lambdaA1=cv.results_m$lambda1[indx], 
                        lambdaA2 = cv.results_m$lambda2[indx],
                        lambdaB1 = cv.results_m$lambdaB1[indx],
                        lambdaB2 = cv.results_m$lambdaB2[indx],
                        max.iter=max.iter, 
                        conv=conv, 
                        solver= solver)

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
  
  lambdax.opt<-cv.results_m$lambda1[indx]
  lambday.opt<-cv.results_m$lambda2[indx]
  
  ### OUTPUT
  out<-list(fit.gcc = Fit.genCCA,
            lambdax.opt=lambdax.opt,lambday.opt=lambday.opt,
            cv.results = cv.results)
}



genCCA.cv.lambda1_2<-function(U,Xtrain,Ytrain,D, Xtest,Ytest,lambday,
                              lambdaB1seq = c(0),
                              lambdaB2seq = c(0),
                              rank=3,
                              A.initial=NULL, B.initial=NULL,
                              max.iter=20,conv=10^(-2),
                              solver = "CGD",
                              verbose=FALSE){ #AUXILIARY FUNCTION
  testcorrelations<-c()
  A.init = A.initial
  B.init = B.initial
  for (V in lambday){
    print(c(V, "lambda1_2"))
    res <- genCCA.cv.lambda2_2(V, Xtrain,Ytrain, D, 
                               Xtest,Ytest, U,
                               rank=rank,
                               lambdaB1 = lambdaB1seq,
                               lambdaB2 = lambdaB2seq,
                               A.initial=A.init,
                               B.initial=B.init)
    A.init <- res$xcoef
    B.init <- res$ycoef
    testcorrelations<-c(testcorrelations,
                        res$score)
  }
  return(data.frame("score" = testcorrelations,
              "lambda1" = U,
              "lambda2" = as.numeric(lambday) ))
}


genCCA.cv.lambda2_2<-function(V,Xtrain,Ytrain, D, Xtest,Ytest,lambdaxfixed,
                              rank=3,
                              A.initial=NULL,
                              B.initial=NULL,
                              max.iter=20,conv=1e-2,
                              solver = "CGD",
                              lambdaB1seq = c(0),
                              lambdaB2seq = c(0),
                              verbose=FALSE){ #AUXILIARY FUNCTION
  #X, Y, Da, Db, lambdaA1=NULL, lambdaB1=NULL,lambdaA2=1.,lambdaB2=1.,rank,A.initial=NULL,B.initial=NULL,max.iter=20,conv=10^-2,mode = c("glmnet", "ECOS", "CD")
  print(c(V,lambdaxfixed ))
  ##### mc-apply over the set of differnt combinations
  combs = merge(data.frame("lambdaB1"= lambdaB1seq), data.frame("lambdaB2"= lambdaB2seq),by=NULL)
  res <- as.data.frame(do.call(rbind, mclapply(1:nrow(combs), FUN=function(i){
    Fit.genCCA <- genCCA2(X=Xtrain,Y=Ytrain, Da = D, Db=NULL, 
                          A.initial = A.initial,
                          B.initial = B.initial,
                          rank=rank, lambdaA1=lambdaxfixed, 
                          lambdaA2 = V, 
                          lambdaB1 = combs$lambdaB1[i],
                          lambdaB2 = combs$lambdaB2[i],
                          max.iter=max.iter, 
                          conv=conv, 
                          solver= solver)
    return(list("score" = mean((Xtest%*%Fit.genCCA$xcoef - Ytest%*%Fit.genCCA$ycoef)^2), #trace(abs(cor(Xtest%*%Fit.genCCA$xcoef,Ytest%*%Fit.genCCA$ycoef))),
                "prediction error" = mean((Xtest%*%Fit.genCCA$xcoef - Ytest%*%Fit.genCCA$ycoef)^2),
                #"xcoef"= Fit.genCCA$xcoef,
                #"ycoef"= Fit.genCCA$ycoef,
                "lambdaA1" = lambdaxfixed,
                "lambdaA2" = V,
                "lambdaB1" = combs$lambdaB1[i],
                "lambdaB2" = combs$lambdaB2[i]))
  })))
  return(res)
}  


