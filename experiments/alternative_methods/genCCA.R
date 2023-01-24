# LIBRARIES
library(glmnet)
library(MASS)
library(pracma)
#library(ParallelLogger)
#library(rsample)


source("/Users/cdonnat/Documents/group-CCA/src/glmgen.R")
source("/Users/cdonnat/Documents/group-CCA/experiments/alternative_methods/SAR.R")
genCCA<-function(X, Y,
                 Da, Db,
                 lambdaA1=NULL,
                 lambdaB1=NULL,
                 lambdaA2=1.,
                 lambdaB2=1.,
                 rank,
                 A.initial=NULL,B.initial=NULL,
                 max.iter=20,conv=10^-2,
                 mode = c("glmnet", "ECOS", "CD")){
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
    #cancor_regul <- cancor(as.matrix(X),as.matrix(Y))
    #A.initial_ridge<-matrix(cancor_regul$xcoef[,1:rank],ncol=rank,nrow=ncol(X))
    #B.initial_ridge<-matrix(cancor_regul$ycoef[,1:rank],ncol=rank,nrow=ncol(Y))
    A.initial_ridge<- matrix(rnorm(rank *ncol(X)),ncol=rank,nrow=ncol(X))
    B.initial_ridge<- matrix(rnorm(rank *ncol(Y)),ncol=rank,nrow=ncol(Y))
    A.initial<-apply(A.initial_ridge,2, NORMALIZATION_UNIT) 
    B.initial<-apply(B.initial_ridge,2,NORMALIZATION_UNIT) 
  } 
  
  # Sequentially extract the r canonical vector pairs
  for (i.r in 1:rank){
    print(paste0("Dimension is " ,i.r))
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
    diff.obj<-conv*10
    
    # FROM i until convergence: canonical vectors 
    while( (it<max.iter) & (diff.obj>conv) ){
      print(it)
      
      # Estimating A conditional on B
      if (is.null(Da) & is.null(lambdaA1) & (p<n)){
        #### don't use regularisation %%%% simple OLS solutioj
        FIT.A <- inv(t(X)%*%X) %*% t(X) %*% Y_data%*%B.STARTING
        
      }else{
        if(is.null(Da)){
          glm.model = glmnet(X_data, Y_data%*%B.STARTING, lambda = lambdaA1, 
                        intercept = FALSE )
          FIT.A <- as.matrix(glm.model$beta)
          #FIT.A <- alternating.regression(Xreg=X_data,Yreg=Y_data%*%B.STARTING,lambdaseq=c(lambdaA1),
          #                                selection.criterion=1, n.cv=1)
        }
        else{
          #### add cross validation later
          X_tilde = X_data %*% pinv(Da)
          glm.model = glmnet(X_tilde, Y_data%*%B.STARTING, 
                        lambda = lambdaA1,
                        intercept = FALSE )
          FIT.A <- pinv(Da) %*% as.matrix(glm.model$beta)
          # FIT.A<-genglm(X_data, Y_data%*%B.STARTING, family="gaussian",
          #               D=Da,lambda1=lambdaA1, lambda2=lambdaA2,
          #               solver="ECOS")
        }
      }
      
      AHAT_FINAL <- FIT.A
      # Estimating B conditional on A
      
      # Estimating A conditional on B
      if (is.null(Db) & is.null(lambdaB1) & (q<n)){
        #### don't use regularisation %%%% simple OLS solutioj
        FIT.B <- inv(t(Y)%*%Y) %*% t(Y) %*% X_data%*%AHAT_FINAL
        
      }else{
        
        if(is.null(Da)){
          glm.model = glmnet(Y_data, X_data%*%AHAT_FINAL, lambda = lambdaB1, 
                        intercept = FALSE )
          FIT.B <- as.matrix(glm.model$beta)
          #FIT.A <- alternating.regression(Xreg=X_data,Yreg=Y_data%*%B.STARTING,lambdaseq=c(lambdaA1),
          #                                selection.criterion=1, n.cv=1)
        }
        else{
          #### add cross validation later
          Y_tilde = Y_data %*% pinv(Db)
          glm.model = glmnet(Y_tilde, X_data%*%AHAT_FINAL, 
                             lambda = lambdaA1,
                             intercept = FALSE )
          FIT.B <- pinv(Db) %*% as.matrix(glm.model$beta)
          # FIT.A<-genglm(X_data, Y_data%*%B.STARTING, family="gaussian",
          #               D=Da,lambda1=lambdaA1, lambda2=lambdaA2,
          #               solver="ECOS")
        }
      }
      
      BHAT_FINAL<-FIT.B
      

      #lambdaB_ALL[it,i.r]<-FIT.B$lambda2
      
      # Check convergence
      obj.new<-mean((X_data %*% matrix(AHAT_FINAL,ncol=1)-Y_data%*%matrix(BHAT_FINAL,ncol=1))^2)
      obj.it[it+1,i.r]<-obj.new
      diff.obj<-abs(obj.new-obj.initial)/obj.initial
      
      # Updated starting values 
      B.STARTING<-BHAT_FINAL
      A.STARTING<-AHAT_FINAL
      obj.initial<-obj.new
      it<-it+1
      if((mean(BHAT_FINAL) == BHAT_FINAL[1]) || (mean(AHAT_FINAL) == AHAT_FINAL[1])){
        it = max.iter + 1
      }
    } # end while-loop
    
    # Number of ITERATIONS
    iterations[1,i.r]<-it
    
    # CANONICAL VARIATES after convergence
    Uhat<-X_data%*%AHAT_FINAL
    Vhat<-Y_data%*%BHAT_FINAL
    
    
    # Express canonical vectors in terms of ORIGINAL DATA MATRICES
    # if (i.r==1){# FIRST DIMENSION
    #   print("First dimension")
    #   # Final estimates of canonical vectors, variates and canonical correlation
    #   ALPHA_ALL[,i.r]<-AHAT_FINAL
    #   BETA_ALL[,i.r]<-BHAT_FINAL
    #   U_ALL[,i.r]<-Uhat
    #   V_ALL[,i.r]<-Vhat
    #   ALPHA_trans[,i.r]<-AHAT_FINAL
    #   BETA_trans[,i.r]<-BHAT_FINAL
    #   U_trans[,i.r]<-Uhat
    #   V_trans[,i.r]<-Vhat
    #   cancors[1,i.r]<-cor(Uhat,Vhat)
    #   cancors_trans[1,i.r]<-cor(Uhat,Vhat)
    #   
    #   # Deflated data matrices
    #   X_data<- round(X_data  - Uhat%*%solve(t(Uhat)%*%Uhat)%*%t(Uhat)%*%X_data,10)
    #   Y_data<- round(Y_data - Vhat%*%solve(t(Vhat)%*%Vhat)%*%t(Vhat)%*%Y_data,10)
    #   
    #   # Sparsity parameters
    #   lambdaA_FINAL<- lambdaA1
    #   lambdaaA2_FINAL<- lambdaA2
    #   lambdaB_FINAL<- lambdaB1
    #   lambdaaB2_FINAL<- lambdaB2
    #   
    #   
    # } else {# HIGHER ORDER DIMENSIONS
    #   
    #   ALPHA_ALL[,i.r]<-ALPHAhat
    #   BETA_ALL[,i.r]<-BETAhat
    #   Uhat<-X%*%ALPHAhat
    #   Vhat<-Y%*%BETAhat
    #   U_ALL[,i.r]<-Uhat
    #   V_ALL[,i.r]<-Vhat
    #   cancors[1,i.r]<-cor(Uhat,Vhat)
    #   # A expressed in terms of original data set X
    #   
    #   # Estimating A conditional on B
    #   if (sum(abs(Uhat)) == 0){
    #     FIT.Aorig <- matrix(0, nrow = p, ncol = 1)
    #   }else{
    #   if (is.null(Da) & is.null(lambdaA1) & p<n){
    #     #### don't use regularisation
    #     FIT.Aorig  <- inv(t(X)%*%X) %*% t(X) %*% Uhat
    #   }else{
    #     if(is.null(Da)){
    #       glm.model = glmnet(X, Uhat, lambda = lambdaA1, 
    #                          intercept = FALSE )
    #       FIT.Aorig <- as.matrix(glm.model$beta)
    #     }
    #     else{
    #       #### add cross validation later
    #       X_tilde = X %*% pinv(Da)
    #       glm.model = glmnet(X_tilde, Uhat, 
    #                          lambda = lambdaA1,
    #                          intercept = FALSE )
    #       FIT.Aorig <- pinv(Da) %*% as.matrix(glm.model$beta)
    #       
    #     }
    #   }
    #   }
      # ALPHAhat<-FIT.Aorig
      # 
      # # B expressed in terms of original data set Y
      # if (sum(abs(Vhat)) == 0){
      #   FIT.Borig <- matrix(0, nrow = q, ncol = 1)
      # }else{
      #   if (is.null(Db) & is.null(lambdaB1) & q<n){
      #     #### don't use regularisation
      #     FIT.Borig  <- inv(t(Y)%*%Y) %*% t(Y) %*% Vhat
      #   }else{
      #     if(is.null(Db)){
      #       glm.model = glmnet(Y, Vhat, lambda = lambdaB1, 
      #                          intercept = FALSE )
      #       FIT.Borig <- as.matrix(glm.model$beta)
      #     }
      #     else{
      #       #### add cross validation later
      #       Y_tilde = Y %*% pinv(Db)
      #       glm.model = glmnet(Y_tilde, Vhat, 
      #                          lambda = lambdaB1,
      #                          intercept = FALSE )
      #       FIT.Borig <- pinv(Db) %*% as.matrix(glm.model$beta)
      #       
      #     }
      #   }
      # }
      # 
      # BETAhat<-FIT.Borig
      # 
      # 
      # # Final estimates of canonical vectors, variates and canonical correlation
      # 
      # ALPHA_trans[,i.r]<-ALPHAhat
      # BETA_trans[,i.r]<-BETAhat
      # Uhat_trans<-X%*%ALPHAhat
      # Vhat_trans<-Y%*%BETAhat
      # U_trans[,i.r]<-Uhat_trans
      # V_trans[,i.r]<-Vhat_trans
      # cancors_trans[1,i.r]<-cor(Uhat_trans,Vhat_trans)
      # 

      
      # Deflated data matrices: regress original data sets on all previously found canonical variates
      #X_data<- X  - U_ALL[,1:i.r]%*%solve(t(U_ALL[,1:i.r])%*%U_ALL[,1:i.r])%*%t(U_ALL[,1:i.r])%*%X
      #Y_data<- Y -  V_ALL[,1:i.r]%*%solve(t(V_ALL[,1:i.r])%*%V_ALL[,1:i.r])%*%t(V_ALL[,1:i.r])%*%Y   
    #}
    #print("Here")
  } # END FOR-LOOP
  ##OUTPUT
  out<-list(ALPHA=ALPHA_ALL,BETA=BETA_ALL,cancors=cancors,U_ALL=U_ALL,V_ALL=V_ALL,
            it=iterations)
}



alternating.gen.regression<-function(Xreg,Yreg,D, 
                                     lambda1,
                                     lambda2){
  ### Function to perform sparse alternating regression
  
  ### INPUT
  #Xreg               : design matrix
  #Yreg               : response
  #lambdaseq          : sequence of sparsity parameters
  #selection.criterion: 1 for BIC, 2 for Cross-validation to select sparsity parameter
  #n.cv               : n.cv-fold cross-validation
  
  ### OUTPUT
  #COEF_FINAL         : estimated coefficients
  #LAMBDA_FINAL       : optimal sparsity parameter
  
  
  ##Standardize
  Xreg_st<-matrix(stdize(Xreg),ncol=ncol(Xreg))
  for (i.variable in 1:ncol(Xreg)){
    #### replace missing variables
    if (is.na(apply(Xreg_st,2,sum)[i.variable])==T) {
      Xreg_st[,i.variable]<-0}    
  }
  
  ##GEN FIT
  GENFIT<-genglm(y=Yreg,x=Xreg_st,D=D, family="gaussian",
                 lambda1=lambda1, lambda2=lambda2, 
                 solver = "ECOS")
  COEFhat<-matrix(GENFIT$beta, nrow=ncol(Xreg_st))
  # if (is.integer(which(GENFIT$df!=0)) && length(which(GENFIT$df!=0)) == 0L) {
  #   # Smaller lambda sequence necessary
  #   GENFIT<-genglm(y=Yreg,x=Xreg_st, D=D, family="gaussian",lambda1=lambda1, lambda2=lambda2, 
  #                    solver = "ECOS")
  #   COEFhat<-matrix(GENFIT$beta[,which(GENFIT$df!=0)[1:length(lambdaseq)]],nrow=ncol(Xreg_st)) # estimated coefficients
  #   LAMBDA<-GENFIT$lambda[which(GENFIT$df!=0)[1:length(lambdaseq)]] # lambda values
  # } else {
  #   COEFhat<-matrix(GENFIT$beta[,which(GENFIT$df!=0)],nrow=ncol(Xreg_st)) # estimated coefficients
  #   LAMBDA<-GENFIT$lambda[which(GENFIT$df!=0)] # lambda values
  # }
  
  ##OUTPUT
  out<-list(COEF_FINAL=COEFhat,LAMBDA1_FINAL=lambda1,  LAMBDA2_FINAL=lambda2)
}

library(parallel)
library(MASS)

cv.genCCA <- function(X, Y, D, rank, 
                      nb.folds, lambdaA1seq,
                      lambdaB1seq, lambdaA2seq, lambdaB2seq){
  
  
  numCores <- detectCores()
  fx <- function(i){
    res = genCCA(X, Y, Da, Db, lambdaA1=lambdaA1seq[i],
                     lambdaB1=lambdaB1seq[i],
                     lambdaA2=lambdaB2seq[i],
                     lambdaB2=lambdaB2seq[i],
                     rank, selection.criterion=1,n.cv=5,
                     A.initial=NULL,B.initial=NULL,max.iter=20,conv=10^-2)
  }
  results <- mclapply(starts, fx, mc.cores = numCores)
  
  # cluster <- makeCluster(numberOfThreads = 3)
  # clusterApply(cluster, 1:10, fun)
  # stopCluster(cluster)
  
}
NORMALIZATION_UNIT<-function(U){
  ### Normalize a vector U to have norm one
  
  length.U<-as.numeric(sqrt(t(U)%*%U))
  if(length.U==0){length.U<-1}
  Uunit<-U/length.U
}


BIC<-function(U,Y.data,X.data){
  ### Calculate value of Bayesian information Criterion
  
  NEGLOGLIK<-sum(diag((1/nrow(Y.data))*(t(Y.data-X.data%*%U)%*%(Y.data-X.data%*%U))))
  BIC<-2*NEGLOGLIK+ length(which(U!=0))*log(nrow(Y.data))
  return(BIC)
}

testsample.correlation<-function(U,Xdata,yscore){
  ### Calculate correlation in test sample
  xscore = Xdata%*%U
  if (all(U==0)) {
    return(0)
  } else {
    return(abs(cor(xscore,yscore)))
  } 
}


############################################################
# CANONICAL RIDGE (used as starting value in SAR algorithm #
############################################################

estim.regul_crossvalidation<-function (X, Y, lambda1grid = NULL, lambda2grid = NULL,n.cv=5){
  # Function to perform Canonical Ridge
  
  ### INPUT
  # X                   : (nxp) data matrix
  # Y                   : (nxq) data matrix
  # lambda1grid         : grid of tuning parameters for lambda1
  # lambda2grid         : grid of tuning parameters for lambda2
  # n.cv                : n.cv-fold cross-validation
  
  ### OUTPUT
  #lambda1.optim        : value of the sparsity parameter lambda1
  #lambda2.optim        : value of the sparsity parameter lambda2
  #cv.optim             : value of cross-validation score
  
  
  if (is.null(lambda1grid)) {
    lambda1grid = matrix(seq(0.001, 1, length = 5),nrow=1)
  }
  if (is.null(lambda2grid)) {
    lambda2grid = matrix(seq(0.001, 1, length = 5),nrow=1)
  }
  
  lambda1.matrix<-matrix(rep(lambda1grid,length(lambda1grid)),ncol=length(lambda2grid),byrow=T)
  lambda2.matrix<-matrix(sort(rep(lambda2grid,length(lambda2grid))),ncol=length(lambda1grid),byrow=T)
  
  
  cvscores<-apply(lambda1grid,2,l1function,Xmatrix=X,Ymatrix=Y,lambda2grid=lambda2grid,n.cv=n.cv) #cv-score
  cv.optim=cvscores[which.max(cvscores)]
  lambda1.optim=lambda1.matrix[which.max(cvscores)]
  lambda2.optim=lambda2.matrix[which.max(cvscores)]
  
  ##OUTPUT
  out=list(lambda1.optim=lambda1.optim,lambda2.optim=lambda1.optim,cv.optim=cv.optim)
}

l1function<-function(U,Xmatrix,Ymatrix,lambda2grid,n.cv){ # AUXILIARY FUNCTION CANONICAL RIDGE
  testcor=apply(lambda2grid, 2, l2function, Xmatrix=Xmatrix,Ymatrix=Ymatrix,lambda1fixed=U,n.cv=n.cv)
  return(testcor)
}

l2function<-function(V,Xmatrix,Ymatrix,lambda1fixed,n.cv){ # AUXILIARY FUNCTION CANONICAL RIDGE
  RCCFIT<-RCC_crossvalidation(X=Xmatrix, Y=Ymatrix, lambda1=lambda1fixed, lambda2=V,n.cv=n.cv) 
  return(RCCFIT$cv)
}
