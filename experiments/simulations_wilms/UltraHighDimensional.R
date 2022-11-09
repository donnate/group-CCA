rm(list=ls())

# LOAD LIBRARIES
source("http://bioconductor.org/biocLite.R") 
biocLite("impute") #After executing this line, indicate if you want to update all/some/none of the packages by typing in a/s/n.
library(PMA)
library(mvtnorm)
library(glmnet)
library(CCA)
library(pls)

# LOAD Functions
source('../Functions/SAR.R')
source('../Functions/Parkhomenko.R')
source('../Functions/Witten_CrossValidation.R')
source('../Functions/Waaijenborg.R')


###########################################################################
# Simulation study Canonical Correlations: Sparse High-dimensional design #
###########################################################################

# Design
Nsim<-200
p<-1000
q<-1000
n<-50
rank<-3

# Covariance matrices
SIGMA_XX<-matrix(0,ncol=p,nrow=p)
SIGMA_XX[1:2,1:2]<-0.009
SIGMA_XX[4:5,4:5]<-0.009
SIGMA_XX[7:8,7:8]<-0.009
diag(SIGMA_XX)<-0.01
SIGMA_YY<-matrix(0,ncol=q,nrow=q)
SIGMA_YY[1:3,1:3]<-0.008
SIGMA_YY[4:6,4:6]<-0.008
SIGMA_YY[7:9,7:9]<-0.008
diag(SIGMA_YY)<-0.01

SIGMA_XY<-matrix(0,nrow=p,ncol=q)
SIGMA_XY[1,1]<-0.002
SIGMA_XY[4,4]<-0.002
SIGMA_XY[7,7]<-0.002
SIGMA_BIG<-rbind(cbind(SIGMA_XX,SIGMA_XY),cbind(t(SIGMA_XY),SIGMA_YY))


# Calculare true canonical vectors
SigmaXX_decomp<-eigen(SIGMA_XX)
SigmaXX_negsqrt<-SigmaXX_decomp$vectors%*%diag(1/sqrt(SigmaXX_decomp$values),length(SigmaXX_decomp$values))%*%t(SigmaXX_decomp$vectors)
SigmaYY_decomp<-eigen(SIGMA_YY)
SigmaYY_negsqrt<-SigmaYY_decomp$vectors%*%diag(1/sqrt(SigmaYY_decomp$values),length(SigmaYY_decomp$values))%*%t(SigmaYY_decomp$vectors)

trueCCA=svd(SigmaXX_negsqrt%*%SIGMA_XY%*%SigmaYY_negsqrt)
true_a<-SigmaXX_negsqrt%*%trueCCA$u[,1:rank]
true_b<-SigmaYY_negsqrt%*%trueCCA$v[,1:rank]


# Store results Estimation Accuracy
MSEa<-matrix(NA,ncol=8,nrow=Nsim)
colnames(MSEa)<-c("SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author")
MSEb<-matrix(NA,ncol=8,nrow=Nsim)
colnames(MSEb)<-c("SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author")

# Store results Sparsity Recognition Performance
TPRa<-matrix(NA,ncol=7,nrow=Nsim)
colnames(TPRa)<-c("SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author")
TPRb<-matrix(NA,ncol=7,nrow=Nsim)
colnames(TPRb)<-c("SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author")
TNRa<-matrix(NA,ncol=7,nrow=Nsim)
colnames(TNRa)<-c("SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author")
TNRb<-matrix(NA,ncol=7,nrow=Nsim)
colnames(TNRb)<-c("SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author")


## START SIMULATION RUNS
for (isim in 1:Nsim){
  cat("Start simulation Run",isim, "\n")
  
  # Seed
  set.seed(isim)
  
  
  # Generate Data
  DATA<-rmvnorm(n,mean=rep(0,p+q),sigma=SIGMA_BIG)
  X<-DATA[,1:p]
  Y<-DATA[,(p+1):(p+q)]
  # Centering data matrices
  X<-(X-matrix(apply(X,2,mean),byrow=TRUE,ncol=ncol(X),nrow=nrow(X)))
  Y<-(Y-matrix(apply(Y,2,mean),byrow=TRUE,ncol=ncol(Y),nrow=nrow(Y)))
  
  # Canonical Ridge
  RCC_cv<-estim.regul_crossvalidation(X,Y,n.cv=5,lambda1grid=matrix(seq(0.1, 0.01, length = 2),nrow=1),lambda2grid=matrix(seq(0.1, 0.01, length = 2),nrow=1))
  RCCA<-rcc(X,Y,RCC_cv$lambda1.optim,RCC_cv$lambda2.optim)

  MSEa[isim,8]<-principal_angles(true_a,RCCA$xcoef[,1:rank])$angles[1,1]
  MSEb[isim,8]<-principal_angles(true_b,RCCA$ycoef[,1:rank])$angles[1,1]

  ### SAR - BIC
  FIT_SAR_BIC<-SparseCCA(X=X,Y=Y,rank=rank,lambdaAseq=seq(from=10^-4,to=10^-6,length=10),lambdaBseq=seq(from=10^-4,to=10^-6,length=10),max.iter=100,conv=10^-2,selection.criterion=1,n.cv=5,A.initial=matrix(RCCA$xcoef[,1:rank],ncol=rank),B.initial=matrix(RCCA$ycoef[,1:rank],ncol=rank))
  
  MSEa[isim,1]<-principal_angles(true_a,FIT_SAR_BIC$ALPHA)$angles[1,1]
  MSEb[isim,1]<-principal_angles(true_b,FIT_SAR_BIC$BETA)$angles[1,1]
  cancors_SAR_BIC<-FIT_SAR_BIC$cancors
  TPRa[isim,1]<-TPR(true_a,FIT_SAR_BIC$ALPHA[,order(cancors_SAR_BIC,decreasing=T)])
  TPRb[isim,1]<-TPR(true_b,FIT_SAR_BIC$BETA[,order(cancors_SAR_BIC,decreasing=T)])
  TNRa[isim,1]<-TNR(true_a,FIT_SAR_BIC$ALPHA[,order(cancors_SAR_BIC,decreasing=T)])
  TNRb[isim,1]<-TNR(true_b,FIT_SAR_BIC$BETA[,order(cancors_SAR_BIC,decreasing=T)])
  
  
  
  ### SAR - CV
  FIT_SAR_CV<-SparseCCA(X=X,Y=Y,rank=rank,lambdaAseq=seq(from=10^-4,to=10^-6,length=10),lambdaBseq=seq(from=10^-4,to=10^-6,length=10),max.iter=100,conv=10^-2,selection.criterion=2,n.cv=5,A.initial=matrix(RCCA$xcoef[,1:rank],ncol=rank),B.initial=matrix(RCCA$ycoef[,1:rank],ncol=rank))
  
  MSEa[isim,2]<-principal_angles(true_a,FIT_SAR_CV$ALPHA)$angles[1,1]
  MSEb[isim,2]<-principal_angles(true_b,FIT_SAR_CV$BETA)$angles[1,1]
  cancors_SAR_CV<-FIT_SAR_CV$cancors
  TPRa[isim,2]<-TPR(true_a,FIT_SAR_CV$ALPHA[,order(cancors_SAR_CV,decreasing=T)])
  TPRb[isim,2]<-TPR(true_b,FIT_SAR_CV$BETA[,order(cancors_SAR_CV,decreasing=T)])
  TNRa[isim,2]<-TNR(true_a,FIT_SAR_CV$ALPHA[,order(cancors_SAR_CV,decreasing=T)])
  TNRb[isim,2]<-TNR(true_b,FIT_SAR_CV$BETA[,order(cancors_SAR_CV,decreasing=T)])
  
  
  
  ### Witten - Permutation approach
  Witten_Perm<-CCA.permute(x=X,z=Y,typex="standard",typez="standard",nperms=50)
  WittenSCCA_Perm<-CCA(x=X,z=Y,typex="standard",typez="standard",K=rank,penaltyx=Witten_Perm$bestpenaltyx,penaltyz=Witten_Perm$bestpenaltyz,trace=F)
  
  MSEa[isim,3]<-principal_angles(true_a,WittenSCCA_Perm$u)$angles[1,1]
  MSEb[isim,3]<-principal_angles(true_b,WittenSCCA_Perm$v)$angles[1,1]
  cancors_Witten_Perm<-WittenSCCA_Perm$cors
  TPRa[isim,3]<-TPR(true_a,WittenSCCA_Perm$u[,order(cancors_Witten_Perm,decreasing=T)])
  TPRb[isim,3]<-TPR(true_b,WittenSCCA_Perm$v[,order(cancors_Witten_Perm,decreasing=T)])
  TNRa[isim,3]<-TNR(true_a,WittenSCCA_Perm$u[,order(cancors_Witten_Perm,decreasing=T)])
  TNRb[isim,3]<-TNR(true_b,WittenSCCA_Perm$v[,order(cancors_Witten_Perm,decreasing=T)])
  
  
  
  ### Witten  - CV
  Witten_CV<-Witten.CV(X=X,Y=Y,n.cv=5,lambdax=matrix(seq(from=0,to=1,length=20),nrow=1),lambday=matrix(seq(from=0,to=1,length=20),nrow=1))
  WittenSCCA_CV<-CCA(x=X,z=Y,typex="standard",typez="standard",K=rank,penaltyx=Witten_CV$lambdax.opt,penaltyz=Witten_CV$lambday.opt,trace=F)
  
  MSEa[isim,4]<-principal_angles(true_a,WittenSCCA_CV$u)$angles[1,1]
  MSEb[isim,4]<-principal_angles(true_b,WittenSCCA_CV$v)$angles[1,1]
  cancors_Witten_CV<-WittenSCCA_CV$cors
  TPRa[isim,4]<-TPR(true_a,WittenSCCA_CV$u[,order(cancors_Witten_CV,decreasing=T)])
  TPRb[isim,4]<-TPR(true_b,WittenSCCA_CV$v[,order(cancors_Witten_CV,decreasing=T)])
  TNRa[isim,4]<-TNR(true_a,WittenSCCA_CV$u[,order(cancors_Witten_CV,decreasing=T)])
  TNRb[isim,4]<-TNR(true_b,WittenSCCA_CV$v[,order(cancors_Witten_CV,decreasing=T)])
  
  
  
  ### Waaijenborg  - Minimize difference between test and training sample correlation
  Waaijenborg_Delta<-Waaijenborg(X=X,Y=Y,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),rank=rank,selection.criterion=1)
  
  MSEa[isim,5]<-principal_angles(true_a,Waaijenborg_Delta$vhat)$angles[1,1]
  MSEb[isim,5]<-principal_angles(true_b,Waaijenborg_Delta$uhat)$angles[1,1]
  cancors_Waaijenborg_Delta<-Waaijenborg_Delta$cancor
  TPRa[isim,5]<-TPR(true_a,Waaijenborg_Delta$vhat[,order(cancors_Waaijenborg_Delta,decreasing=T)])
  TPRb[isim,5]<-TPR(true_b,Waaijenborg_Delta$uhat[,order(cancors_Waaijenborg_Delta,decreasing=T)])
  TNRa[isim,5]<-TNR(true_a,Waaijenborg_Delta$vhat[,order(cancors_Waaijenborg_Delta,decreasing=T)])
  TNRb[isim,5]<-TNR(true_b,Waaijenborg_Delta$uhat[,order(cancors_Waaijenborg_Delta,decreasing=T)])
  
  
  
  ### Waaijenborg  - Maximize  test sample correlation
  Waaijenborg_Test<-Waaijenborg(X=X,Y=Y,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),rank=rank,selection.criterion=2)
  
  MSEa[isim,6]<-principal_angles(true_a,Waaijenborg_Test$vhat)$angles[1,1]
  MSEb[isim,6]<-principal_angles(true_b,Waaijenborg_Test$uhat)$angles[1,1]
  cancors_Waaijenborg_Test<-Waaijenborg_Test$cancor
  TPRa[isim,6]<-TPR(true_a,Waaijenborg_Test$vhat[,order(cancors_Waaijenborg_Test,decreasing=T)])
  TPRb[isim,6]<-TPR(true_b,Waaijenborg_Test$uhat[,order(cancors_Waaijenborg_Test,decreasing=T)])
  TNRa[isim,6]<-TNR(true_a,Waaijenborg_Test$vhat[,order(cancors_Waaijenborg_Test,decreasing=T)])
  TNRb[isim,6]<-TNR(true_b,Waaijenborg_Test$uhat[,order(cancors_Waaijenborg_Test,decreasing=T)])
  
  
  
  #### Parkhomenko    
  Parkhomenko_SCCA<-SCCA_Parkhomenko(x.data=X,y.data=Y,Krank=rank,lambda.v.seq=seq(0, 0.15, length=20), lambda.u.seq=seq(0, 0.1, length=20))

  
  MSEa[isim,7]<-principal_angles(true_a,Parkhomenko_SCCA$a)$angles[1,1]
  MSEb[isim,7]<-principal_angles(true_b,Parkhomenko_SCCA$b)$angles[1,1]
  cancors_Parkhomenko<-Parkhomenko_SCCA$cancor
  TPRa[isim,7]<-TPR(true_a,-Parkhomenko_SCCA$a[,order(cancors_Parkhomenko,decreasing=T)])
  TPRb[isim,7]<-TPR(true_b,Parkhomenko_SCCA$b[,order(cancors_Parkhomenko,decreasing=T)])
  TNRa[isim,7]<-TNR(true_a,Parkhomenko_SCCA$a[,order(cancors_Parkhomenko,decreasing=T)])
  TNRb[isim,7]<-TNR(true_b,Parkhomenko_SCCA$b[,order(cancors_Parkhomenko,decreasing=T)])

}

## Results
# Results Estimation accucary Estimator A
round(apply(MSEa,2,mean),3)
# Results Estimation accucary Estimator B
round(apply(MSEb,2,mean),3)
# Results Sparsity Recognition Performance Estimator A
round(apply(TPRa,2,mean),2)
round(apply(TNRa,2,mean),2)
# Results Sparsity Recognition Performance Estimator B
round(apply(TPRb,2,mean),2)
round(apply(TNRb,2,mean),2)
