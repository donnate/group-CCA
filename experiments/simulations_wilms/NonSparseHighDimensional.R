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
Nsim<-1000
p<-25
q<-40
n<-50
rank<-1

# Covariance matrices
SIGMA_XX<-matrix(0.1,ncol=p,nrow=p)
diag(SIGMA_XX)<-1
SIGMA_YY<-matrix(0.1,ncol=q,nrow=q)
diag(SIGMA_YY)<-1
SIGMA_XY<-matrix(0,nrow=p,ncol=q)
SIGMA_XY[1,1]<-0.9
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
MSEa<-matrix(NA,ncol=9,nrow=Nsim)
colnames(MSEa)<-c("SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author","CCA")
MSEb<-matrix(NA,ncol=9,nrow=Nsim)
colnames(MSEb)<-c("SAR-Author","SAR-CV","Witten-Author","Witten-CV","Waaijenborg-Author","Waaijenborg-CV","Parkhomenko-Author","Canonical Ridge-Author","CCA")


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
  
  ### SAR - BIC
  FIT_SAR_BIC<-SparseCCA(X=X,Y=Y,rank=rank,lambdaAseq=seq(from=0.2,to=0.02,length=10),lambdaBseq=seq(from=0.2,to=0.02,length=10),max.iter=100,conv=10^-2,selection.criterion=1,n.cv=5)
  
  MSEa[isim,1]<-principal_angles(true_a,FIT_SAR_BIC$ALPHA)$angles[1,1]
  MSEb[isim,1]<-principal_angles(true_b,FIT_SAR_BIC$BETA)$angles[1,1]
  cancors_SAR_BIC<-FIT_SAR_BIC$cancors

  
  
  
  ### SAR - CV
  FIT_SAR_CV<-SparseCCA(X=X,Y=Y,rank=rank,lambdaAseq=seq(from=0.2,to=0.02,length=10),lambdaBseq=seq(from=0.2,to=0.02,length=10),max.iter=100,conv=10^-2,selection.criterion=2,n.cv=5)
  
  MSEa[isim,2]<-principal_angles(true_a,FIT_SAR_CV$ALPHA)$angles[1,1]
  MSEb[isim,2]<-principal_angles(true_b,FIT_SAR_CV$BETA)$angles[1,1]
  cancors_SAR_CV<-FIT_SAR_CV$cancors

  
  
  ### Witten - Permutation approach
  Witten_Perm<-CCA.permute(x=X,z=Y,typex="standard",typez="standard",nperms=50)
  WittenSCCA_Perm<-CCA(x=X,z=Y,typex="standard",typez="standard",K=rank,penaltyx=Witten_Perm$bestpenaltyx,penaltyz=Witten_Perm$bestpenaltyz,trace=F)
  
  MSEa[isim,3]<-principal_angles(true_a,WittenSCCA_Perm$u)$angles[1,1]
  MSEb[isim,3]<-principal_angles(true_b,WittenSCCA_Perm$v)$angles[1,1]
  cancors_Witten_Perm<-WittenSCCA_Perm$cors

  
  
  
  ### Witten  - CV
  Witten_CV<-Witten.CV(X=X,Y=Y,n.cv=5,lambdax=matrix(seq(from=0,to=1,length=20),nrow=1),lambday=matrix(seq(from=0,to=1,length=20),nrow=1))
  WittenSCCA_CV<-CCA(x=X,z=Y,typex="standard",typez="standard",K=rank,penaltyx=Witten_CV$lambdax.opt,penaltyz=Witten_CV$lambday.opt,trace=F)
  
  MSEa[isim,4]<-principal_angles(true_a,WittenSCCA_CV$u)$angles[1,1]
  MSEb[isim,4]<-principal_angles(true_b,WittenSCCA_CV$v)$angles[1,1]
  cancors_Witten_CV<-WittenSCCA_CV$cors

  
  
  
  ### Waaijenborg  - Minimize difference between test and training sample correlation
  Waaijenborg_Delta<-Waaijenborg(X=X,Y=Y,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),rank=rank,selection.criterion=1)
  
  MSEa[isim,5]<-principal_angles(true_a,Waaijenborg_Delta$vhat)$angles[1,1]
  MSEb[isim,5]<-principal_angles(true_b,Waaijenborg_Delta$uhat)$angles[1,1]
  cancors_Waaijenborg_Delta<-Waaijenborg_Delta$cancor

  
  
  ### Waaijenborg  - Maximize  test sample correlation
  Waaijenborg_Test<-Waaijenborg(X=X,Y=Y,lambdaxseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),lambdayseq=matrix(seq(from=0.1,to=5,length=50),nrow=1),rank=rank,selection.criterion=2)
  
  MSEa[isim,6]<-principal_angles(true_a,Waaijenborg_Test$vhat)$angles[1,1]
  MSEb[isim,6]<-principal_angles(true_b,Waaijenborg_Test$uhat)$angles[1,1]
  cancors_Waaijenborg_Test<-Waaijenborg_Test$cancor

  
  
  #### Parkhomenko    
  Parkhomenko_SCCA<-SCCA_Parkhomenko(x.data=X,y.data=Y,Krank=rank)
  
  MSEa[isim,7]<-principal_angles(true_a,Parkhomenko_SCCA$a)$angles[1,1]
  MSEb[isim,7]<-principal_angles(true_b,Parkhomenko_SCCA$b)$angles[1,1]
  cancors_Parkhomenko<-Parkhomenko_SCCA$cancor

  
  
  # Canonical Ridge
  RCC_cv<-estim.regul_crossvalidation(X,Y,n.cv=5)
  RCCA<-rcc(X,Y,RCC_cv$lambda1.optim,RCC_cv$lambda2.optim)
  
  MSEa[isim,8]<-principal_angles(true_a,RCCA$xcoef[,1:rank])$angles[1,1]
  MSEb[isim,8]<-principal_angles(true_b,RCCA$ycoef[,1:rank])$angles[1,1]
  
  
  
  # Traditional CCA
  trad_cancor<-cancor(X,Y,xcenter=FALSE,ycenter=FALSE)
  
  MSEa[isim,9]<-principal_angles(true_a,trad_cancor$xcoef[,1:rank])$angles[1,1]
  MSEb[isim,9]<-principal_angles(true_b,trad_cancor$ycoef[,1:rank])$angles[1,1]
}

## Results
# Results Estimation accucary Estimator A
round(apply(MSEa,2,mean),3)
# Results Estimation accucary Estimator B
round(apply(MSEb,2,mean),3)

