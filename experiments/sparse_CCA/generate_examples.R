library(Matrix)
library(MASS)
library(geigen)
library(pracma)
library(VGAM)
library(mvtnorm)


generate_simple_example <- function(n, p1, p2, 
                                    nnzeros = 5,
                                    r = 2,
                                    theta = diag( c(0.9,  0.8))){
  ###
  # This generates a dataset (X,Y) from a multivariate normal where Sigma_XX  and 
  # Sigma_YY=I, and Sigma_XY = U Lambda V^T is rank r on a set of nnzeros rows, 0 elsewhere.
  # n: number of samples
  # p1: nb of features for X 
  # p2: nb of features for Y
  # nnzeros: number of non zero rows of U and V
  # r: rank of Sigma_XY
  # theta : canonical correlation (must be a list of length r)
  # Returns:
  # S : empirical covariance matrix X^TY
  # Sigma: underlying (population) covariance Sigma_{XY}
  # u: ground truth for u 
  # v: ground truth for v
  ###
  p_tot <- p1 + p2;
  pp <- c(p1,p2);
  print("We have:")
  print(c(p_tot, nnzeros<=min(p1,p2), nnzeros))
  print('--------------------------------------');
  print('Generating data ...');
  #a <- 0.3;
  Sigma <- diag(p1+p2)
  s = (1):(nnzeros)
  
  # generate covariance matrix for X and Y
  u = matrix(0, pp[1], r)
  u[s,1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  u <- u %*%(sqrtm(t(u[s,1:r])%*% u[s,1:r])$Binv)
  
  v = matrix(0, pp[2], r)
  v[s,1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  v <- v %*%(sqrtm(t(v[s,1:r]) %*% v[s,1:r])$Binv)
  
  ###
  Sigma[(p1+1):(p1+p2), 1:p1] =  v  %*% theta %*% t(u);
  Sigma[1:p1, (p1+1):(p1+p2)] = t(Sigma[(p1+1):(p1+p2), 1:p1])
  
  Sigmax = Sigma[1:p1,1:p1];
  Sigmay = Sigma[(p1+1):p_tot,(p1+1):p_tot];
  
  #Generate Multivariate Normal Data According to Sigma
  sqrt_Sigma = sqrtm(Sigma)$B
  Data = matrix(rnorm(n= n * p_tot), ncol=p_tot) %*% sqrt_Sigma;
  
  X = Data[,1:p1];
  Y = Data[,(p1+1):(p1+p2)];
  
  print('Data generated.');
  print('--------------------------------------');
  
  Mask = matrix(0, p_tot, p_tot);
  idx1 = 1:pp[1];
  idx2 = (pp[1]+1):(pp[1]+pp[2]);
  Mask[idx1,idx1] <- matrix(1,pp[1],pp[1]);
  Mask[idx2,idx2] <- matrix(1,pp[2],pp[2]);
  Sigma0 = Sigma * Mask;
  
  
  S <- cov(Data)
  sigma0hat <- S * Mask;
  
  ### Just to make sure we have the GT
  GT = svd(sqrtm(Sigma[1:p1, 1:p1])$Binv %*% Sigma[1:p1, (p1+1):p_tot]%*% sqrtm(Sigma[(p1+1):(p_tot), (p1+1):(p_tot)])$Binv, 
           nu=r, nv=r)
  
  
  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data, u=GT$u, v=GT$v, 
              Sigmax=Sigmax, Sigmay=Sigmay
  ))
}



generate_toeplitz_example <- function(n, p1, p2,  a = 0.4,
                                      nnzeros = 5){
  ###
  # This generates a dataset (X,Y) from a multivariate normal where Sigma_XX = I
  # Sigma_YY = I, and Sigma_XY = toeplitz(a) on a set of nnzeros rows, 0 elsewhere
  # n: number of samples
  # p1: nb of features for X 
  # p2: nb of features for Y
  # nnzeros: number of non zero rows of U and V
  # a : strength of the Toeplitz parameter
  #
  # Returns:
  # S : empirical covariance matrix X^TY
  # Sigma: underlying (population) covariance Sigma_{XY}
  # u: ground truth for u 
  # v: ground truth for v
  ###
  p_tot <- p1 + p2;
  pp <- c(p1,p2);
  print("We have:")
  print(c(p, nnzeros<=min(p1,p2), nnzeros))
  print('--------------------------------------');
  print('Generating data ...');
  Sigma <- diag(p1+p2)
  s = (1):(nnzeros)
  
  # generate covariance matrix for X and Y
  Tss = toeplitz(a^(s))
  Tss[which(Tss<1e-4)]=0
  Sigma[s, p1+s] <- Tss
  Sigma[p1+s, s] <- t(Tss)
  
  Sigmax = Sigma[1:p1,1:p1];  ### No correlations in Sigma_XX and Sigma_YY
  Sigmay = Sigma[(p1+1):p,(p1+1):p_tot];
  
  #Generate Multivariate Normal Data According to Sigma
  Data = mvrnorm(n, rep(0, p_tot), Sigma);
  
  X = Data[,1:p1];
  Y = Data[,(p1+1):(p1+p2)];
  
  print('Data generated.');
  print('--------------------------------------');
  
  Mask = matrix(0, p_tot, p_tot);
  idx1 = 1:pp[1];
  idx2 = (pp[1]+1):(pp[1]+pp[2]);
  Mask[idx1,idx1] <- matrix(1,pp[1],pp[1]);
  Mask[idx2,idx2] <- matrix(1,pp[2],pp[2]);
  Sigma0 = Sigma * Mask;
  
  
  S <- cov(Data)
  sigma0hat <- S * Mask
  
  # Generate ground truth canonical vectors
  GT = svd(sqrtm(Sigma[1:p1, 1:p1])$Binv %*% Sigma[1:p1, (p1+1):p_tot]%*% sqrtm(Sigma[(p1+1):(p_tot), (p1+1):(p_tot)])$Binv)
  
  #scale <- a %*% sqrtm(diag(r)+t(a) %*% Sigma %*% a/lambda)$B;
  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data, u=GT$u, v=GT$v, 
              Sigmax=Sigmax, Sigmay=Sigmay
  ))
}



generate_example_none_trivial_pca <- function(n, p1, p2, 
                                              r_pca=3, 
                                              nnzeros = 5,
                                              theta = diag( c(0.9,  0.8)),
                                              lambda_pca = 1,
                                              r=2, overlapping_amount = 0,
                                              normalize_diagonal = TRUE){
  ###
  # This generates a dataset (X,Y) from a multivariate normal where Sigma_XX  and 
  # Sigma_YY have a correlation structure, and Sigma_XY = U Lambda V^T is rank r on a set of nnzeros rows, 0 elsewhere.
  # The number of rows (resp. Columns) on which Sigma_XX (resp. Sigma_YY)  and Sigma_XY
  # overlap is controlled by overlapping_amount
  # n: number of samples
  # p1: nb of features for X 
  # p2: nb of features for Y
  # nnzeros: number of non zero rows of U and V
  # r: rank of Sigma_XY
  # theta : canonical correlation (must be a list of length r)
  # r_pca: rank of the PCA (used in the creation of Sigma_XX and SigmaYY)
  # lambda_pca: also used to create Sigma_XX as Sigma_{XX} = U_X Lambda_pca U_X^T, and setting diag(Sigma_XX) to 1
  # Returns:
  # S : empirical covariance matrix X^TY
  # Sigma: underlying (population) covariance Sigma_{XY}
  # u: ground truth for u 
  # v: ground truth for v
  ###
  p_tot <- p1 + p2;
  pp <- c(p1,p2);
  print("We have:")
  print(c(p_tot, nnzeros<=min(p1,p2), nnzeros))
  print('--------------------------------------');
  print('Generating data ...');
  #a <- 0.3;
  Sigma <- diag(p1+p2)
  s = (1):(nnzeros)
  print('Number of non zeros is:');
  print(nnzeros)
  start = ceiling((1- overlapping_amount) * nnzeros)
  s_pca = (start+1):(start + nnzeros)
  s_pca2 = (start + nnzeros +1):(start + 2 * nnzeros)
  
  Lambda_pca = rep(lambda_pca, r_pca)
  
  # generate vcovariance matrix for X and Y
  u1 = matrix(0, p1, r_pca)
  u1[s_pca,] <- matrix(runif(n=nnzeros * r_pca,max = 3, min=1), nrow=nnzeros, ncol=r_pca) # * matrix(sample(c(-1,1), nnzeros*r_pca, replace=TRUE), nrow=nnzeros, ncol=r_pca)
  u1[s_pca,1:r_pca] <- u1[s_pca,1:r_pca] %*%(sqrtm(t(u1[s_pca,1:r_pca]) %*% u1[s_pca,1:r_pca])$Binv)
  T1 = u1 %*% diag(Lambda_pca) %*% t(u1)#as.matrix(toeplitz(a^(0:(pp[1]-1))));
  #T1[which(T1<1e-6)] = 0
  if (normalize_diagonal){
    diag(T1) <- 1
    Sigma[1:p1, 1:p1] = T1;
  }else{
    Sigma[1:p1, 1:p1] =Sigma[1:p1, 1:p1] + T1;
  }
  
  
  #T1 = Sigma[1:p1, 1:p1]
  Tss = Sigma[s,s];
  TsMs = T1[s,-s];
  TMsMs =  T1[-s,-s];
  u = matrix(0, pp[1], r)
  u[s,1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  print("Tss is")
  print(dim(Tss))
  print(Tss)
   print("U is")
  print(dim(u))
  u <- u %*%(sqrtm(t(u[s,1:r]) %*% Tss%*% u[s,1:r])$Binv)
  
  
  
  u2 = matrix(0, pp[2], r_pca)
  u2[s_pca2,1:r_pca] <- as.matrix(runif( nnzeros * r_pca,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r_pca, replace=TRUE), nrow=nnzeros)
  u2[s_pca2,1:r_pca] <- u2[s_pca2,1:r_pca] %*%(sqrtm(t(u2[s_pca2,1:r_pca]) %*% u2[s_pca2,1:r_pca])$Binv)
  T2 = u2 %*% diag(Lambda_pca) %*% t(u2)#as.matrix(toeplitz(a^(0:(pp[1]-1))));
  if (normalize_diagonal){
    diag(T2) <- 1
    Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)] = T2;
  }else{
    Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)] = T2 + Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)];
  }
  
  Tss = Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)][s, s];
  v = matrix(0, pp[2], r)
  v[(max(s)+1):(max(s) + length(s)),1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  v <- v %*%(sqrtm(t(v[s,1:r]) %*% Tss %*% v[s,1:r])$Binv)
  
  ###
  T2 = Sigma[(p1+1):p_tot,(p1+1):p_tot]
  T1 = Sigma[1:p1, 1:p1]
  Sigma[(p1+1):(p1+p2), 1:p1] <- T2 %*%  v  %*% theta %*% t(u) %*% T1;
  Sigma[1:p1, (p1+1):(p1+p2)] <- t(Sigma[(p1+1):(p1+p2), 1:p1])
  
  Sigmax = Sigma[1:p1,1:p1];
  Sigmay = Sigma[(p1+1):p_tot,(p1+1):p_tot];
  
  #Generate Multivariate Normal Data According to Sigma
  sqrt_Sigma = sqrtm(Sigma)$B
  Data = matrix(rnorm(n= n * p_tot), ncol=p_tot) %*% sqrt_Sigma;
  
  X = Data[,1:p1];
  Y = Data[,(p1+1):(p1+p2)];
  
  print('Data generated.');
  print('--------------------------------------');
  
  Mask = matrix(0, p_tot, p_tot);
  idx1 = 1:pp[1];
  idx2 = (pp[1]+1):(pp[1]+pp[2]);
  Mask[idx1,idx1] <- matrix(1,pp[1],pp[1]);
  Mask[idx2,idx2] <- matrix(1,pp[2],pp[2]);
  Sigma0 = Sigma * Mask;
  
  
  S <- cov(Data)
  sigma0hat <- S * Mask;

  
  
  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data, u=u, v=v, 
              Sigmax=Sigmax, Sigmay=Sigmay
  ))
}
