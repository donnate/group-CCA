library(fda)
library(combinat)
library(clue)
library(Matrix)
library(MASS)
library(geigen)
library(pracma)
library(VGAM)
library(mvtnorm)

############## simulation #################

generate =  function(n, p, q, s, prop){
  Sigma = diag(p + q)
  X = matrix(rnorm(n*p), n, p)
  R = matrix(rnorm(p*q), p, q)
  Y = X %*% R
  signal = mean(apply(Y, 2, var))
  Y = Y/sqrt(signal)
  E = matrix(rnorm(n*q, 0, s), n, q)
  Y = Y + E
  X = scale(X, center = TRUE, scale = FALSE)
  Y = scale(Y, center = TRUE, scale = FALSE)
  maskX = sample(1:(n*p), n*p*prop)
  Xna = X
  Xna[maskX] = NA
  maskY = sample(1:(n*q), n*q*prop)
  Yna = Y
  Yna[maskY] = NA
  Sigma_hat_sqrt_inv = sqrtm(t(rbind(X, Y)) %*% rbind(X, Y)/n)$Binv
  Sigma0_sqrt_inv = Sigma_hat_sqrt_inv
  return(list(X = X, Y = Y, Xna = Xna, Yna = Yna,
              Sigma_hat_sqrt_inv =Sigma_hat_sqrt_inv,
              Sigma0_sqrt_inv = Sigma0_sqrt_inv
              ))
}


generate_example_non_trivial_pca <- function(n, p1, p2,
                                             r_pca = 3,
                                             nnzeros = 5,
                                             theta = diag(c(0.9,  0.8)),
                                             lambda_pca = 1,
                                             r = 2, overlapping_amount = 0,
                                             normalize_diagonal = TRUE,
                                             prop_missing = 0) {
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
  p_tot <- p1 + p2
  pp <- c(p1, p2)
  print("We have:")
  print(c(p_tot, nnzeros <= min(p1, p2), nnzeros))
  print('--------------------------------------')
  print('Generating data ...')
  Sigma <- diag(p1 + p2)
  s <- (1):(nnzeros)
  print(paste0('Number of non zeros is: ', nnzeros))
  start <- ceiling((1 - overlapping_amount) * nnzeros)
  
  if (r_pca >0){
    s_pca  <- (start + 1) : (start + nnzeros)
    s_pca2  <- s_pca
    Lambda_pca <- rep(lambda_pca, r_pca)
    # generate vcovariance matrix for X and Y
    u1 = matrix(0, p1, r_pca)
    u1[s_pca, ] <- matrix(runif(n = nnzeros * r_pca, max = 3, min=1), nrow = nnzeros, ncol = r_pca) * matrix(sample(c(-1,1), nnzeros * r_pca, replace=TRUE), nrow=nnzeros, ncol=r_pca)
    # Normalize u1
    u1[s_pca, 1:r_pca] <- u1[s_pca,1:r_pca] %*% (sqrtm(t(u1[s_pca,1:r_pca]) %*% u1[s_pca, 1:r_pca])$Binv)
    T1 <- u1 %*% diag(Lambda_pca) %*% t(u1)
    if (normalize_diagonal){
      diag(T1) <- 1
      Sigma[1:p1, 1:p1] <- T1
    }else{
      Sigma[1:p1, 1:p1] <- Sigma[1:p1, 1:p1] + T1 
    }
    
    ### Same for Sigma_y
    u2 <- matrix(0, pp[2], r_pca)
    u2[s_pca2, 1:r_pca] <- matrix(runif( nnzeros * r_pca,max = 3, min=1), nrow=nnzeros)  * matrix(sample(c(-1,1), nnzeros*r_pca, replace=TRUE), nrow=nnzeros, ncol=r_pca)
    u2[s_pca2, 1:r_pca] <- u2[s_pca2,1:r_pca] %*% (sqrtm(t(u2[s_pca2,1:r_pca]) %*% u2[s_pca2,1:r_pca])$Binv)
    T2 <- u2 %*% diag(Lambda_pca) %*% t(u2)
    if (normalize_diagonal){
      diag(T2) <- 1
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2
    }else{
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2 + Sigma[(p1 + 1):(p1 + p2), (p1 + 1) : (p1 + p2)]
    }
    
  }
  Sigmax = Sigma[1:p1,1:p1];
  Sigmay = Sigma[(p1+1):p_tot,(p1+1):p_tot];
  
  
  ### Generate cross covariance
  Tss <- Sigma[s,s]
  u <- matrix(0, pp[1], r)
  u[s, 1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  u <- u %*% (sqrtm(t(u[s, 1:r]) %*% Tss %*% u[s, 1:r])$Binv)
  Tss_v <- Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)][s, s]
  v <- matrix(0, pp[2], r)
  v[s, 1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  v <- v %*% (sqrtm(t(v[s, 1:r]) %*% Tss_v %*% v[s, 1:r])$Binv)
  Sigma[(p1 + 1) :(p1 + p2), 1:p1] <- Sigmay %*%  v  %*% theta %*% t(u) %*% Sigmax
  Sigma[1:p1, (p1 + 1):(p1 + p2)] <- Sigmax %*%  u  %*% theta %*% t(v) %*% Sigmay
  
  
  
  #Generate Multivariate Normal Data According to Sigma
  sqrt_Sigma <- sqrtm(Sigma)$B
  Data <- matrix(rnorm(n = n * p_tot), ncol=p_tot) %*% sqrt_Sigma
  X <- Data[, 1:p1]
  Y <- Data[, (p1 + 1):(p1 + p2)]
  print('Data generated.')
  print('--------------------------------------')
  Mask <- matrix(0, p_tot, p_tot)
  idx1 <- 1:pp[1]
  idx2 <- (pp[1] + 1):(pp[1] + pp[2])
  Mask[idx1, idx1] <- matrix(1, pp[1], pp[1])
  Mask[idx2, idx2] <- matrix(1, pp[2], pp[2])
  Sigma0 <- Sigma * Mask
  S <- cov(Data)
  sigma0hat <- S * Mask
  # Generate ground truth canonical vectors
  Sigma_X_inv <- solve(Sigma[1:p1, 1:p1])
  Sigma_Y_inv <-  solve(Sigma[(p1 + 1):(p_tot), (p1 + 1):(p_tot)])
  GT = svd(Sigma_X_inv %*% Sigma[1:p1, (p1 + 1):p_tot] %*% Sigma_Y_inv, nu = r, nv = r)
  
  if (prop_missing >0){
    maskX = sample(1:(n*p1), n * p1 * prop_missing)
    Xna = X
    Xna[maskX] = NA
    maskY = sample(1:(n*p2), n * p2 * prop_missing)
    Yna = Y
    Yna[maskY] = NA
  }
  
  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data, u=GT$u, v=GT$v, 
              Xna = Xna, Yna = Yna
  ))
}


generate_example_non_trivial_pca <- function(n, p1, p2,
                                              r_pca = 3,
                                              nnzeros = 5,
                                              theta = diag(c(0.9,  0.8)),
                                              lambda_pca = 1,
                                              r = 2, overlapping_amount = 0,
                                              normalize_diagonal = TRUE,
                                              prop_missing = 0,
                                              noise = 1) {
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
  p_tot <- p1 + p2
  pp <- c(p1, p2)
  print("We have:")
  print(c(p_tot, nnzeros <= min(p1, p2), nnzeros))
  print('--------------------------------------')
  print('Generating data ...')
  Sigma <- diag(x = noise, nrow = p1 + p2)
  s <- (1):(nnzeros)
  print(paste0('Number of non zeros is: ', nnzeros))
  start <- ceiling((1 - overlapping_amount) * nnzeros)
  
  if (r_pca >0){
    s_pca  <- (start + 1) : (start + nnzeros)
    s_pca2  <- s_pca
    Lambda_pca <- rep(lambda_pca, r_pca)
    # generate vcovariance matrix for X and Y
    u1 = matrix(0, p1, r_pca)
    u1[s_pca, ] <- matrix(runif(n = nnzeros * r_pca, max = 3, min=1), nrow = nnzeros, ncol = r_pca) * matrix(sample(c(-1,1), nnzeros * r_pca, replace=TRUE), nrow=nnzeros, ncol=r_pca)
    # Normalize u1
    u1[s_pca, 1:r_pca] <- u1[s_pca,1:r_pca] %*% (sqrtm(t(u1[s_pca,1:r_pca]) %*% u1[s_pca, 1:r_pca])$Binv)
    T1 <- u1 %*% diag(Lambda_pca) %*% t(u1)
    if (normalize_diagonal){
      diag(T1) <- 1
      Sigma[1:p1, 1:p1] <- T1
    }else{
      Sigma[1:p1, 1:p1] <- Sigma[1:p1, 1:p1] + T1 
    }
    
    ### Same for Sigma_y
    u2 <- matrix(0, pp[2], r_pca)
    u2[s_pca2, 1:r_pca] <- matrix(runif( nnzeros * r_pca,max = 3, min=1), nrow=nnzeros)  * matrix(sample(c(-1,1), nnzeros*r_pca, replace=TRUE), nrow=nnzeros, ncol=r_pca)
    u2[s_pca2, 1:r_pca] <- u2[s_pca2,1:r_pca] %*% (sqrtm(t(u2[s_pca2,1:r_pca]) %*% u2[s_pca2,1:r_pca])$Binv)
    T2 <- u2 %*% diag(Lambda_pca) %*% t(u2)
    if (normalize_diagonal){
      diag(T2) <- 1
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2
    }else{
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2 + Sigma[(p1 + 1):(p1 + p2), (p1 + 1) : (p1 + p2)]
    }
  
  }
  Sigmax = Sigma[1:p1,1:p1];
  Sigmay = Sigma[(p1+1):p_tot,(p1+1):p_tot];
  
  
  ### Generate cross covariance
  Tss <- Sigma[s,s]
  u <- matrix(0, pp[1], r)
  u[s, 1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  u <- u %*% (sqrtm(t(u[s, 1:r]) %*% Tss %*% u[s, 1:r])$Binv)
  Tss_v <- Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)][s, s]
  v <- matrix(0, pp[2], r)
  v[s, 1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  v <- v %*% (sqrtm(t(v[s, 1:r]) %*% Tss_v %*% v[s, 1:r])$Binv)
  Sigma[(p1 + 1) :(p1 + p2), 1:p1] <- Sigmay %*%  v  %*% theta %*% t(u) %*% Sigmax
  Sigma[1:p1, (p1 + 1):(p1 + p2)] <- Sigmax %*%  u  %*% theta %*% t(v) %*% Sigmay
  
  
  
  #Generate Multivariate Normal Data According to Sigma
  sqrt_Sigma <- sqrtm(Sigma)$B
  Data <- matrix(rnorm(n = (2 *n) * p_tot), ncol=p_tot) %*% sqrt_Sigma
  X <- Data[1:n, 1:p1]
  Y <- Data[1:n, (p1 + 1):(p1 + p2)]
  newX <- Data[(n+1):(2*n), 1:p1]
  newY <- Data[(n+1):(2*n), (p1 + 1):(p1 + p2)]
  print('Data generated.')
  print('--------------------------------------')
  Mask <- matrix(0, p_tot, p_tot)
  idx1 <- 1:pp[1]
  idx2 <- (pp[1] + 1):(pp[1] + pp[2])
  Mask[idx1, idx1] <- matrix(1, pp[1], pp[1])
  Mask[idx2, idx2] <- matrix(1, pp[2], pp[2])
  Sigma0 <- Sigma * Mask
  S <- cov(Data)
  sigma0hat <- S * Mask
  # Generate ground truth canonical vectors
  Sigma_X_inv <- sqrtm(Sigma[1:p1, 1:p1])$Binv
  Sigma_Y_inv <-  sqrtm(Sigma[(p1 + 1):(p_tot), (p1 + 1):(p_tot)])$Binv
  GT = svd(Sigma_X_inv %*% Sigma[1:p1, (p1 + 1):p_tot] %*% Sigma_Y_inv, nu = r, nv = r)
  
  if (prop_missing >0){
    maskX = sample(1:(n*p1), n * p1 * prop_missing)
    Xna = X
    Xna[maskX] = NA
    maskY = sample(1:(n*p2), n * p2 * prop_missing)
    Yna = Y
    Yna[maskY] = NA
  }else{
    Xna = X
    Yna = Y
  }

  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data, u=Sigma_X_inv %*% GT$u,
              v=Sigma_Y_inv %*% GT$v, 
              Xna = Xna, Yna = Yna, newX = newX,
              newY=newY
  ))
}


