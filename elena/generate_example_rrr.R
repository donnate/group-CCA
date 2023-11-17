library(Matrix)
library(pracma)


generate_example_sparse_U <- function(n, p1, p2,
                                      r_pca = 3,
                                      nnzeros = 5,
                                      theta = diag(c(0.9,  0.8)),
                                      lambda_pca = 1,
                                      r = 2, overlapping_amount = 0,
                                      normalize_diagonal = TRUE,
                                      n_new = 50000) {
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
  s_pca  <- (start + 1) : (start + nnzeros)
  s_pca2  <- 1:q
  Lambda_pca <- rep(lambda_pca, r_pca)
  # generate vcovariance matrix for X and Y
  
  if (r_pca > 0){
    u1 = matrix(0, p1, r_pca)
    u1[s_pca, ] <- matrix(runif(n = nnzeros * r_pca, max = 3, min=0), nrow = nnzeros, ncol = r_pca)
    # Normalize u1
    u1[s_pca,] <- 0.5 * diag(1/ sqrt(apply(u1[s_pca, ]^2, 1, sum))) %*% u1[s_pca,]
    u1[s_pca, 1:r_pca] <- u1[s_pca,1:r_pca] #%*% (sqrtm(t(u1[s_pca,1:r_pca]) %*% u1[s_pca, 1:r_pca])$Binv)
    T1 <- u1 %*% diag(Lambda_pca) %*% t(u1)
    if (normalize_diagonal){
      diag(T1) <- 1
      Sigma[1:p1, 1:p1] <- T1
    }else{
      Sigma[1:p1, 1:p1] <- Sigma[1:p1, 1:p1] + T1 
    } 
    ### Same for Sigma_y
    u2 <- matrix(runif( p2 * r_pca, max = 1, min=0), nrow=p2)
    u2 <- 0.5 * diag(1/ sqrt(apply(u2^2, 1, sum))) %*% u2
    #u2[s_pca2, 1:r_pca] <- u2[s_pca2, 1:r_pca] #%*% (sqrtm(t(u2[s_pca2, 1:r_pca]) %*% u2[s_pca2, 1:r_pca])$Binv)
    T2 <- u2 %*% t(u2)
    if (normalize_diagonal){
      diag(T2) <- 1
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2
    }else{
      Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] <- T2 + Sigma[(p1 + 1):(p1 + p2), (p1 + 1) : (p1 + p2)]
    }
  }
  Sigmax <- Sigma[1:p1,1:p1];
  Sigmay <- Sigma[(p1+1):p_tot,(p1+1):p_tot];
  
  ### Generate cross covariance
  precision = solve(Sigmax)
  Tss <- precision[s,s]
  prod <- matrix(0, pp[1], r)
  prod[s, 1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  prod <- prod %*% (sqrtm(t(prod[s, 1:r]) %*% Tss %*% prod[s, 1:r])$Binv)
  u = precision %*% prod
  
  precision_y = solve(Sigma[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)] )
  Tss_v <- precision_y
  prod_v <- matrix(0, pp[2], r)
  prod_v[, 1:r] <- as.matrix(runif( p2 * r , max = 3, min=1), nrow=p2)  * as.matrix(sample(c(-1,1), p2 * r, replace=TRUE), nrow=p2)
  prod_v <- prod_v %*% (sqrtm(t(prod_v[, 1:r]) %*% Tss_v %*% prod_v[, 1:r])$Binv)
  v = precision_y %*% prod_v
  
  #Sigma[(p1 + 1) :(p1 + p2), 1:p1] <- Sigmax %*% u %*% theta %*% t(v) %*% Sigmay
  Sigma[1:p1, (p1 + 1):(p1 + p2)] <- prod %*% theta %*% t(prod_v)
  Sigma[(p1 + 1) :(p1 + p2), 1:p1]  <- prod_v%*% theta %*% t(prod)
  
  #Generate Multivariate Normal Data According to Sigma
  sqrt_Sigma <- sqrtm(Sigma)$B
  Data <- matrix(rnorm(n = (n + n_new) * p_tot), ncol=p_tot) %*% sqrt_Sigma
  X <- Data[1:n, 1:p1]
  Y <- Data[1:n, (p1 + 1):(p1 + p2)]
  Xnew <- Data[(n+1):(n_new + n), 1:p1]
  Ynew <- Data[(n+1):(n_new + n), (p1 + 1):(p1 + p2)]
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
  #Sigma_X_inv <- sqrtm(Sigma[1:p1, 1:p1])$Binv
  #Sigma_Y_inv <-  sqrtm(Sigma[(p1+1):(p_tot), (p1+1):(p_tot)])$Binv
  #GT = svd(Sigma_X_inv %*% Sigma[1:p1, (p1 + 1):p_tot] %*% Sigma_Y_inv, nu = r, nv = r)
  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data, 
              u=u,#sqrtm(Sigma[1:p1, 1:p1])$Binv %*% GT$u, 
              v=v, #sqrtm(Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)])$Binv %*% GT$v, 
              Xnew = Xnew, Ynew=Ynew,
              Sigmax=Sigmax, Sigmay=Sigmay
  ))
}