library(Matrix)



iterative_CCA <- function(X, Y, k){
  
  # Initializations
  n <- dim(X)[1]
  p1 <- dim(X)[2]
  p2 <- dim(Y)[2]
  
  # Center the variables
  X = scale(X, center = TRUE, scale=FALSE)
  Y = scale(Y, center = TRUE, scale=FALSE)
  
  # These are the U and V matrices we want to estimate
  U <- matrix(0, p1, k)
  V <- matrix(0, p2, k)
  ZU <- matrix(0, p1, k)
  ZV <- matrix(0, p2, k)
  
  # ADMM variables
  Lambda1 <- matrix(0, k, k)
  Lambda2 <- matrix(0, k, k)
  
  I <- diag(k)
  
  for (i in 1:max_iter) {
    # U-update
    term1 <- 1/n * t(X) %*% Y %*% V +  1/n * t(X) %*% (rho1 * t(ZU)  + Lambda1)
    #term1 <- threshold()
    term2 <- 1/n * t(X) %*% X + 1/n * rho1 * I
    U <- solve(term2, term1)
    
    ZU = (Lambda1 + X %*% U)
    ZU = ZU %*% sqrtm(t(ZU) %*% ZU)
    
    # Lambda1-update
    Lambda1 <- Lambda1 + rho1/n * (X %*% U - ZU)
    
    
    
    # V-update
    term2 <- 1/n * t(Y) %*% X %*% U +  1/n * t(Y) %*% (rho2 * t(ZV)  + Lambda2)
    #term1 <- threshold()
    term2 <- 1/n * t(Y) %*% Y + 1/n * rho2 * I
    V <- solve(term2, term1)
    
    ZV = (Lambda2 + Y %*% V)
    ZV = ZV %*% sqrtm(t(ZV) %*% ZV)
    
    # Lambda2-update
    Lambda2 <- Lambda2 + rho2/n * (Y %*% V - ZV)
    
    
    # Convergence check
    primal_res_U <- norm(t(U) %*% X %*% X %*% U - I, type = "F")
    primal_res_V <- norm(t(V) %*% Y %*% Y %*% V - I, type = "F")
    
    if (primal_res_U < tol && primal_res_V < tol) {
      break
    }
  }
  
  return(list(U = U, V = V))
}
