library(Matrix)

admm_cca <- function(X, Y, rho = 1.0, max_iter = 1000, tol = 1e-6) {
  # Initializations
  n <- dim(X)[1]
  d1 <- dim(X)[2]
  d2 <- dim(Y)[2]
  
  # These are the a and b vectors we want to estimate
  a <- matrix(0, d1, 1)
  b <- matrix(0, d2, 1)
  
  # ADMM variables
  lambda <- matrix(0, n, 1)
  u <- matrix(0, n, 1)
  
  for (i in 1:max_iter) {
    # a-update
    a <- solve(t(X) %*% X + rho * diag(d1)) %*% (t(X) %*% Y %*% b + rho * t(X) %*% (Y %*% b - u + lambda / rho))
    
    # b-update
    b <- solve(t(Y) %*% Y + rho * diag(d2)) %*% (t(Y) %*% X %*% a + rho * t(Y) %*% (X %*% a + u - lambda / rho))
    
    # u-update
    u <- X %*% a - Y %*% b + lambda / rho
    
    # Dual variable update
    lambda <- lambda + rho * (X %*% a - Y %*% b - u)
    
    # Convergence check
    primal_res <- norm(X %*% a - Y %*% b - u, type = "F")
    dual_res <- rho * norm(a - b, type = "F")
    if (primal_res < tol && dual_res < tol) {
      break
    }
  }
  
  return(list(a = a, b = b))
}

# Example
X <- matrix(rnorm(100*5), 100, 5)
Y <- matrix(rnorm(100*6), 100, 6)
result <- admm_cca(X, Y)
