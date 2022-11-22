library(fda)

D <- function(alpha){
  as.matrix(dist(alpha, diag = T, upper = T))
}

kernel <- function(alpha, lambda, gamma, A){
  W <- lambda * A * (gamma/pmax(D(alpha), 1e-12) + (1 - gamma))
  W[A < 1e-12] <- 0
  w <- rowSums(W)
  2 * (diag(w) - W)
}

rho <- function(X, Y, alpha, beta, K){
  alpha <- as.matrix(alpha)
  beta <- as.matrix(beta)
  Cxx = var(X, use = "pairwise") + K
  Cyy = var(Y, use = "pairwise")
  Cxy = cov(X, Y, use = "pairwise")
  c(t(alpha) %*% Cxy %*% beta/sqrt(t(alpha) %*% Cxx %*% alpha) / sqrt(t(beta) %*% Cyy %*% beta))
}

CCA <- function(X, Y, K){
  Cxx = var(X, use = "pairwise") + K
  Cyy = var(Y, use = "pairwise")
  Cxy = cov(X, Y, use = "pairwise")
  cc <- geigen(Cxy, Cxx, Cyy)
  names(cc) <- c("cor", "xcoef", "ycoef")
  list(alpha = cc$xcoef[,1], beta = cc$ycoef[,1], rho = cc$cor[1])
}


CCA_cluster <- function(X, Y, lambda, gamma = 1, A = 1, maxiter = 100){
  #lambda: penalty magnitude
  #gamma: if 1 then l1, if 0 then l2
  #A: weights in the penalty
  p <- ncol(X)
  q <- ncol(Y)
  iter <- 0
  alpha <- rnorm(p)
  beta <- rnorm(q)
  cc <- list(alpha = alpha, beta = beta)
  cor0 <-  rho(X, Y, cc$alpha, cc$beta, kernel(alpha, lambda, gamma, A))
  delta <- Inf
  cat("iter ", iter, "cor", cor0, "\n")
  while(delta > 1e-6 & iter < maxiter){
    iter <- iter + 1
    cc <- CCA(X, Y, kernel(cc$alpha, lambda, gamma, A))
    cor <-  rho(X, Y, cc$alpha, cc$beta, kernel(alpha, lambda, gamma, A))
    delta <- abs((cor - cor0)/cor0)
    cor0 <- cor
    cat("iter ", iter, "cor", cor, "delta", delta, "\n")
  }
  cc
}