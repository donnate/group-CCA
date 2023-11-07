library(Matrix)
library(CVXR)

get.subset <-function (x, method = c("theory", "method"), alpha.method = 0.05,
            alpha.theory = 1.5, sigma = NA, df = Inf) {
    if (method == "theory") {
      ans <- which(x/sigma^2/df > (1 + alpha.theory * sqrt(log(df)/df)))
    }
    if (method == "method") {
      ans <- holm.robust(x = x, alpha = alpha.method)
    }
    ans
}

init_cca <- function(X, Y, r, lambdax, lambday, Gamma_u, Gamma_v){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  if (lambdax >0  & is.null(Gamma_u)){
    rownorm2 <- apply(atanh(cov(X, Y))^2, 1, sum) *n
    all_I_p = p.adjust(1-pchisq(rownorm2, q))
    I0 = which(all_I_p < 0.1)
  }else{
    I0  = 1:p
    if (is.null(Gamma_u) == FALSE & lambdax >0){
      L = t(Gamma_u) %*% Gamma_u
      L <- L -diag(diag(L))
      A = -L 
      X = A %*% X
    }
  }
  
  if (lambday >0  & is.null(Gamma_v)){
    colnorm2 <- apply(atanh(cov(X, Y))^2, 2, sum) * n
    all_J_p = p.adjust(1-pchisq(colnorm2, p))
    J0 = which(all_J_p < 0.1)
  }else{
    J0  = 1:q
    if (is.null(Gamma_v) == FALSE & lambday >0){
      L = t(Gamma_v) %*% Gamma_v
      L <- L -diag(diag(L))
      A = -L 
      Y = A %*% Y
    }
  }
  
  if(length(I0)<2){
    I0 = 1:p
  }
  if(length(J0)<2){
    J0 = 1:p
  }
  
  init = CCA::cc(X[, I0], Y[, J0])
  U = matrix(0, p, r)
  U[I0, ] <- init$xcoef[,1:r]
  V = matrix(0, q, r)
  V[J0, ] <- init$ycoef[,1:r]
  return(list(U=U, V=V))

}

alternating_cca <- function(X, Y, r, init_coef = NULL, lambdax = 0,
                            lambday=0, Gamma_u=NULL, Gamma_v=NULL,
                            thres=1e-5, max_iter=1000){
  # Initializations
  n <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Y)[2]
  
  # Check for missingness
  missY = is.na(Y)
  missX = is.na(X)
  
  if(max(missX) >0){
    X = data.frame(X) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))) %>% as.matrix()
  } 
  if(max(missY) >0){
    Y = data.frame(Y) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))) %>% as.matrix()
  } 
  
  if (is.null(init_coef)){
    init_coef <- init_cca(X, Y, r, lambdax, lambday, Gamma_u, Gamma_v)
  }
  U <- init_coef$U
  V <- init_coef$V
  ZV <- Y %*% V
  ZU <- X %*% U
  U.old = U
  V.old = V
  U_init = U
  V_init = V
  it = 1
  eps = 1e10
  while(it<max_iter & eps > thres) {
    if (lambdax >0){
      if (is.null(U_init)){
         U <- Variable(p, r)
      }else{
         U <- Variable(p, r)
      }
     
      if (is.null(Gamma_u)){
        objective <- Minimize(1/n * sum_squares(X %*% U - ZV) + lambdax * sum(norm2(U, axis=1)))
        problem <- Problem(objective)
        result <- solve(problem)
        U <- result$getValue(U)
        U_init <- U
      }else{
        objective <- Minimize(1/n * sum_squares(X %*% U - ZV) + lambdax * sum(norm2(Gamma_u %*% U, axis=1)))
        problem <- Problem(objective)
        result <- solve(problem)
        U <- result$getValue(U)
        U_init <- U
      }

    }else{
      U = solve(t(X) %*%X) %*% t(X) %*% ZV
    }
    #U = U %*% diag(1/sqrt(apply(U^2, 2, sum)))
    if (norm(U) < 1e-5){
      ### Solution is essentially 0
      return(list(U=U, V = matrix(0, n, q)))
    }
    U = U %*% sqrtm(t(U) %*% cov(X) %*% U)$Binv
    ZU <- X %*% U

    if (lambday >0){
      if (is.null(V_init)){
         V <- Variable(p, r)
      }else{
         V <- Variable(p, r)
      }
      if (is.null(Gamma_v)){
        objective <- Minimize(1/n * sum_squares(ZU - Y %*% V) + lambday * sum(norm2(V, axis=1)))
        problem <- Problem(objective)
        result <- solve(problem)
        V <- result$getValue(V)
        V_init <- V
      }else{
        objective <- Minimize(1/n * sum_squares(ZU - Y %*% V) + lambday * sum(norm2(Gamma_v %*% V, axis=1)))
        problem <- Problem(objective)
        result <- solve(problem)
        V <- result$getValue(V)
        V_init <- V
      }

    }else{
      V = solve(t(Y) %*%Y) %*% t(Y) %*% ZU
    }
    #V = V %*% diag(1/sqrt(apply(V^2, 2, sum)))
    if (norm(V) < 1e-5){
      ### Solution is essentially 0
      return(list(U=U, V = V))
    }
    V = V %*% sqrtm(t(V) %*% cov(Y) %*% V)$Binv
    #### Solve for Y and X
    if(sum(missY)>0){
     Y[missY] <- ((X %*% U) %*%  t(V) %*% pinv((V) %*% t(V)))[missY]
    }
    if(sum(missX)>0){
     X[missX] <- ((Y %*% V) %*%  t(U)  %*% pinv((U) %*% t(U)))[missX]
    }
    ZV <- Y %*% V
    eps = subdistance( U, U.old) + subdistance(V, V.old)
    U.old = U
    V.old = V
    it = it + 1
    print(it)
  }
  print(it)
  return(list(U=U, V = V))
}
