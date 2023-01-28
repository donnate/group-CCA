library(roxygen2)
library(pracma)
library(roperators)
library(gear)
library(tictoc)


cgd_solver = function(X,y, D, lambda1, lambda2, 
                      eps = 1e-4, max_it = 10000,
                      threshold = 1e-6,
                      verbose=TRUE){
  m = dim(D)[1]
  p = dim(D)[2]
  X_til = rbind(X, sqrt(2*lambda2)*D)
  y_til = rbind(y, matrix(0, m, 1))
  #tic("Pinv")
  X_til_pinv = pinv(X_til)
  #toc()
  #tic("Compute y_v")
  y_v = X_til %*% (X_til_pinv %*% y_til)
  #toc()
  #tic("Compute D_v")
  D_v = D %*% X_til_pinv
  #toc()
  #tic("Compute Q and b")
  Q = D_v %*% t(D_v)
  b = D_v %*% y_v
  #toc()

  u = matrix(0, m, 1)
  n_iter = 0
  prev_u = matrix(runif(m), m, 1)/m
  print("starting")
  print(paste0("lambda1 is ", lambda1, " lambda2 is ", lambda2))
  print(u)

  while(TRUE){
    if(verbose){
      tic(paste0("Iteration ", n_iter))
      }
    n_iter %+=% 1
    if(n_iter > max_it){
      if (verbose){print("Iterations exceed max_it")}
      return(X_til_pinv %*% (y_v - t(D_v) %*% u))
    }

    for (i in 1:m){
      t = 1/Q[i,i] * (b[i] - dot(Q[i,][-c(i)], u[-c(i)]))
      u[i] = sign(t)*min(abs(t), lambda1)
    }
    u[which(abs(u) < threshold )] = 0
    print(max(u), min(u))
    print(norm(u - prev_u, "2"))
    print(norm(u, '2'))
    if (norm(u) < threshold){
      break;
    }
    if ((norm(u - prev_u, '2')/(threshold + norm(prev_u, '2')) <= eps) | (norm(u, '2') < threshold) ){
      break
    }
    prev_u <- u
    if(verbose){
       toc()
       }
  }
  beta = X_til_pinv %*% (y_v - t(D_v) %*% u)
  return (beta)
}
