# Function for Intialization via generalized fantope projection
# Inputs:
# =======
# A, B:       Pair of matrix to calculate generalized eigenspace (sigma, sigma0)
# nu:         Parameter of ADMM, default set to 1
# K:          nuclear norm constraint, equal to r
# rho:     penalty parameter on the l_1 norm of the solution, scaled by
#             sqrt(log(max(p1,p2))/n)
# epsilon:    tolerance level for convergence in ADMM
# maxiter:    maximum number of iterations in ADMM
# trace:      if set to True will print all iterations 

# Outputs:
# ========
# $Pi:     optimum of the convex program

fantope_init <-
  function(Sigma_xy, rho, r, eta=1, epsilon=5e-3, maxiter=1000, trace.logs=FALSE){
    criteria <- 1e10
    i <- 1
    # Initialize parameters
    G <- F <- oldF <-  diag(1,p,p)
    H <- matrix(0,p,p)
    # While loop for the iterations
    while(criteria > epsilon && i <= maxiter){
       # Update F
       F <- soft_threshold(G - 1/eta * H  + Sigma_xy, rho/eta)
       W <- 1/eta * H + F
       G <- updateG(W, r)
       H <- H +  (F - G) * eta	
       criteria <- sqrt(sum((F-oldF)^2))
       oldPi <- Pi
       i <- i+1
      if(trace.logs==TRUE)
      {
        print(i)
        print(criteria)
      }
    }
    return(list(Pi=Pi,H=H,Gamma=Gamma,iteration=i,convergence=criteria))
    
  }



soft_threshold <- function(A, lambda){
  return( A * ((abs(A) >lambda) * (1-lambda * sign(A))))
}


updateG <-
  function(W, r){
    #### We should convert everything to sparse matrices.
    temp <- (W+t(W))/2
    J = which(apply(temp, 2, sum) >0)
    temp_r = temp[J, J]
    svd_W <- svdr(temp, k = 3 * r)

    if(sum(pmin(1,pmax(svd_W$d,0)))<=K){
      dfinal <- pmin(1,pmax(svd_W$d,0))
      return(svd_W$u%*%diag(dfinal)%*%t(svd_W$v))
    }

    fr <- function(x){
      sum(pmin(1,pmax(d-x,0)))
    }
    # Vincent Vu Fantope Projection
    dfinal <- pmax(svd_W$d,0)
    gammas <-  dfinal[which(dfinal > 0)]
    choices <- sapply(1:length(gammas), function(i){
         sum(fr(gammas[1:i], gammas[i])) <= r})
    index <- argmax((1:length(gammas))[choices])
    dfinal <- pmin(1,pmax(gammas[1:index]-gammas[index], 0))
    res <- svd_W$u %*% dfinal %*% t(svd_W$v)
    return(res)
  }


RBKI <-(A, q, k){
    ### Generate Y
    ###  Needs to be appropriately implemented
    p <- ncol(A)
    Omega <- matrix(runif(n * r), ncol = r)
    for (i in 1:q){
        X = A %*% Y
        X = X - 
    }
    return(list(U=U, V=V, d=d))
} 

