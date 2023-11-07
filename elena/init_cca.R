library(pracma)

init_cca <- function(X,Y, r, penalize_X=TRUE,
                     penalize_Y=TRUE,
                     Zu=NULL, A_u=NULL,
                     Zv=NULL, A_v=NULL){
  Sigma_X <- cov(X)
  Sigma_Y <- cov(Y)
  Sigma_XY <- cov(X, Y)
  n = nrow(X)
  p = nrow(Sigma_X)
  q = nrow(Sigma_Y)
  u_init = matrix(0, p, r)
  v_init = matrix(0, q, r)
  Sigma_XY_new <- Sigma_XY
  new_X <- X
  new_Y <- Y
  ind.x <- NULL
  ind.y <- NULL
  if (penalize_X){
    ### Select the top rows
    if (is.null(A_u) == FALSE){
      #### Smooth on X
      new_X <- X %*% A_u
      
    }else{
      if(is.null(Zu) == FALSE){
        unique_Z_u = unique(Zu)
        new_X = matrix(0, n, length(unique_Z_u))
        for (i in 1:length(unique_Z_u)){
          new_X[,i] = apply(X[, which(Zu == unique_Z_u[i])],1,mean)
        }
      }else{
        norm_u <- (apply(Sigma_XY^2, 1, sum))
        #ind <- sort(norm_u, decreasing=TRUE, index.return = TRUE)
        sigma.hat <- mad(as.vector(norm_u))
        df = p
        alpha.theory = 1.5
        ind.x <- which(norm_u/sigma.hat^2/df > (1 + alpha.theory * sqrt(log(df)/df)))
        if (length(ind.x) < r) {
          ind.x <- (order(norm_u, decreasing = TRUE))[1:min(r + 10, p)]
        }
        print(ind.x)
        new_X <- X[, ind.x]
      }
    }
  }
  if (penalize_Y){
    if (is.null(A_v) == FALSE){
      #### Smooth on Y
      new_Y <- Y %*% A_v
    }else{
      if(is.null(Zv) == FALSE){
        unique_Z_v = unique(Zv)
        new_Y = matrix(0, n, length(unique_Z_v))
        for (i in 1:length(unique_Z_v)){
          new_Y[,i] = apply(Y[, which(Zv == unique_Z_v[i])],1,mean)
        }
      }else{
        norm_v <- (apply(Sigma_XY_new^2, 2, sum))
        sigma.hat <- mad(as.vector(norm_v))
        df = q
        alpha.theory = 1.5
        ind.y <- which(norm_v/sigma.hat^2/df > (1 + alpha.theory * sqrt(log(df)/df)))
        if (length(ind.y) < r) {
          ind.y <- (order(norm_v, decreasing = TRUE))[1:min(r + 10, q)]
        }
        print(ind.y)
        new_Y <- Y[, ind.y]
        
      }
    }
  }
  
  init_res = CCA::cc(new_X, new_Y)
  
  if (penalize_X){
    if (is.null(Zu) && is.null(A_u)){
      u_init[ind.x, ] <- init_res$xcoef[, 1:r]
    }else{
      for (i in 1:length(unique_Z_u)){
        u_init[which(Zu == unique_Z_u[i]), ] = init_res$xcoef[i, 1:r]
      }
    }
  }
  if (penalize_Y){
    if (is.null(Zv)&& is.null(A_v)){
      v_init[ind.y, ] <- init_res$ycoef[, 1:r]
    }else{
      for (i in 1:length(unique_Z_v)){
        v_init[which(Zv == unique_Z_v[i]), ] = init_res$ycoef[i, 1:r]
      }
    }
  }
  return(list(u_init = u_init,
              v_init = v_init))
}
