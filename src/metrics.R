#########################################
# AUXILIARLY FUNCTIONS SIMULATION STUDY #
#########################################


principal_angles <- function(a, b){
  ### Calculate principal angles between subspace spanned by the columns of a and the subspace spanned by the columns of b

  angles=matrix(0, ncol=ncol(a), nrow=1)
  qa= qr.Q(qr(a))
  qb= qr.Q(qr(b))
  C=svd(t(qa)%*%qb)$d
  rkA=qr(a)$rank;rkB=qr(b)$rank
  if (rkA<=rkB){
    B = qb - qa%*%(t(qa)%*%qb);
  } else {B = qa - qb%*%(t(qb)%*%qa)}
  S=svd(B)$d
  S=sort(S)

  for (i in 1:min(rkA, rkB)){
    if (C[i]^2 < 0.5) {angles[1, i]=acos(C[i])}
    else if (S[i]^2 <=0.5) {angles[1, i]=asin(S[i])}
  }
  angles=t(angles)

  ##OUTPPUT
  out=list(angles=angles)
}

TPR  <-  function(A, B){
  # This is a function that compares the structure of two matrixes A and B
  # It outputs the number of entries that A and B have in common that are different from zero
  # A and B need to have the same number of rows and columns
  out  <-  sum(A!=0&B!=0)/sum(A!=0)
}

TNR  <-  function(A, B){
  # This is a function that compares the structure of two matrixes A and B
  # It outputs the number of entries that A and B have in common that are zero #
  # A and B need to have the same number of rows and columns
  out  <-  sum(A==0&B==0)/sum(A==0)
}

l1_smoothness  <-  function(A, D){
  out  <-  sum( abs(D %*% A))
}

l2_smoothness  <-  function(A, D){
  out  <-  sum( (D %*% A)^2)
}


evaluate_method <- function(xcoef, ycoef, 
                            Sigma_x, Sigma_y,
                            X_train, Y_train,
                            X_test, Y_test,
                            D,
                            normalize=TRUE){
  if (normalize){
    svd_sol = svd(t(xcoef) %*% Sigma_x %*% xcoef )
    inv_sqrt_sol = svd_sol$u %*% diag(sapply(svd_sol$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol$v)
    svd_sol2 = svd(t(ycoef) %*% Sigma_y %*% ycoef )
    inv_sqrt_sol2 = svd_sol2$u %*% diag(sapply(svd_sol2$d, function(x){ifelse(x<1e-5, 0, 1/sqrt(x))})) %*% t(svd_sol2$v)
    results.x <- xcoef %*%  inv_sqrt_sol
    results.y <- ycoef %*% inv_sqrt_sol2
  }else{
    results.x <- xcoef 
    results.y <- ycoef
  }
  MSEa <- principal_angles(trueA, xcoef)$angles[1,1]
  MSEb <- principal_angles(trueB, ycoef)$angles[1,1]
  TPRa <- TPR(trueA, xcoef)
  TPRb <- TPR(trueB, ycoef)
  TNRa <- TNR(trueA, xcoef)
  TNRb <- TNR(trueB, ycoef)
  l1_smooth <- l1_smoothness(xcoef, D)
  l2_smooth <- l2_smoothness(xcoef, D)
  l2loss  <- sqrt(mean((X_train %*% results.x - Y_train %*% results.y )^2))
  l2loss_test <- sqrt(mean((X_test %*% results.x - Y_test %*% results.y )^2))
  correl = cor(trueA[,1], xcoef[,1])
  
  return(list(MSEa =MSEa , MSEb= MSEb, TPRa=TPRa, TPRb=TPRb,
              TNRa=TNRa, TNRb=TNRb, l2loss=l2loss, l2loss_test=l2loss_test,
              l1_smooth=l1_smooth, l2_smooth=l2_smooth,
              correl=correl
              ))

}