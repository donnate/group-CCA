############## evaluation #################

principal_angles = function(A, B){
  angles = rep(0, ncol(A))
  QRa = qr(A)
  QRb = qr(B)
  Qa = qr.Q(QRa)
  Qb = qr.Q(QRb)
  
  d = svd(t(Qa) %*% Qb)$d
  if(QRa$rank <= QRb$rank) Q = Qb - Qa %*% (t(Qa) %*% Qb)
  else Q = Qa - Qb %*% (t(Qb) %*% Qa)
  s = sort(svd(Q)$d)
  
  for(i in 1:min(QRa$rank, QRb$rank)){
    if(d[i]^2 < 0.5) angles[i] = acos(d[i])
    else if(s[i]^2 <= 0.5) angles[i] = asin(s[i])
  }
  angles
}

angles = function(A, B){
  p = ncol(A)
  angs = rep(0, p)
  for(i in 1:p){
    angs[i] = principal_angles(A[,1:i,drop = F], B[,1:i,drop = F])[1]
  } 
  angs
}

subdistance <- function(A, B){
  svdresult = svd(t(A) %*% B)
  U = svdresult$u
  V = svdresult$v
  O = U %*% t(V);
  l = norm(A %*% O-B, type = c('F'))
  return(l)
}

sinTheta<- function(U, V){
  l = 1/sqrt(2) * norm(U %*% t(U)-V %*% t(V), type = c('F'));
  return(l)
}


TPR  <-  function(A, B, tol=1e-4){
  # This is a function that compares the structure of two matrices A and B
  # It outputs the number of entries that A and B have in common that are different from zero
  # A and B need to have the same number of rows and columns
  A[which(abs(A) <tol)] =0
  B[which(abs(B) <tol)] =0
  out  <-  sum((A!=0)*(B!=0))/max(1,sum(A!=0))
  return(out)
}

FPR  <-  function(A, B, tol=1e-4){
  # This is a function that compares the structure of two matrices A and B
  # It outputs the number of entries that A and B have in common that are different from zero
  # A and B need to have the same number of rows and columns
  A[which(abs(A) <tol)] =0
  B[which(abs(B) <tol)] =0
  out  <-  sum((A!=0)*(B==0))/max(1,sum(A!=0))
  return(out)
}

TNR  <-  function(A, B, tol=1e-4){
  # This is a function that compares the structure of two matrices A and B
  # It outputs the number of entries that A and B have in common that are zero #
  # A and B need to have the same number of rows and columns
  A[which(abs(A) <tol)] =0
  B[which(abs(B) <tol)] =0
  out  <-  sum((A==0)*(B==0))/max(1,sum(A==0))
  return(out)
}

evaluate = function(X, Y, U, V, U0, V0, 
                    Sigma_hat_sqrt, Sigma0_sqrt,
                    thres = 1e-4){
  U_tot = rbind(U,V)
  U0_tot = rbind(U0,V0)
  p1 = ncol(X)
  p2 = ncol(Y)
  p = p1 + p2
  n_new = nrow(X)
  r = ncol(U0)
  Y = scale(Y, scale = FALSE)
  X = scale(X, scale = FALSE)
  silly_benchmark = subdistance(matrix(0, p1 + p2, r), U0_tot)
  all_cors = diag(cor(as.matrix(X) %*% U, as.matrix(Y) %*% V))
  data.frame(
    "n_new" = n_new,
    #"exp" = seed,
    "p1" = p1,
    "p2" = p2,
    "r" = r,
    "zero_benchmark" = silly_benchmark,
    "nb_discoveries" = sum(apply(U_tot^2, 1, sum)>0),
    "nb_real_discoveries" = sum(apply(U_tot^2, 1, sum)>thres),
    "distance_tot" = subdistance(U_tot, U0_tot),
    "distance_U" = subdistance(U, U0),
    "distance_V" = subdistance(V, V0),
    "sinTheta_tot" = sinTheta(Sigma_hat_sqrt %*% U_tot, Sigma0_sqrt %*%U0_tot),
    "sinTheta_U" = sinTheta(Sigma_hat_sqrt[1:p1, 1:p1] %*% U, 
                            Sigma0_sqrt[1:p1, 1:p1] %*% U0),
    "sinTheta_V" = sinTheta(Sigma_hat_sqrt[(p1+1):p, (p1+1):p] %*%V, 
                            Sigma0_sqrt[(p1+1):p, (p1+1):p] %*% V0),
    "prediction_tot" = mean((X %*% U - Y %*% V)^2),
    "prediction_U" = mean((X %*% U - X %*% U0)^2),
    "prediction_V" = mean((Y %*% V - Y %*% V0)^2),
    "prediction_U_sub" = subdistance(X %*% U, X %*% U0)^2,  ### need to be rotated
    "prediction_V_sub" = subdistance(Y %*% V, Y %*% V0)^2,
    "avg_corr" = mean(all_cors),
    "max_corr" = max(all_cors),
    "min_corr" = min(all_cors),
    "med_corr" = median(all_cors),
    "TPR" =TPR(apply(U_tot^2, 1, sum), apply(U0_tot^2, 1, sum), tol=thres),
    "TNR" = TNR(apply(U_tot^2, 1, sum), apply(U0_tot^2, 1, sum), tol=thres),
    "FPR" = FPR(apply(U_tot^2, 1, sum), apply(U0_tot^2, 1, sum), tol=thres),
    "FNR" = FPR(apply(U0_tot^2, 1, sum),apply(U_tot^2, 1, sum), tol=thres)
  )
}

evaluate_nogold = function(X, Y, U, V){
  XU = X %*% U
  YV = Y %*% V
  data.frame(comp = 1:ncol(V), 
             cors = diag(cor(XU, YV)),
             mses = colMeans((XU - YV)^2),
             angles = angles(XU, YV),
             cors_perm = compmatch(XU, YV, metric = "cor"),
             mses_perm = compmatch(XU, YV, metric = "mse"))
  
}

compmatch = function(A, B, metric){
  K = ncol(A)
  n = nrow(A)
  if(metric == "cor"){
    pwerrors = cor(A, B)
    perm = solve_LSAP(pwerrors - min(pwerrors), maximum = T)
  }
  if(metric == "mse"){
    pwerrors = (outer(diag(t(A) %*% A), rep(1, K)) + outer(rep(1, K), diag(t(B) %*% B)) - 2 * t(A) %*% B)/n
    perm = solve_LSAP(pwerrors, maximum = F)
  } 
  #print(pwerrors)
  return(list(permutation = perm, errors = pwerrors[cbind(1:K, perm)]))
}

############## graphics #################

gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

