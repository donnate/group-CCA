library(MASS)
library(geigen)
library(pracma)
library(VGAM)
library(mvtnorm)
n= 500
p1=210
p2=230


example <- generate_other_example( n, p1, p2, a=0.4,  nnzeros = 10,
                                   theta = diag( c(1,  1)),
                                   r=2, signal_strength="normal")


example <- generate_example( n, p1, p2, nnzeros = 10,
                                   theta = diag( c(1,  1)),
                                   r=2, signal_strength="normal")

example <- generate_example_none_trivial_pca( n, p1, p2, nnzeros = 10,
                                            theta = diag( c(0.8,  0.9)), non_overlapping = TRUE,
                                            r=2, r_pca = 3, signal_strength="normal")
plot(example$X[,1], example$Y[,2])
### Understand what is going on
p_tot = p1+p2
GT = svd(example$Sigma[1:p1, (1+p1):(p1+p2)])


norm_Ys = sqrt(apply(example$Y^2, 2, sum))
norm_Ys = (1-(p2-2)/sum(norm_Ys^2)) * norm_Ys
norm_Xs = sqrt(apply(example$X^2, 2, sum))
norm_Xs = (1-(p1-2)/sum(norm_Xs^2))* norm_Xs
stat_X = apply((example$S[1:p1, (1+p1):(p_tot)])^2, 1, sum) * n
plot(stat_X)
quantile(stat_X, 0.8)
I = which(p.adjust(1-pchisq(stat_X, df=p2), "bonferroni") < 0.1)
stat_Y = apply((example$S[(1+p1):(p_tot), 1:p1])^2, 1, sum) * n 
plot(stat_Y)
J = which(p.adjust(1-pchisq(stat_Y, df=p1), "BH") < 0.1)
plot(stat_Y* n)

q = 100
S1 =  mvrnorm(q, mu =rep(0, p1), Sigma = diag(rep(1, p1)))
plot(apply((S1 %*% example$S[1:p1, (1+p1):(p_tot)])^2, 2, sum))

##### Simple Tests
GT = svd(sqrtm(example$Sigma[1:p1, 1:p1])$Binv %*% example$Sigma[1:p1, (p1+1):p_tot]%*% sqrtm(example$Sigma[(p1+1):(p_tot), (p1+1):(p1+1):(p_tot)])$Binv, 
         nu=2, nv=2)
res = svd(sqrtm(example$S[1:p1, 1:p1])$Binv %*% example$S[1:p1, (p1+1):p_tot] %*% sqrtm(example$S[(p1+1):(p_tot), (p1+1):(p1+1):(p_tot)])$Binv , nu=2, nv=2) #### This is not working out

hist(example$S[1:p1, (p1+1):p_tot])
res2 = svd(sqrtm(example$S[1:p1, 1:p1])$Binv %*% example$S[1:p1, (p1+1):p_tot] %*% sqrtm(example$S[(p1+1):(p_tot), 
                                                                                                   (p1+1):(p1+1):(p_tot)])$Binv , nu=2, nv=2) #### This is not working outmax()
subdistance(res2$u, GT$u)
subdistance(res2$v, GT$v)

res = svd(example$S[I, J], nu=2, nv=2) #### This is not working outmax()
u_test = matrix(0, p1, 2)
v_test = matrix(0, p2, 2)
u_test[I, ] = res$u
v_test[J, ] = res$v
u_test = u_test %*%(sqrtm(t(u_test[I,1:r]) %*% example$S[I, I] %*% u_test[I,1:r])$Binv) 
v_test = v_test %*%(sqrtm(t(v_test[J,1:r]) %*% example$S[(p1+J), (p1+J)] %*% v_test[J,1:r])$Binv) 
subdistance(u_test, GT$u)
subdistance(v_test, GT$v)

subdistance(res$v[, 1:2], example$v)
subdistance(res$u[, 1:2], example$u)
print(mean((example$Y %*% GT$v- example$Y %*% res$v)^2))
subdistance(example$Y %*% res$v[, 1:2], example$Y %*% GT$v)
##### Rnadomized version
q = 20
S1 = 1/sqrt(q) * matrix(rnorm(p1 * q), ncol = p1)
S2 = 1/sqrt(q) *  matrix(rnorm(p2 * q), ncol = p2)
B =  S1 %*% example$S[1:p1, (p1+1):p_tot] #%*% t(S2)
row_space = svd(B, nu=r, nv=r)
subdistance(row_space$v[, 1:2], example$v)
print(mean((example$Y %*% row_space$v- example$Y %*% GT$v)^2))
subdistance(example$Y %*% row_space$v[, 1:2], example$Y %*% GT$v)


res = ssvd(example$S[1:p1, (p1+1):p_tot] , r=2, method = "theory")
V_ssvd = res$v
V_ssvd = res$v %*%(sqrtm(t(res$v) %*% example$S[(p1+1):p_tot, (p1+1):p_tot] %*% res$v)$Binv)
U_ssvd = res$u
U_ssvd = res$u %*%(sqrtm(t(res$u) %*% example$S[1:p1, 1:p1] %*% res$u)$Binv)
print(mean((example$Y %*% res$v- example$Y %*% GT$v)^2))
#subdistance(V_ssvd[, 1:2], example$v)
subdistance(V_ssvd[, 1:2], GT$v)
subdistance(U_ssvd[, 1:2], GT$u)


res = my.ssvd(example$S[1:p1, (p1+1):p_tot],
              Sigma_u = example$S[1:p1, 1:p1],
              Sigma_v=example$S[(p1+1):p_tot, (p1+1):p_tot],
              r=2, method = "theory")
subdistance(res$v, GT$v)
subdistance(res$u, GT$u)
#subdistance(example$Y %*% res$v[, 1:2], example$Y %*% GT$v)


library(CAPIT)

resq = ssvd(example$S[1:p1, (p1+1):p_tot] , r=2, method = "method")
print(mean((example$Y %*% res$v- example$Y %*% GT$v)^2))
subdistance(resq$v[, 1:2], example$v)
subdistance(example$Y %*% res$v[, 1:2], example$Y %*% GT$v)

test <- CAPIT(example$X, example$Y)



estimate1 =matrix(0, p2, 2)
for (it in 1:100){
  S1 = 1/sqrt(q) * matrix(rnorm(p1 * q), ncol = p1)
  S2 = 1/sqrt(q) *  matrix(rnorm(p2 * q), ncol = p2)
  B =  S1 %*% example$S[1:p1, (p1+1):p_tot] #%*% t(S2)
  row_space = svd(B, nu=r, nv=r)
  subdistance(row_space$v[, 1:2], example$v)
  print(mean((example$Y %*% row_space$v- example$Y %*% res$v)^2))
  
  estimate1 = t(S2)%*% row_space$v[, 1:2]
  estimate1 = estimate1 %*% sqrtm(t(estimate1_m)%*% estimate1_m)$Binv
  estimate = estimate + estimate1
  estimate_m = estimate/it
  #estimate1_m = estimate1/it
  #estimate1 = estimate1 + t(S2)%*% row_space$v[, 1:2]
  #estimate1_m = estimate1/it
  #estimate = estimate1_m %*% sqrtm(t(estimate1_m)%*% estimate1_m)$Binv
  print(subdistance(estimate_m, example$v))
  print(mean((example$Y %*% GT$v - example$Y %*% estimate_m)^2))
}



#plot(example$Y %*% GT$v[,1], example$Y %*% estimate[,1])
#BV = B %*% row_space$v
#low_rank_approx = svd(BV, nu=r, nv=r)
#Vtilde = row_space$v %*% low_rank_approx$v
ag <- sgca_init(A= example$S, B=diag(1, p_tot), 
                rho=0.5*sqrt(log(p_tot)/n),K=r,
                nu=1,trace=FALSE)
ainit <- init_process(ag$Pi, r)
uest <- ainit
V =  uest[(p1+1):p_tot,]
V[which(abs(V) <1e-5)]=0
print(subdistance(GT$v,  V))

U =  uest[1:p1,]
U[which(abs(U) <1e-5)]=0
print(subdistance(GT$u,  U))
print(mean((example$Y %*% GT$v - example$Y %*% uest[(p1+1):p_tot,])^2))
subdistance(example$Y %*% V, example$Y %*% GT$v)


plot(as.numeric(example$Y %*% GT$v[,1]), 
      as.numeric(example$Y %*% uest[(p1+1):p_tot,1]))


# Load necessary libraries
library(Matrix)

# Define the original dimension d and the target dimension k
d <- p1
d2 <- p2
k <- q

# Step 1: Random Signs Matrix D
D <- Diagonal(x = sample(c(-1, 1), d, replace = TRUE))
D2 <- Diagonal(x = sample(c(-1, 1), d2, replace = TRUE))
# Step 2: Hadamard Transform
n <- nextn(d, 2)  # Find the smallest power of 2 greater than or equal to d
H <- hadamard(n)  # Create a Hadamard matrix of size n
H_d <- H[1:d, 1:d]  # Truncate to the original dimension d

n2 <- nextn(d2, 2)  # Find the smallest power of 2 greater than or equal to d
H2 <- hadamard(n2)  # Create a Hadamard matrix of size n
H2_d <- H2[1:d2, 1:d2]  # Truncate to the original dimension d

transformed_matrix <- D %*% H_d
transformed_matrix2 <- D2 %*% H2_d
library(CAPIT)



# Step 3: Projection
P <- matrix(0, nrow = k, ncol = d)
P2 <- matrix(0, nrow = k, ncol = d2)
for (i in 1:k) {
  row_idx <- sample(1:d, 1)
  P[i, row_idx] <- 1
  row_idx2 <- sample(1:d2, 1)
  P2[i, row_idx] <- 1
}
reduced_matrix <- P %*% transformed_matrix
reduced_matrix2 <- P2 %*% transformed_matrix2
# Your input data matrix A
A <- example$S[1:p1, (p1+1):p_tot]

# Apply FJLT to your data
reduced_data <- reduced_matrix %*% A %*% t(reduced_matrix2)

row_space = svd(reduced_data, nu=r, nv=r)
estimate1 = t(reduced_matrix2)%*% row_space$v[, 1:2]
estimate1 = estimate1 %*% sqrtm(as.matrix(t(estimate1)%*% estimate1))$Binv
#estimate = estimate + estimate1
estimate_m = estimate/it
#estimate1_m = estimate1/it
#estimate1 = estimate1 + t(S2)%*% row_space$v[, 1:2]
#estimate1_m = estimate1/it
#estimate = estimate1_m %*% sqrtm(t(estimate1_m)%*% estimate1_m)$Binv
print(subdistance(estimate1, example$v))
print(mean((example$Y %*% GT$v - example$Y %*% estimate1)^2))



plot(example$X[,1], example$Y[,1])
q = 20
r=2
p_tot = p1+p2
Yv =matrix(0, p2, r)
nb_it = 500
for (it in 1:nb_it){
  S1 =  mvrnorm(q, mu =rep(0, p1), Sigma = diag(rep(1/q, p1)))
  S2 =  mvrnorm(q, mu =rep(0, p2), Sigma = diag(rep(1/q, p2)))
  
  B =  S1 %*% example$S[1:p1, (p1+1):p_tot] %*% t(S2)
  row_space = svd(B, nu=r, nv=r)
  #BV = B %*% row_space$v
  #low_rank_approx = svd(BV, nu=r, nv=r)
  #Vtilde = row_space$v %*% low_rank_approx$v
  YV = YV + example$Y %*% t(S2) %*% row_space$v/sqrt(n) #nb_it
}

Yv =matrix(0, p2, r)
for(it in 1:100){
   S1 = matrix(rnorm(p1 * q), ncol = p1)
   S2 = matrix(rnorm(p2 * q), ncol = p2)
   B =  S1 %*% example$S[1:p1, (p1+1):p_tot] %*% t(S2)
   low_rank = svd(B, nu=2, nv=2)
   YV = YV +  example$Y%*% t(S2) %*% low_rank$v/(sqrt(n))
}
library(FastC)
test = fjlt(example$S[1:p1, (p1+1):p_tot] %, k = 100)

q = 100
S1 = matrix(rnorm(p1 * q), ncol = p1)
S2 = matrix(rnorm(p2 * q), ncol = p2)
B =  S1 %*% example$S[1:p1, (p1+1):p_tot] %*% t(S2)
low_rank = svd(B, nu=2, nv=2)
YV = example$Y%*% t(S2) %*% low_rank$v/(sqrt(n))/sqrt(q)
print(cbind(YV[1:12,]/100, YV_GT[1:12,]))

res = ssvd(example$S[1:p1, (p1+1):p_tot] , r=2, method = "theory") #### This is not working out
subdistance(res$v[, 1:2], example$v)
subdistance(res$u[, 1:2], example$u)

example$S[1:p1, (p1+1):p_tot]
AAT = example$S[1:p1, (p1+1):p_tot] %*% t(example$S[1:p1, (p1+1):p_tot])
plot(apply(AAT^2, 1, sum))

GT = svd(example$Sigma[1:p1, (p1+1):p_tot], nu=r, nv=r)
GT2 = svd(example$S[1:p1, (p1+1):p_tot], nu=r, nv=r)

YV_GT = example$Y %*% GT$v
YV_GT2 = example$Y %*% GT2$v
plot(YV_GT, YV)



ag <- sgca_init(A= example$S, B=example$sigma0hat, 
                rho=0.5*sqrt(log(p_tot)/n),K=r,
                nu=1,trace=FALSE)
ainit <- init_process(ag$Pi, r)
uest <- ainit
uest[which(abs(uest)<1e-5)]=0
subdistance(uest[1:p1, 1:2], example$v)
subdistance(uest[(1+p1):p_tot, 1:2], example$u)



#' @export
rsvd_product <- function(A, k=NULL, nu=NULL, nv=NULL, p=10, q=2, sdist="normal") {
  #*************************************************************************
  #***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
  #***                              <2015>                               ***
  #***                       License: BSD 3 clause                       ***
  #*************************************************************************
  
  #Dim of input matrix
  n <- nrow(A)
  m <- ncol(A)
  
  #Flip matrix, if wide
  if(m < n){
    A <- H(A)
    n <- nrow(A)
    m <- ncol(A)
    flipped <- TRUE
  } else flipped <- FALSE
  
  #Set target rank
  if(is.null(k)) k = n
  if(k > n) k <- n
  if(is.character(k)) stop("Target rank is not valid!")
  if(k < 1) stop("Target rank is not valid!")
  
  #Set oversampling parameter
  l <- round(k) + round(p)
  if(l > n) l <- n
  if(l < 1) stop("Target rank is not valid!")
  
  #Check if array is real or complex
  if(is.complex(A)) {
    isreal <- FALSE
  } else {
    isreal <- TRUE
  }
  
  #Set number of singular vectors
  if(is.null(nu)) nu <- k
  if(is.null(nv)) nv <- k
  if(nu < 0) nu <- 0
  if(nv < 0) nv <- 0
  if(nu > k) nu <- k
  if(nv > k) nv <- k
  if(flipped==TRUE) {
    temp <- nu
    nu <- nv
    nv <- temp
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Generate a random sampling matrix O
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  O <- switch(sdist,
              normal = matrix(stats::rnorm(l*n), n, l),
              unif = matrix(stats::runif(l*n), n, l),
              rademacher = matrix(sample(c(-1,1), (l*n), replace = TRUE, prob = c(0.5,0.5)), n, l),
              stop("Selected sampling distribution is not supported!"))
  
  if(isreal==FALSE) {
    O <- O + switch(sdist,
                    normal = 1i * matrix(stats::rnorm(l*n), n, l),
                    unif = 1i * matrix(stats::runif(l*n), n, l),
                    rademacher = 1i * matrix(sample(c(-1,1), (l*n), replace = TRUE, prob = c(0.5,0.5)), n, l),
                    stop("Selected sampling distribution is not supported!"))
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Build sample matrix Y : Y = O * A 
  #Note: Y should approximate the range of A
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Y <- O %*% A
  remove(O)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Orthogonalize Y using economic QR decomposition: Y=QR
  #If q > 0 perfrom q subspace iterations
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if( q > 0 ) {
    for( i in 1:q) {
      Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
      Z <- crossprod_help(A , Y )
      Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
      Y <- A %*% Z
    }#End for
    remove(Z)
  }#End if
  
  Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
  remove(Y)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Project the data matrix a into a lower dimensional subspace
  #B := Q.T * A
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  B <- crossprod_help(Q , A )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Singular Value Decomposition
  #Note: B =: U * S * Vt
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rsvdObj <- svd(B, nu=nu, nv=nv) # Compute SVD
  rsvdObj$d <- rsvdObj$d[1:k] # Truncate singular values
  
  if(nu != 0) rsvdObj$u <- Q %*% rsvdObj$u # Recover left singular vectors
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Flip SVD back
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(nu == 0){ rsvdObj$u <- NULL}
  if(nv == 0){ rsvdObj$v <- NULL}
  if(flipped == TRUE) {
    u_temp <- rsvdObj$u
    rsvdObj$u <- rsvdObj$v
    rsvdObj$v <- u_temp
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Return
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  class(rsvdObj) <- "rsvd"
  return(rsvdObj) 
  
} # End rsvd


H <- function( X ) {
  if(is.complex(X)) {
    return( Conj(t(X)) )
  } else {
    return( t(X) )
  }
}

#Helper function for conjugate crossprod
#' @importFrom Matrix crossprod
crossprod_help <- function( A , B ) {
  if(is.complex(A)) {
    return( crossprod( Conj(A) , B) )
  } else {
    return( crossprod( A , B ) )
  }
}

#Helper function for conjugate tcrossprod
#' @importFrom Matrix tcrossprod
tcrossprod_help <- function( A , B ) {
  if(is.complex(B)) {
    return( tcrossprod( A , Conj(B) ) )
  } else {
    return( tcrossprod( A , B ) )
  }
}

#Helper function for Moore Penrose pseudoinverse
pinv <- function(A){
  s <- svd(A)
  nz <- s$d > s$d[1] * .Machine$double.eps
  if(any(nz)){
    return(s$v[, nz] %*% (H(s$u[, nz]) / s$d[nz]))
  } else {
    return(A)
  }
}

