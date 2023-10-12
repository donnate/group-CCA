library(MASS)
library(stats)
library(geigen)
library(pracma)
library(expm)
library(ramify)
setwd("~/Documents/Sparse-Generalized-Correlation-Analysis/src/R Code/")
source("sgca_init.R")
source("gca_to_cca.R")
source("init_process.R")
source("sgca_tgd.R")
source("subdistance.R")
source("utils.R")
#Example of TGD on Sparse CCA
#n = 500, p1 = p2 = 100, s_u = s_v = 5
#k = 20, eta = 0.0025, lambda =0.01, T = 12000


print('--------------------------------------');
print('SCCA Example: Toeplitz covariance matrix'); 
print('n = 500, p1 = p2 = 100, s_u = s_v = 5');
print('lambda1 = 0.9, lambda2 = 0.8');

n  <- 500;
p1 <- 1000;
p2 <- 1000;
p <- p1 + p2;
pp <- c(p1,p2);
theta <- diag( c(0.9,  0.8) );
r <- 2
print('--------------------------------------');
print('Generating data ...');
a <- 0.3;
Sigma <- diag(p1+p2)


# generate covariance matrix for X and Y
s  <- sample(1:min(p1,p2),5);
s2  <- sample(1:min(p1,p2),5);
T1 = toeplitz(a^(0:(pp[1]-1)));
Sigma[1:p1, 1:p1] = T1;
Tss = T1[s, s];
u = matrix(0, pp[1], r)
u[s,1:r] <- matrix(floor(runif(min=1, max=5, n = length(s) * r)) - 3, ncol=r)
u <- u %*%(pracma::sqrtm(t(u[s,1:r]) %*% Tss %*% u[s,1:r])$Binv)

T2 = toeplitz(a^(0:(pp[2]-1)));
Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)] = T2;
Tss = T2[s2, s2];
v = matrix(0, pp[2], r)
v[s2,1:r] <-  matrix(floor(runif(min=1, max=5, n = length(s2) * r)) - 3, ncol=r)
v <- v %*%(pracma::sqrtm(t(v[s2,1:r]) %*% Tss %*% v[s2,1:r])$Binv)

Sigma[(p1+1):(p1+p2), 1:p1] = T2 %*%  v  %*% theta %*% t(u) %*% T1;
Sigma[1:p1, (p1+1):(p1+p2)] = t(Sigma[(p1+1):(p1+p2), 1:p1])

Sigmax = Sigma[1:p1,1:p1];
Sigmay = Sigma[(p1+1):p,(p1+1):p];

#Generate Multivariate Normal Data According to Sigma
Data = mvrnorm(n, rep(0, p), Sigma);

X = Data[,1:p1];
Y = Data[,(p1+1):(p1+p2)];
#z = matrix(rnorm(n * r), nrow=n)
#epsilon = matrix(rnorm(n * p1), nrow=n)
#X = z %*% t(u)   + epsilon

C = cov(X)
col = rep("black", p1)
subset_x = which(apply(u, 1, function(x){sum(x^2)})>0)
subset_y = which(apply(v, 1, function(x){sum(x^2)})>0)
col[subset_x] = "red"
plot(diag(C), col=col)

C = cov(X, Y)
plot(apply(C, 1, function(x){sum(x^2)}), col=col)

plot(apply(C, 1, function(x){max(abs(x))}), col=col)

hist(apply(C, 1, function(x){max(abs(x))}), bins=30
     )

library(tidyverse)
ggplot(data.frame(x=(apply(C, 1, function(x){max(abs(x))}))),
                  aes(x=x)) + geom_histogram(
                  )
library(limma)
alpha=0.70
subset_x_hat = which(apply(C, 1, function(x){max(abs(x))}) > (1+ alpha) * sqrt(log(p1)/n))
subset_y_hat = which(apply(C, 2, function(x){max(abs(x))}) > (1+ alpha) * sqrt(log(p2)/n))
#### Selet them out, then:
XtY = cov(X[ ,subset_x_hat], Y[, subset_y_hat])
init = CCA::cc(X[ ,subset_x_hat],  Y[, subset_y_hat])
test_cca = CCA::cc(X[ ,subset_x],  Y[, subset_y])

init$xcoef[,1:2]
U_mine = matrix(0, p1, r)
U_mine[subset_x_hat, ] = init$xcoef[,1:2]
V_mine = matrix(0, p2, r)
V_mine[subset_y_hat, ] = init$ycoef[,1:2]

U_test = matrix(0, p1, r)
U_test[subset_x, ] = test_cca$xcoef[,1:2]
V_test = matrix(0, p2, r)
V_test[subset_y, ] = test_cca$ycoef[,1:2]

sqx <- pracma::sqrtm(Sigmax)$B
sqy <- pracma::sqrtm(Sigmay)$B

print('Initial prediction error of U is')
print( subdistance(sqx %*% U_mine, sqx %*% u)^2)
print('Initial prediction error of V is')
print( subdistance(sqy %*% V_mine, sqy %*% v)^2)

print('Test prediction error of U is')
print( subdistance(sqx %*% U_test, sqx %*% u)^2)
print('Test prediction error of V is')
print( subdistance(sqy %*% V_test, sqy %*% v)^2)
print('Data generated.');
print('--------------------------------------');

Mask = matrix(0, p, p);
idx1 = 1:pp[1];
idx2 = (pp[1]+1):(pp[1]+pp[2]);
Mask[idx1,idx1] <- matrix(1,pp[1],pp[1]);
Mask[idx2,idx2] <- matrix(1,pp[2],pp[2]);
Sigma0 = Sigma * Mask;

S <- cov(Data)
sigma0hat <- S * Mask

# Estimate the subspace spanned by the largest eigenvector using convex relaxation and TGD
# First calculate ground truth
result = geigen::geigen(Sigma,Sigma0)
evalues <- result$values
evectors <-result$vectors
evectors <- evectors[,p:1]
a <- evectors[,1:r]

## Running initialization using convex relaxation

#ag <- sgca_init(A=S, B=sigma0hat, rho=sqrt(log(p)/n),K=r ,nu=1,trace=FALSE, maxiter = 30)
ainit <- init_process(ag$Pi, r)
sqx <- sqrtm(Sigmax)$B
sqy <- sqrtm(Sigmay)$B

print('Initial prediction error of U is')
print( subdistance(sqx %*% initu, sqx %*% u)^2)
print('Initial prediction error of V is')
print( subdistance(sqy %*% initv, sqy %*% v)^2)

## Perform TGD

lambda <- 0.1
k <- 20

ainit = rbind(U_mine, V_mine)
scale <- a %*% pracma::sqrtm(diag(r)+t(a) %*% Sigma %*% a/lambda)$B;
final <- sgca_tgd(A=S, B=sigma0hat,r,ainit,k,lambda = 0.01, eta=0.00025,convergence=1e-6,maxiter=12000, plot = TRUE)

init <- gca_to_cca(ainit, S, pp)
final <- gca_to_cca(final, S, pp)
initu<- init$u
initv <- init$v
finalu <- final$u
finalv <- final$v
sqx <- pracma::sqrtm(Sigmax)$B
sqy <- pracma::sqrtm(Sigmay)$B

print('Initial prediction error of U is')
print( subdistance(sqx %*% initu, sqx %*% u)^2)
print('Initial prediction error of V is')
print( subdistance(sqy %*% initv, sqy %*% v)^2)

print('The final prediction error on U is')
print( subdistance(sqx %*% finalu , sqx %*% u)^2)

print('The final prediction error on V is')
print( subdistance(sqy %*% finalv, sqy %*% v)^2)

##### Suppose I do an adaptive lasso

library(glmnet)
cv.glmnet(X, Y%*% V_mine[,1])
#### Turn it into a regression setting
res <- glmnet(X, Y%*% V_mine[,1], lambda=0.03662, intercept = FALSE)
y_flat = as.vector(t( Y%*% V_mine))
sol = group_lasso(X, Y%*% V_mine, lambda=32.5)
print( subdistance(sqx %*% sol$Uhat , sqx %*% u)^2)

new_z = matrix(0, ncol(X), r)
for(i in 1:ncol(X)){
  for(j in 1:r){
    E = matrix(0,nrow(X), r)
    E[i,j]=1
    new_z[i,j] = t(X)%*% E
  }
}
grps <-  sapply(1:ncol(Y), function(x){rep(x,r)})

print( subdistance(sqx %*% as.vector(coef(res))[2:(p1+1)] , sqx %*% u[,1])^2)
print( subdistance(sqx %*% U_mine[,1], sqx %*% u[,1])^2)
print( subdistance(sqx %*% finalu[,1], sqx %*% u[,1])^2)

library(CVXR)
group_lasso <- function(X, Y,
                        lambda, max.iter=5000,
                         verbose = FALSE,
                        ZERO_THRESHOLD=1e-5){
  p <- ncol(X)
  k <- ncol(Y)
  ## Define_Parameters
  Uhat <- Variable(p, k)
  ## Build_Penalty_Terms
  penalty_term2 <- sum(cvxr_norm(Uhat, 2, axis = 1))
  
  ## Compute_Fitted_Value
  ## Build_Objective
  objective <- 1 / 2 * sum_squares( X %*% Uhat - Y) + lambda  * penalty_term2
  ## Define_and_Solve_Problem
  prob <- Problem(Minimize(objective))
  result <- psolve(prob, verbose = TRUE, num_iter =max.iter)
  ## Return_Values
  Uhat <- result$getValue(Uhat)
  ## Zero out stuff before returning
  Uhat[abs(Uhat) < ZERO_THRESHOLD] <- 0.0
  list(
    Uhat = Uhat,
    criterion = result$value)
}
sol = group_lasso(X, Y%*% V_mine, lambda=32.5)
norms = apply(sol$Uhat, 1, function(x){sum(x^2)})
norms[which(norms==0)] = 1e-4
sol2 = group_lasso(X%*% diag(norms), Y%*% V_mine, lambda=5)
print( subdistance(sqx %*% sol2$Uhat, sqx %*% u)^2)
