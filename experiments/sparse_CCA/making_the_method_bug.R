library(MASS)
library(geigen)
library(pracma)
library(VGAM)
library(mvtnorm)


generate_difficult_example <- function(n, p1, p2,   nnzeros = 5,
                             theta = diag( c(0.9,  0.8)),
                             r=2, 
                             signal_strength="normal"){
  # n  <- 200;
  # p1 <- 100;
  # p2 <- 100;
  #  r <- 2
  # a = 0.3
  p <- p1 + p2;
  pp <- c(p1,p2);
  print("We have:")
  print(c(p, nnzeros<=min(p1,p2), nnzeros))
  print('--------------------------------------');
  print('Generating data ...');
  #a <- 0.3;
  Sigma <- diag(p1+p2)
  s = (1):(nnzeros)
  s_cross = (3*nnzeros+1):(4*nnzeros)
  
  # generate covariance matrix for X and Y
  u1 = matrix(0, pp[1], r)
  u1[s,1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  u1 <- u1 %*%(sqrtm(t(u1[s,1:r]) %*% u1[s,1:r])$Binv)
  #T1 = eye(p1)#
  T1 = diag(rep(0.05, pp[1])) +  u1 %*% diag(rep(0.9, r)) %*% t(u1)
  
  v1 = matrix(0, pp[2], r)
  v1[s,1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  v1 <- v1 %*%(sqrtm(t(v1[s,1:r]) %*% v1[s,1:r])$Binv)
  #T1 = eye(p1)#
  T2 = diag(rep(0.05, pp[2])) +  v1 %*% diag(rep(0.8, r)) %*% t(v1)
  
  #T1 = as.matrix(toeplitz(a^(0:(pp[1]-1))));
  
  #T1[which(T1<1e-6)] = 0
  Sigma[1:p1, 1:p1] = T1;
  #T1 = Sigma[1:p1, 1:p1]
  Tss = T1[s,s];
  TsMs = T1[s,-s];
  TMsMs =  T1[-s,-s];
  u = matrix(0, pp[1], r)
  
  u[s_cross,1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  u <- u %*%(sqrtm(t(u[s_cross,1:r]) %*% Tss%*% u[s_cross,1:r])$Binv)
  
  #T2 = toeplitz(a^(0:(pp[2]-1)));
  Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)] = T2;
  Tss = Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)][s_cross, s_cross];
  v = matrix(0, pp[2], r)
  v[s_cross,1:r] <- as.matrix(runif( nnzeros * r,max = 3, min=1), nrow=nnzeros)  * as.matrix(sample(c(-1,1), nnzeros*r, replace=TRUE), nrow=nnzeros)
  v <- v %*%(sqrtm(t(v[s_cross,1:r]) %*% Tss %*% v[s_cross,1:r])$Binv)
  
  Sigma[(p1+1):(p1+p2), 1:p1] = T2 %*%  v  %*% theta %*% t(u) %*% T1;
  Sigma[1:p1, (p1+1):(p1+p2)] = t(Sigma[(p1+1):(p1+p2), 1:p1])
  
  Sigmax = Sigma[1:p1,1:p1];
  Sigmay = Sigma[(p1+1):p,(p1+1):p];
  
  #Generate Multivariate Normal Data According to Sigma
  Data = mvrnorm(n, rep(0, p), Sigma);
  
  X = Data[,1:p1];
  Y = Data[,(p1+1):(p1+p2)];
  
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
  #result = geigen::geigen(Sigma, Sigma0)
  #evalues <- result$values
  #evectors <-result$vectors
  #evectors <- evectors[,p:1]
  #a <- evectors[,1:r]
  #scale <- a %*% sqrtm(diag(r)+t(a) %*% Sigma %*% a/lambda)$B;
  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data, u=u, v=v, 
              Sigmax=Sigmax, Sigmay=Sigmay#,#a=a
  ))
}

p1=100
p2=100

example <- generate_example( n, p1, p2,   nnzeros = 10,
                                       theta = diag( c(1,  0.6)),
                                       r=2, 
                                       signal_strength="normal")



####
q = 50
S =  mvrnorm(q, mu =rep(0, p2), Sigma = diag(rep(1/q, p2)))
S1 =  mvrnorm(q, mu =rep(0, p1), Sigma = diag(rep(1/q, p1)))

transformed_covariance =  S1 %*% cov(example$X, example$Y) %*% t(S)
transf = rsvd(cov(example$X, example$Y) %*% t(S), k=2)
estimate_V = transf$v


###
### We have:
q = 30
for (it in 1:100){
  S =  mvrnorm(q, mu =rep(0, p2), Sigma = diag(rep(1/q, p2)))
  S1 =  mvrnorm(q, mu =rep(0, p1), Sigma = diag(rep(1/q, p1)))
  
  transformed_covariance =  S1 %*% cov(example$X, example$Y) %*% t(S)
  transf = svd(S1 %*% cov(example$X, example$Y) %*% t(S), nu=8, nv=8)
  estimate_V = transf$v
  
  library(CVXR)
  # Define variables
  U <- Variable(p1, 8)
  # Define the objective function
  Ytilde = example$Y %*% t(S) %*% estimate_V
  lambda = 1
  objective <- 1/nrow(example$X) * sum_squares(example$X %*% U - Ytilde) + lambda * sum(norm2(U, 1))
  # Define the problem
  problem <- Problem(Minimize(objective))
  
  # Solve the problem
  result <- solve(problem)
  
  # Extract the optimized U matrix
  optimal_U_temp <- result$getValue(U)
  #optimal_U_temp = optimal_U_temp %*% sqrtm(t(optimal_U_temp) %*% optimal_U_temp)$Binv
  optimal_U_temp[which(abs(optimal_U_temp)<1e-5)]=0
  print(subdistance(optimal_U_temp, example$u[1:p1,]))
  if(it==1){
    U_res = optimal_U_temp/100
  }else{
    U_res = U_res + optimal_U_temp/100
  }
  
    
}

U_res = U_res %*% sqrtm(t(U_res) %*% U_res)$Binv
subdistance(U_res, example$u[1:p1,])

subdistance(U_res, example$u[1:p1,])
subdistance(svd(cov(example$X, example$Y), nu=2, nv=2)$u, example$u[1:p1,])
transf = rsvd(cov(example$X, example$Y) %*% t(S), k=2)
new_mat = cov(example$X, example$Y) %*% t(S)

ag <- sgca_init(A= (new_mat) %*% t(new_mat), B=diag(rep(p2,1)), 
                rho=0.5*sqrt(log(p2)/q),K=r,
                nu=1,trace=FALSE)
ainit <- init_process(ag$Pi, r)
uest <- ainit
uest[which(abs(uest)<1e-5)]=0
uest = uest %*% sqrtm(t(uest) %*% uest)$Binv
subdistance(ainit[1:p1,], example$u[1:p1,])
# Running initialization using convex relaxation
p=p1+p2
ag <- sgca_init(A= example$S, B=example$sigma0hat, 
                rho=0.5*sqrt(log(p)/n),K=r,
                nu=1,trace=FALSE)
ag <- sgca_init(A= example$S, B=diag(rep(1, nrow(example$S))), 
                rho=0.5*sqrt(log(p1+p2)/n),K=r,
                nu=1,trace=FALSE)


ainit <- init_process(ag$Pi, r)
uest <- ainit
uest[which(abs(uest)<1e-5)]=0
uest = uest %*% sqrtm(t(uest) %*% uest)$Binv
subdistance(uest[1:p1,], example$u[1:p1,])



dest <- diag(ainit$d)
if (r == 1){
  ainit <- uest[,1] * sqrt(dest[1:r,1:r])
} else
  ainit <- uest[,1:r] %*% sqrtm(dest[1:r,1:r])$B

ainit <- init_process(ag$Pi, r)
ainit[which(abs(ainit)<1e-6)] = 0
fant = data.frame(ainit)
fant["x"] = 1:nrow(fant)
ggplot(pivot_longer(fant,
                    cols=-c("x")), aes(x=name, y=x, fill=value)) + 
  geom_tile()



a = rbind(example$u, example$v)
real = data.frame(a)
real["x"] = 1:nrow(real)
ggplot(pivot_longer(real,
                    cols=-c("x")), aes(x=name, y=x, fill=value)) + 
  geom_tile()


#### So the method sucks and failed
#### Idea: treat this as a combinatorial problem. I have a bunch of observations 
#### that I want to denoise





print('The initial error is')
print(subdistance(ainit, a))

## Running Thresholded Gradient Descent
final <- sgca_tgd(A=example$S, B=example$sigma0hat,r,ainit,k,
                  lambda = 1, 
                  eta=0.001,convergence=1e-6,maxiter=15000, plot = TRUE)
print('The final error is')
print(subdistance(final, a))

q = length(s) * 3
Omega1 <- matrix(rnorm(p1 * q, sd = 1/sqrt(q)), nrow=p1)
Omega2 <- matrix(rnorm(p2 * q, sd = 1/sqrt(q)), nrow=p2)

sketch = t(Omega1) %*% t(example$X)%*%example$Y %*% (Omega2)

example <- generate_example( n, p1, p2,   nnzeros = 5,
                                       theta = diag( c(0.8,  0.8)),
                                       r=2, 
                                       signal_strength="normal")
example <- generate_difficult_example( n, p1, p2,   nnzeros = 5,
                             theta = diag( c(0.8,  0.8)),
                             r=2, 
                             signal_strength="normal")
test = t(scale(example$X)) %*% scale(example$Y)/n
test =  t(scale(example$X)) %*% scale(example$Y)/n
test = cov(X, Y)
principal_a = test %*% t(test)
ag <- sgca_init(A=principal_a, B=diag(rep(1, nrow(principal_a))), 
                rho=0.5*sqrt(log(p1)/p2),K=r,
                nu=1,trace=FALSE)
svd_sketch = svd( test)
plot(svd_sketch$d)
subdistance(svd_sketch$u[, 1:2], example$u)
subdistance(svd_sketch$v[, 1:2], example$v)

plot(apply(svd_sketch$u[,1:20], 1, function(x){sum(x^2)}))
plot(apply(svd_sketch$v[,1:20], 1, function(x){sum(x^2)}))

test = sqrtm(example$Sigma[1:p1, 1:p1])$Binv %*% t(example$X)%*%example$Y %*% sqrtm(example$Sigma[(p1+1):(p1+p2), (p1+1):(p1+p2)])$Binv/n
res = ssvd(example$Sigma[1:p1, (p1+1):p_tot] , r=2, method = "method") #### This is not working out
res = ssvd(test , r=2, method = "theory") #### This is not working out
subdistance(res$u[, 1:2], example$u)
subdistance(res$v[, 1:2], example$v)

res2 = rsvd(test , k=10)
subdistance(res2$u[, 1:2], example$u)
subdistance(res2$u[, 1:2], example$u)
res2 = rcur(test , k=10)
test = data.frame(test)

ag <- sgca_init(A=example$S, B=example$sigma0hat, 
                rho=0.5*sqrt(log(p)/n),K=r,
                nu=1,trace=FALSE)
ainit = svd(ag$Pi)
uest <- ainit$u
dest <- diag(ainit$d)
if (r == 1){
  ainit <- uest[,1] * sqrt(dest[1:r,1:r])
} else
  ainit <- uest[,1:r] %*% sqrtm(dest[1:r,1:r])$B

ainit <- init_process(ag$Pi, r)
subdistance(ainit[1:p1, 1:2], example$u)


test["x"] = 1:nrow(test)
ggplot(pivot_longer(test,
                    cols=-c("x")), aes(x=name, y=x, fill=value)) + 
  geom_tile()

res = data.frame(res$u)
res["x"] = 1:nrow(res)
ggplot(pivot_longer(res,
                    cols=-c("x")), aes(x=name, y=x, fill=value)) + 
  geom_tile()


res2 = data.frame(res2$u)
res2["x"] = 1:nrow(res2)
ggplot(pivot_longer(res2,
                    cols=-c("x")), aes(x=name, y=x, fill=value)) + 
  geom_tile()
library(rsvd)






generate_grid_example<- function(n, p,    nnzeros = 5,
                             theta = diag( c(0.9,  0.8)),
                             a = 0, r=2, signal_strength="normal"){
  #n  <- 200;
  #p1 <- 100;
  #p2 <- 100;
  # r <- 2
  library(igraph)
  # Create a 2D grid graph
  
  g <- make_lattice(c(p, p))
  p1 = vcount(g)
  Sigma <- diag( 2 * vcount(g))
  ## Plot the grid graph
  #plot(g)
  # Select a random node
  random_node <- sample(V(g), r  * 2)
  # Find the 2-hop neighborhood of the selected node
  two_hop_neighborhood <- ego(g, order = 2, nodes = random_node, mindist = 0)
  
  
  p_tot <- 2* vcount(g)
  u = matrix(0, vcount(g), r)
  v = matrix(0, vcount(g), r)
  
  
  T1 = eye(p1) #as.matrix(toeplitz(a^(0:(pp[1]-1))));
  T1[which(T1<1e-6)] = 0
  Sigma[1:p1, 1:p1] = T1;
  T2 = eye(p1) #as.matrix(toeplitz(a^(0:(pp[1]-1))));
  T2[which(T1<1e-6)] = 0
  Sigma[(p1+1):(2*p1), (p1+1):(2*p1)] = T2;
  #T1 = Sigma[1:p1, 1:p1]
  
  
  for (i in 1:length(two_hop_neighborhood)){
    nnzeros = as.numeric(two_hop_neighborhood[[i]])
    print( i%% r)
    u[nnzeros, i%% r+1] <- runif( 1, max = 3, min=1)  * sample(c(-1,1), 1, replace=TRUE)  #runif( length(nnzeros), max = 3, min=1)  * sample(c(-1,1), length(nnzeros), replace=TRUE)
    v[nnzeros, i%% r+1] <- runif( 1, max = 3, min=1)  * sample(c(-1,1), 1, replace=TRUE)
  }
  s = unique(which(u!=0, arr.ind = TRUE))[,1]
  Tss = T1[s,s];
  T2ss = T2[s,s];
  u[s,] <- u[s,] %*%(sqrtm(t(u[s,]) %*% Tss%*% u[s,])$Binv)
  v[s,] <- v[s,] %*%(sqrtm(t(v[s,]) %*% T2ss%*% v[s,])$Binv)
  
  ###
  #V(g)$label = u[,2]
  #pal <- colorRampPalette(c("blue", "red"))
  #V(g)$color <- pal(100)[floor((V(g)$label - min(V(g)$label)) / diff(range(V(g)$label)) * 99) + 1]
  # Plot the graph with labels and node colors
  #plot(g, vertex.label = round(V(g)$label, 2), vertex.color = V(g)$color, main = "Graph with Colored Nodes")
  #plot(g,  vertex.color = V(g)$label)
  
  Sigma[(p1+1):(p1*2), 1:p1] = T2 %*%  v  %*% theta %*% t(u) %*% T1;
  Sigma[1:p1, (p1+1):(p1*2)] = t(Sigma[(p1+1):(2 * p1), 1:p1])
  
  Sigmax = Sigma[1:p1,1:p1];
  Sigmay = Sigma[(p1+1):p,(p1+1):p];
  
  #Generate Multivariate Normal Data According to Sigma
  Data = mvrnorm(n, rep(0, 2 * p1), Sigma);
  
  X = Data[,1:p1];
  Y = Data[,(p1+1):(p_tot)];
  
  print('Data generated.');
  print('--------------------------------------');
  
  Mask = matrix(0, p_tot, p_tot);
  idx1 = 1:p1;
  idx2 = (p1+1):(p_tot);
  Mask[idx1,idx1] <- matrix(1, p1,p1);
  Mask[idx2,idx2] <- matrix(1, p1,p1);
  Sigma0 = Sigma * Mask;
  
  
  S <- cov(Data)
  sigma0hat <- S * Mask
  
  
  all_edges = rbind(as_edgelist(g), cbind(as_edgelist(g)[,2], as_edgelist(g)[,1]))
  
  all_edges = data.frame(all_edges)
  all_edges["e"] = 1:nrow(all_edges)
  all_edges["fill"] = 1
  source = pivot_wider(all_edges[c("X1", "e", "fill")], id_cols = c("e"), 
                       values_from = "fill", 
                       names_from  = "X1",
                       values_fill=0)
  target = pivot_wider(all_edges[c("X2", "e", "fill")], id_cols = c("e"), 
                       values_from = "fill", 
                       names_from  = "X2",
                       values_fill=0)
  edge_incidence = source[,2:ncol(source)] - target[,2:ncol(target)]
  
  
  return(list(Sigma=Sigma, Sigma0=Sigma0,
              S = S, sigma0hat =  sigma0hat, Mask= Mask,
              X=X, Y = Y, Data=Data, u=u, v=v, 
              Sigmax=Sigmax, Sigmay=Sigmay,
              edge_incidence=edge_incidence#,#a=a
  ))

}

  
 


example <-generate_grid_example(n, p=11, theta = diag( c(0.9,  0.8)),  a = 0, r=2, signal_strength="normal")

#### This has all the edges

### Try
apply(example$S, 2, sum)

U =  matrix(rep(1, nrow(example$S)), nrow=nrow(example$S))/ sqrt(nrow(example$S))
AU =  t(U) %*% (example$S)
AU  = apply(example$S[1:p1,], 1, sum)/ sqrt(p1)
#svd_AU = svd(AU, 1)
#AU_k = svd_AU$d * (svd_AU$u  %*% t(svd_AU$v)) %*% t(U) #### principal subspace 1
#AU_k = AU_k/sqrt(sum(AU_k^2))
#solution = svd_AU$d * t(svd_AU$v)  %*% matrix(rep(1, nrow(example$S)), nrow=nrow(example$S))

####
test = as.matrix(example$edge_incidence) %*% example$S[1:p1, p1:p_tot]
res = ssvd(test, r=2, method = "theory") #### This is not working out
Utilde = res$u
Utilde2 = t(edge_incidence) %*% pinv((edge_incidence) %*% t(edge_incidence)) %*% res$u
Utilde2 = Utilde2 %*% sqrtm(t(Utilde2) %*% Utilde2)$Binv


edge_incidence = as.matrix(edge_incidence)

res = t(edge_incidence) %*% pinv((edge_incidence) %*% t(edge_incidence)) %*% res$u



ag <- sgca_init(A=test2, B=diag(nrow(test2)), 
                rho=0.5*sqrt(log(nrow(test2))/n),K=r,
                nu=1,trace=FALSE)
ainit = svd(ag$Pi)
uest <- ainit$u
dest <- diag(ainit$d)
if (r == 1){
  ainit <- uest[,1] * sqrt(dest[1:r,1:r])
} else
  ainit <- uest[,1:r] %*% sqrtm(dest[1:r,1:r])$B

ainit <- init_process(ag$Pi, r)


matrix(rep(1, nrow(example$S)), nrow=nrow(example$S)) %*% AU_k




