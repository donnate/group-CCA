library(ggplot2)
library(dplyr)
library(fields)
source("CCCA_reweighting.R")

set.seed(0)
n = 1000
p = 100
G = 5
s2 = 10
mu = matrix(rnorm(n * G, 0, s2), n, G)
X = c()
for(g in 1:G){
  X = cbind(X, replicate(p/G, rnorm(n = n, mean = mu[,g], s2)))
}
cluster = rep(1:G, times = rep(p/G, G))
R = matrix(rnorm(G * G), G, G)
Y = mu %*% R + matrix(rnorm(n * G, 0, s2), n, G)
image.plot(cor(X, Y))

coef_paths <- function(X, Y, lambdas, gammas, A, cluster){
  p = ncol(X)
  alphas <- c()
  alpha <- rnorm(p)
  for(gamma in gammas){
    for(lambda in lambdas){
      cc <- CCA_cluster(X, Y, lambda, gamma, A)
      alpha <- cc$alpha * sign(cor(alpha, cc$alpha))
      alphas <- rbind(alphas, data.frame(alpha, gamma, feature = 1:p, cluster = cluster, lambda))
    }
  }
  plt <- alphas %>% mutate(cluster = factor(cluster)) %>%
    ggplot(aes(x = log(lambda, 10), y = alpha, group = feature, color = cluster))+
    theme(legend.position = "none")+
    geom_line()+
    facet_wrap(~gamma)
  return(list(alphas = alphas, plt = plt))
}



#unknown clusters
lambdas = 10^(seq(-7,2))
gammas =  c(0, 0.25, 0.5, 0.75, 1)
uc = coef_paths(X, Y, lambdas, gammas, 1, cluster)
uc$plt 

#known clusters
A = matrix(0, p, p)
for(g in 1:G){
  A[(p/G*(g-1)+1):(p/G*g),(p/G*(g-1)+1):(p/G*g)] = 1
}
image.plot(A)
kc = coef_paths(X, Y, lambdas, gammas, A, cluster)
kc$plt 

#unknown clusters, weighted
#W = as.matrix(dist(t(X), diag = T, upper = T))
A = abs(cor(X))
diag(A) = 0
image.plot(A)

ucw = coef_paths(X, Y, lambdas, gammas, A, cluster)
ucw$plt 
