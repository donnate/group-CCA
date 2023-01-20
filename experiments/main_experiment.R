setwd("~/Documents/group-CCA/")

library(tidyverse)
library(CCA)
library(VGAM)
library(matlib)
library(igraph)
library(recipes)
library(parallel)
library(tictoc)
library(PMA)
library(mvtnorm)
library(glmnet)
library(CCA)
library(pls)
library(igraph)

numCores <- detectCores()


#### parameters
n=250
p= 337
q=10 
sigma=0.1 
k=3
sigma_noise=0.1
power=1.6
lambdaA1seq= c(0.0001, 0.001, 0.01, 0.1, 1)
lambdaA2seq= c(0.0001, 0.001, 0.01, 0.1, 1)
conv=10^{-3} 
max.iter=200
egonet_size= 2
n.cv= 5
power=1
rank=k
lambdaA1=0.5
lambdaA2=1
plot.result = TRUE

vfn <- function(x){
  ifelse(x=="x", 1, -1)
}


G <- sample_pa(p, power = power)
E = data.frame(as_edgelist(G))
colnames(E)  = c("x", "y")
E["e"] = 1:nrow(E)
E = pivot_longer(E, cols=-c("e"))
E["fill_value"] = sapply(E$name, vfn)
D = pivot_wider(E, id_cols =c("e"), names_from = "value", values_from  =  fill_value)
D[is.na(D)] <- 0
D = as.matrix(D[, 2:(p+1)])

indices <- sample(1:p, 1)
#### Make sure the selected clusters are independent
not_all_indices = TRUE
while(not_all_indices){
  print("here")
  for (j in 2:k){
    found_vertex = FALSE
    iter = 1
    while(found_vertex==FALSE){
      iter  = iter + 1
      index <- sample(1:p, 1)
      d = sapply(indices, FUN=function(x){distances( G, v =x, to=index)})
      print(d)
      if(min(d)> 2 * egonet_size + 1){
        indices <- c(indices, index)
        found_vertex = TRUE
      }
      if(iter > 100){
        break;
      }
    }
  }
  not_all_indices = (length(indices) < k)
}
print(indices)


X = matrix(rnorm(p * n, mean=0, sd=sigma_noise), n, p)
Y = matrix(rnorm(q*n,  mean=0, sd=sigma_noise), n, q)
source = matrix(0, k, p)
colors = matrix(0, k, p)
colors_l = rep(0, p)
trueA = matrix(0, p, k)
trueB = matrix(0, q, k)
true_corr = rep(0, k)
for (i in 1:k){
  print(i)
  idx = indices[i]
  subg <- ego(G, order=egonet_size, nodes = idx, 
              mode = "all", mindist = 0)[[1]]
  nodes_in_network <- as.numeric(subg)
  colors_l[nodes_in_network] = i
  colors[i, nodes_in_network]  = i
  source[i, nodes_in_network] = 1
  #mean_value <- sapply(1:n, function(i){rnorm(1, mean=effect_size, sd=sigma)})
  #X[, nodes_in_network] <- mean_value + X[, nodes_in_network] 
  Y[, i] =Y[, i] + X[, nodes_in_network] %*% matrix(rep(1, length(nodes_in_network)), nrow=length(nodes_in_network), ncol=1)
  true_corr[i] = cor(Y[, i], X[, nodes_in_network] %*% matrix(rep(1, length(nodes_in_network)), nrow=length(nodes_in_network), ncol=1))
  trueA[nodes_in_network, i] =  1/sqrt(sum(X[, source[i,]]^2))
  trueB[i,i] = 1/ sqrt(sum(Y[,i]^2))
}


if(plot.result==TRUE){
  source_df = data.frame(source)
  source_df["component"] = 1:k
  ggplot(pivot_longer(source_df, cols=-c("component"))) +
    geom_tile(aes(x=name,y=component, fill=as.factor(value)))
  #### we can also simply plot the graph
  g <- simplify(G)
  V(g)$color= "black"
  V(g)$color[which(source[1,]>0)] =  "lightblue"
  V(g)$color[which(source[2,]>0)] =   "orange"
  V(g)$color[which(source[3,]>0)] =  "red"
  plot(g)
  
}
df.x <- data.frame(X) %>% mutate_all(~(scale(.) %>% as.vector))
df.y <- data.frame(Y) %>% mutate_all(~(scale(.) %>% as.vector))
X <- as.matrix(df.x)
Y <- as.matrix(df.y)


# Store results Estimation Accuracy
Nsim = 1
results.x <-list()
results.y <-list()
recovered_corr <- list()

#### divide data into train and test
train  = sample(1:n, size=as.integer(0.8 *n), replace = FALSE)
test = setdiff(1:n, train)
X_train = X[train, ]
X_test = X[test, ]
Y_train = Y[train, ]
Y_test = Y[test, ]


mysavedir = "~/Downloads/"
namefile ="testCCA_"

save.image(file=paste0(mysavedir, namefile, 'environment.RData'))

system(paste0("Rscript ~/Documents/group-CCA/experiments/run_experiment.R  ",namefile, 'environment', "   ", 'genCCA', '   ', mysavedir, '  ', rank,
              ' ', 0.01, ' ', 0.01, ' ', 0.0,  ' ', 200 ))
name_experiment  = args[1]
type_experiment = args[2]
mysavedir = args[3]
rank = args[4]
lambda1 = args[5]
lambda2 = args[6]
lambda3 = args[7]
max.iter = args[8]






