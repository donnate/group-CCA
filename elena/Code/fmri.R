library(dplyr)
library(plotly)
source("CCCA_reweighting.R")

dir = '../Data/HCP Disordered Emotional States/'
X = read.csv(paste0(dir, 'activation_avg.csv'), header = FALSE, row.names = 1)
Y = read.csv(paste0(dir, 'behavior.csv'), header = TRUE, row.names = 1) %>% 
  dplyr::select(-demo_age, -bio_sex)
coord = read.csv(paste0(dir, 'parcellation_coordinates.csv'), header = TRUE) %>% dplyr::select(x, y, z)
Y = na.omit(Y)
subs = intersect(rownames(X), rownames(Y))
X = as.matrix(X[subs,])
Y = as.matrix(Y[subs,])
p = ncol(X)


perf_paths <- function(Xtrain, Ytrain, Xtest, Ytest, type = "CCCA", lambdas, gammas = NULL, A = 1, cluster = 1){
  p = ncol(X)
  rhos <- c()
  alpha <- rnorm(p)
  if(type == "RCCA"){
    for(lambda in lambdas){
      cat("\nlambda:", lambda, "\n")
      cc <- CCA(Xtrain, Ytrain, K = diag(lambda, p))
      rhotrain <- cor(Xtrain %*% cc$alpha, Ytrain %*% cc$beta)
      rhotest <- cor(Xtest %*% cc$alpha, Ytest %*% cc$beta)
      rhos <- rbind(rhos, data.frame(rho = c(rhotrain, rhotest), type = c("train", "test"), lambda))
    }
  }
  if(type == "CCCA"){
    for(gamma in gammas){
      for(lambda in lambdas){
        cat("\nlambda:", lambda, "gamma:", gamma, "\n")
        cc <- CCA_cluster(Xtrain, Ytrain, lambda, gamma, A)
        rhotrain <- cor(Xtrain %*% cc$alpha, Ytrain %*% cc$beta)
        rhotest <- cor(Xtest %*% cc$alpha, Ytest %*% cc$beta)
        rhos <- rbind(rhos, data.frame(rho = c(rhotrain, rhotest), type = c("train", "test"), lambda, gamma))
      }
    }
  }
  return(rhos)
}


coef_paths <- function(X, Y, type = "CCCA", lambdas, gammas = NULL, A = 1, cluster = 1){
  p = ncol(X)
  alphas <- c()
  alpha <- rnorm(p)
  if(type == "RCCA"){
    for(lambda in lambdas){
        cat("\nlambda:", lambda, "\n")
        cc <- CCA(X, Y, K = diag(lambda, p))
        alpha <- cc$alpha * sign(cor(alpha, cc$alpha))
        alphas <- rbind(alphas, data.frame(alpha, lambda, feature = 1:p, cluster = cluster))
    }
  }
  if(type == "CCCA"){
    for(gamma in gammas){
      for(lambda in lambdas){
        cat("\nlambda:", lambda, "gamma:", gamma, "\n")
        cc <- CCA_cluster(X, Y, lambda, gamma, A)
        alpha <- cc$alpha * sign(cor(alpha, cc$alpha))
        alphas <- rbind(alphas, data.frame(alpha, lambda, gamma, feature = 1:p, cluster = cluster))
      }
    }
  }
  return(alphas)
}


############performance################
set.seed(1)
lambdas = 10^seq(-6,2)
gammas =  c(0, 0.25, 0.5, 0.75, 1)
nfold = 10
split = sample(1:nfold, nrow(X), replace = TRUE)


perf = c()
for(fold in 1:nfold){
  train = (split != fold)
  #RCCA
  rcca_perf = perf_paths(X[train,], Y[train,], X[!train,], Y[!train,], type = "RCCA", lambdas)
  #CCCA
  ccca_perf = perf_paths(X[train,], Y[train,], X[!train,], Y[!train,], type = "CCCA", lambdas, gammas)
  perf = rbind(perf, data.frame(rcca_perf, gamma = NA, fold), data.frame(ccca_perf, fold)) 
}

saveRDS(perf, "performance.rds")

perf %>% group_by(type, lambda, gamma) %>%
  summarise(rho_cv = mean(rho), se = sd(rho)/sqrt(nfold)) %>%
  mutate(method = ifelse(is.na(gamma), "RCCA", paste("CCCA, gamma =", gamma))) %>%
  ggplot(aes(log(lambda, 10), rho_cv, color = method, fill = method))+
  geom_ribbon(aes(ymin = rho_cv-se, ymax = rho_cv+se), alpha = 0.1, color = NA)+
  geom_line(size = 1)+
  facet_wrap(~type, ncol = 1, scale = "free")

############plots################
#RCCA
rcca_coef = coef_paths(X, Y, type = "RCCA", lambdas)
data.frame(coord, rcca_coef) %>%
  plot_ly(x = ~x, y = ~y, z = ~z, color = ~alpha, frame = ~lambda, type = "scatter3d", mode = "markers")

rcca_coef %>% mutate(cluster = factor(cluster)) %>%
  ggplot(aes(x = log(lambda, 10), y = alpha, group = feature, color = cluster))+
  theme(legend.position = "none")+
  geom_line()

#CCCA
ccca = coef_paths(X, Y, type = "CCCA", lambdas, gammas)
data.frame(coord, ccca) %>%
  filter(gamma == 1) %>%
  plot_ly(x = ~x, y = ~y, z = ~z, color = ~alpha, frame = ~lambda, type = "scatter3d", mode = "markers")

plt = ccca %>% mutate(cluster = factor(cluster)) %>%
  ggplot(aes(x = log(lambda, 10), y = alpha, group = feature, color = cluster, frame = gamma))+
  theme(legend.position = "none")+
  geom_line()
ggplotly(plt)

alpha_opt = ccca %>% filter(gamma == 1, lambda == 10^(-5)) %>% pull(alpha)
d = dist(alpha_opt, method = "euclidean")
hcl = hclust(d)
plot(hcl, cex = 0.6, hang = -1)

clusters = c()
for(k in 1:8) clusters = rbind(clusters, data.frame(cluster = cutree(hcl, k = k), k))

data.frame(coord, clusters) %>%
  filter(k == 3) %>%
  plot_ly(x = ~x, y = ~y, z = ~z, split = ~cluster, type = "scatter3d", mode = "markers")


