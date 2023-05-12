library(dplyr)
library(plotly)
library(fields)
source("functions.R")

dir = 'Data/'
coord = read.csv(paste0(dir, 'loci_3d_coord_reduced10.csv'), header = TRUE) 
X = t(read.csv(paste0(dir, 'activations_reduced10.csv'), header = TRUE))
Y = t(read.csv(paste0(dir, 'behavior_scores.csv'), header = TRUE, row.names = 1)) %>% 
  data.frame() %>%
  dplyr::select(starts_with("wn")) %>%
  as.matrix()
keep = names(which(rowSums(is.na(Y)) == 0))
keep = intersect(rownames(X), keep)
X = X[keep,]
Y = Y[keep,]
X = scale(X, center = T, scale = F)
Y = scale(Y, center = T, scale = F)

p = ncol(X)
r = 1
A = (as.matrix(dist(coord, upper = T)) <= r) * 1
image.plot(A)

############CCCA cross-validation################

set.seed(1)
nfold = 5
split = sample(1:nfold, nrow(X), replace = TRUE)

lambdas = 10^seq(-6,2)
mus =  10^seq(-6,2)
gammas =  c(0, 0.25, 0.5, 0.75, 1)

perf = c()
for(fold in 1:nfold){
  train = (split != fold)
  cat("FOLD:", fold, "\n")
  #RCCA
  rcca = CCCA_gridsearch(X[train,], Y[train,], X[!train,], Y[!train,], method = "RCCA", lambdas, gammas = 0, A = 0, mus)
  #CCCA
  ccca = CCCA_gridsearch(X[train,], Y[train,], X[!train,], Y[!train,], method = "CCCA", lambdas, gammas, A = 1, mus, eps = 1e-5, verbose = F)
  cccaf = CCCA_gridsearch(X[train,], Y[train,], X[!train,], Y[!train,], method = "CCCA", lambdas, gammas, A = A, mus, eps = 1e-5, verbose = F)
  perf = rbind(perf, data.frame(rbind(rcca, ccca, cccaf), fold))
  write.csv(result, "Results/brain-CCCA.csv", row.names = F)
}
