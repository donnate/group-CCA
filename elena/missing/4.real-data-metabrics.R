library(ggplot2)
library(dplyr)
library(fields)
library(mixOmics)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
library(readr)
library(matrixStats)
source("helper.R")
source("CCAimpute.R")


data = read_csv("Data/METABRIC_RNA_mutation.csv")
X = data %>% dplyr::select(-(1:31)) %>% 
  dplyr::select(-contains("_mut")) %>% 
  as.matrix()
Y = data %>% dplyr::select(age_at_diagnosis, nottingham_prognostic_index, neoplasm_histologic_grade, lymph_nodes_examined_positive, mutation_count, overall_survival_months, tumor_size, tumor_stage) %>%
  as.matrix()
colMeans(is.na(Y))
complete = which(rowSums(is.na(Y)) == 0)
X = X[complete,]
Y = Y[complete,]

n = nrow(X)
p = ncol(X)
q = ncol(Y)

props = seq(0.05, 0.3, 0.05)
lambda1s = 10^seq(-3,3,1) 
set.seed(0)
nfold = 5
folds = sample(rep(1:nfold, rep(length(complete)/nfold)))
seeds = 1

compare = function(Xna, Yna, Xtest, Ytest, lambda1s){
  mX = colMeans(Xna, na.rm = T)
  Xna = scale(Xna, center = mX, scale = F)
  mY = colMeans(Yna, na.rm = T)
  sY = colSds(as.matrix(Y), na.rm = T)
  Yna = scale(Yna, center = mY, scale = sY)
  Ximp = data.frame(Xna) %>% mutate_all(~replace_na(., 0))
  Yimp = data.frame(Yna) %>% mutate_all(~replace_na(., 0))
  Xtest = scale(Xtest, center = mX, scale = F)
  Ytest = scale(Ytest, center = mY, scale = sY)
  
  result = c()
  for(lambda1 in lambda1s){
    cat("\nlambda1:", lambda1)
    rrr = CCAimpute(Xna, Yna, Kx = diag(lambda1, p), Ky = NULL, eps = 1e-4, verbose = F)
    #Xtestrrr = scale(Xtest, center = colMeans(rrr$X), scale = F)
    #Ytestrrr = scale(Ytest, center = colMeans(rrr$Y), scale = F)
    result = rbind(result, data.frame(evaluate_nogold(Xtest, Ytest, rrr$U, rrr$V),
                                      lambda1, method = "RRR"))
    
    cca = CCA::rcc(Ximp, Yimp, lambda1 = lambda1, lambda2 = 0)
    result = rbind(result, data.frame(evaluate_nogold(Xtest, Ytest, cca$xcoef, cca$ycoef),
                                      lambda1, method = "RCCA"))
  }
  result
}


result = c()
for(fold in 1:nfold){
  cat("\n\nfold", fold, "prop:", 0, "\n")
  test = (folds == fold)
  Xtrain = X[-test,]
  Ytrain = Y[-test,]
  Xtest = X[test,]
  Ytest = Y[test,]
  result = rbind(result, data.frame(compare(Xtrain, Ytrain, Xtest, Ytest, lambda1s), prop = 0, fold, seed = 0))
  for(prop in props){
    for(seed in seeds){
      cat("\n\nfold", fold, "seed:", seed, "prop:", prop, "\n")
      
      maskX = sample(1:(nrow(Xtrain)*p), nrow(Xtrain)*p*prop)
      Xna = Xtrain
      Xna[maskX] = NA
      maskY = sample(1:(nrow(Ytrain)*q), nrow(Ytrain)*q*prop)
      Yna = Ytrain
      Yna[maskY] = NA
      
      result = rbind(result, data.frame(compare(Xna, Yna, Xtest, Ytest, lambda1s), prop, fold, seed))
      write.csv(result, "Fits/real-data-metabrics-RRR-cv-randomized.csv", row.names = F)
    }
  }
}


