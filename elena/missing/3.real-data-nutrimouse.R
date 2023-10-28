library(ggplot2)
library(dplyr)
library(fields)
library(mixOmics)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
source("helper.R")
source("CCAimpute.R")

#Load X and Y from nutrimouse data
data(nutrimouse)
X = as.matrix(nutrimouse$gene %>% mutate_all(~rank(.)))
Y = as.matrix(nutrimouse$lipid %>% mutate_all(~rank(.)))
#X = scale(X, center = T, scale = F)
#Y = scale(Y, center = T, scale = F)
n = nrow(X)
p = ncol(X)
q = ncol(Y)

props = seq(0, 0.3, 0.05)
seeds = 1:30
lambda1s = c(0.001, 0.01, 0.1, 1) 
lambda2s = c(0.001, 0.01, 0.1, 1)
set.seed(0)
nfold = 10
#folds = sample(rep(1:nfold, rep(n/nfold)))

result = c()
for(prop in props){
  #for(fold in 1:nfold){ # can do regular CV or bootstrap-type
    #test = which(folds == fold)
    for(seed in seeds){
      set.seed(seed)
      folds = sample(rep(1:nfold, rep(n/nfold)))
      fold = 1
      test = which(folds == fold)
      cat("\n\nfold", fold, "seed:", seed, "prop:", prop, "\n")
      Xna = X[-test,]
      Yna = Y[-test,]
      Xtest = X[test,]
      Ytest = Y[test,]
      
      maskX = sample(1:(nrow(Xna)*p), nrow(Xna)*p*prop)
      Xna[maskX] = NA
      mX = colMeans(Xna, na.rm = T)
      Xna = scale(Xna, center = mX, scale = F)
      
      maskY = sample(1:(nrow(Yna)*q), nrow(Yna)*q*prop)
      Yna[maskY] = NA
      mY = colMeans(Yna, na.rm = T)
      Yna = scale(Yna, center = mY, scale = F)
      
      Ximp = data.frame(Xna) %>% mutate_all(~replace_na(., 0))
      Yimp = data.frame(Yna) %>% mutate_all(~replace_na(., 0))
      
      Xtest = scale(Xtest, center = mX, scale = F)
      Ytest = scale(Ytest, center = mY, scale = F)
        
      for(lambda1 in lambda1s){
        for(lambda2 in lambda2s){
          cat("\nlambda1:", lambda1, "lambda2:", lambda2)
          rrr = CCAimpute(Xna, Yna, Kx = diag(lambda1, p), Ky = diag(lambda2, q), eps = 1e-4, verbose = F)
          Xtestrrr = scale(Xtest, center = colMeans(rrr$X), scale = F)
          Ytestrrr = scale(Ytest, center = colMeans(rrr$Y), scale = F)
          
          result = rbind(result, data.frame(evaluate_nogold(Xtestrrr, Ytestrrr, rrr$U, rrr$V), 
                                            lambda1, lambda2, method = "RRR",  prop, seed, fold))
          cca = CCA::rcc(Ximp, Yimp, lambda1 = lambda1, lambda2 = lambda2)
          result = rbind(result, data.frame(evaluate_nogold(Xtest, Ytest, cca$xcoef, cca$ycoef), 
                                            lambda1, lambda2, method = "RCCA", prop, seed, fold))
          write.csv(result, "Fits/real-data-nutrimouse-ranks-RRR.csv", row.names = F)
        }
      }
    }
  }
#}
  
