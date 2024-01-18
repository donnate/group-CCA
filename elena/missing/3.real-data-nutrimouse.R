library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
setwd("~/Documents/group-CCA/")

source('experiments/sparse_CCA/experiment_functions.R')
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source("elena/missing/helper.R")
source("elena/missing/evaluation.R")
source("elena/missing/original_CCA_impute.R")
source("elena/gradient_descent.r")
source("elena/iterative_cca.R")
source("elena/reduced_rank_regression.R")
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

r = 5
#### Apply CCA
X = scale(X)
Y =scale(Y)
lambda_Kx = NULL
param_lambda=c(10^seq(-3, 3, length.out = 100))
kfolds=5
penalty="l21"
solver="rrr"
LW_Sy = TRUE

res_alt = CCA_rrr.CV(X, Y, 
                     r=r, Kx = NULL, lambda_Kx = 0,
                     param_lambda=c(10^seq(-3, 1, length.out = 10)),
                     kfolds=14, penalty="l21", solver="CVX", LW_Sy = TRUE)

res_alt$rmse

res_alt = CCA_rrr.CV(X, Y, 
                     r=r, Kx = NULL, lambda_Kx = 0,
                     param_lambda=c(10^seq(-3, 0.5, length.out = 10)),
                     kfolds=20, penalty="l21", solver="rrr", LW_Sy = TRUE)

res_alt = CCA_rrr.CV(X, Y, 
                     r=r, Kx = NULL, lambda_Kx = 0,
                     param_lambda=c(10^seq(-3, 0.5, length.out = 10)),
                     parallelize = FALSE,
                     kfolds=10, penalty="l21", solver="CVX", LW_Sy = TRUE)

library(mixOmics)
shrink.rcc.nutrimouse <- rcc(X,Y, ncomp=5, method = 'shrinkage') 
shrink.rcc.nutrimouse2_comp <- rcc(X,Y, method = 'shrinkage') 
# examine the optimal lambda values after shrinkage 
shrink.rcc.nutrimouse$lambda 


test = CCA_rrr(X, Y, Sx=NULL, Sy=NULL,
        lambda =3, Kx=NULL, r, highdim=TRUE, 
        penalty = "l21", lambda_Kx=0, solver="CVX",
        LW_Sy = TRUE)


Uhat = matrix(0, p, r)
index_zeros = which(apply(abs(test$U), 1, sum) >0)
#t = (varimax(test$U[index_zeros, ], normalize=FALSE))$loadings
for (i in 1:r){
  print(i)
  Uhat[,i] = test$U[,i]
}
###
XUhat = matrix(X, nrow=40) %*% Uhat
XU = matrix(X, nrow=40) %*% Uhat
#test$U <- data.frame(test$U)
XUhat = as.data.frame(XUhat)

Uhat = data.frame(Uhat)
rownames(Uhat) = colnames(X)
shrink.rcc.nutrimouse$loadings$X = Uhat
shrink.rcc.nutrimouse$variates$X = XU


colnames(XUhat) = sapply(1:r, function(rr){paste0("X", rr)})
#colnames(Uhat) = sapply(1:r, function(rr){paste0("X", rr)})
#Uhat["variable"]= colnames(X)
XUhat["diet"]= nutrimouse$diet
XUhat["genotype"]= nutrimouse$genotype

ggplot(XUhat, aes(x = X1, y = X2,
                  shape = genotype, colour = diet))+
  geom_point() +
  geom_jitter()

#ggplot(Uhat, aes(x = X1, y = X2))+
#  geom_point() +
#  geom_jitter()

Vhat = matrix(0, q, r)
index_zeros = which(apply(abs(test$V), 1, sum) >0)
t = (varimax(test$V[index_zeros, ], normalize=FALSE))$loadings
for (i in 1:r){
  print(i)
  Vhat[,i] = test$V[,i]
}
YVhat = matrix(Y, nrow=40) %*% Vhat
YV = matrix(Y, nrow=40) %*% Vhat
Vhat = as.data.frame(Vhat)
rownames(Vhat) = colnames(Y)
shrink.rcc.nutrimouse$loadings$Y = Vhat
shrink.rcc.nutrimouse$variates$Y = YV

YVhat = as.data.frame(YVhat)
colnames(Vhat) = sapply(1:r, function(rr){paste0("X", rr)})
colnames(YVhat) = sapply(1:r, function(rr){paste0("X", rr)})
#Vhat["variable"]= colnames(Y)
YVhat["diet"]= nutrimouse$diet
YVhat["genotype"]= nutrimouse$genotype

ggplot(YVhat, aes(x = X1, y = X2,
                  shape = genotype, colour = diet))+
  geom_point() +
  geom_jitter()

#ggplot(Vhat, aes(x = X1, y = X2))+
#  geom_point() +
#  geom_jitter()



test1<- additional_checks(X,
                          Y, S=NULL, 
                          rank=r, kfolds=10, 
                          method.type = "FIT_SAR_CV",
                          lambdax= 10^seq(-3,0.5, length.out = 10),
                          lambday = c(0, 0))

test1$u = test1$u %*% sqrtm(t(test1$u) %*% cov(X) %*%test1$u )$Binv
test1$v = test1$v %*% sqrtm(t(test1$v) %*% cov(Y) %*%test1$v )$Binv
Uhat_comp = matrix(0, p, r)
index_zeros = which(apply(abs(test1$u), 1, sum) >0)
t = (varimax(test1$u[index_zeros, ], normalize=FALSE))$loadings
for (i in 1:r){
  print(i)
  Uhat_comp[index_zeros,i] = t[,i]
}
XU_comp = matrix(X, nrow=40) %*% Uhat_comp
Vhat_comp = matrix(0, q, r)
index_zeros = which(apply(abs(test1$v), 1, sum) >0)
t = (varimax(test1$v[index_zeros, ], normalize=FALSE))$loadings
for (i in 1:r){
  print(i)
  Vhat_comp[index_zeros,i] = t[,i]
}
YV_comp = matrix(Y, nrow=40) %*% Vhat_comp
Uhat_comp = data.frame(Uhat_comp)
rownames(Uhat_comp) = colnames(X)
shrink.rcc.nutrimouse2_comp$loadings$X = Uhat_comp
shrink.rcc.nutrimouse2_comp$variates$X = XU_comp
shrink.rcc.nutrimouse2_comp$variates$Y = YV_comp
Vhat_comp = data.frame(Vhat_comp)
rownames(Vhat_comp) = colnames(Y)
shrink.rcc.nutrimouse2_comp$loadings$Y = Vhat_comp


ggplot()
plotIndiv(shrink.rcc.nutrimouse, comp = c(1,2), 
          ind.names = nutrimouse$genotype,
          group = nutrimouse$diet, rep.space = "XY-variate", 
          legend = TRUE, title = '(a) Our method  in XY-space',
          point.lwd = 1, cex=5)
plotIndiv(shrink.rcc.nutrimouse2_comp, comp = 1:2, 
          ind.names = nutrimouse$genotype,
          group = nutrimouse$diet, rep.space = "XY-variate", 
          legend = TRUE, title = '(b) Nutrimouse, SAR  XY-space')
ggplot()
network(shrink.rcc.nutrimouse, comp = 3, interactive = FALSE,
        lwd.edge = 2,
        cutoff = 0.45)
network(shrink.rcc.nutrimouse2_comp, comp = 3, interactive = FALSE,
        lwd.edge = 2,
        cutoff = 0.5)

cim(shrink.rcc.nutrimouse, comp = 1:5, xlab = "genes", ylab = "lipids")
cim(shrink.rcc.nutrimouse2_comp, comp = 1:5, xlab = "genes", ylab = "lipids")

comps = c(2, 5)
plotIndiv(shrink.rcc.nutrimouse, comp = comps, 
          ind.names = nutrimouse$genotype,
          ellipse = TRUE,  # plot using the ellipses
          group = nutrimouse$diet, rep.space = "XY-variate", 
          legend = TRUE, 
          title = '(a) Our method in XY-space',
          point.lwd = 1, cex=5)

# plot the arrow plot of samples for CV rCCA data
plotArrow(shrink.rcc.nutrimouse, group = nutrimouse$diet, 
          col.per.group = color.mixo(1:5),
          title = '(a) Nutrimouse, CV method',
          pch.size = 3, cex=5, ind.names.size = 3,
          legend.title = "Diet")

# plot the arrow plot of samples for shrinkage rCCA data
plotArrow(shrink.rcc.nutrimouse2_comp, group = nutrimouse$diet, 
          col.per.group = color.mixo(1:5),
          title = '(b) Nutrimouse, shrinkage method',
          pch.size = 3, cex=5, ind.names.size = 3,
          legend.title = "Diet",  pch = 1)


plotVar(shrink.rcc.nutrimouse, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.5,
        title = '(b) Nutrimouse, rCCA shrinkage comp 1 - 2')

ggplot(data.frame(pivot_longer(rbind(Uhat, Vhat), cols = sapply(1:r, function(rr){paste0("X", rr)})
                               )))+
  geom_tile(aes(x=variable, y=name, fill =value)) + 
  scale_fill_gradient2()

ggplot(data.frame(pivot_longer(Vhat, cols = c("X1", "X2", "X3"))))+
  geom_tile(aes(x=lipid, y=name, fill =value)) + 
  scale_fill_gradient2()


# plot the projection of samples for shrinkage rCCA data
# run the rCCA method using shrinkage
shrink.rcc.nutrimouse <- rcc(X,Y, method = 'shrinkage') 
# examine the optimal lambda values after shrinkage 
shrink.rcc.nutrimouse$lambda 

plotIndiv(shrink.rcc.nutrimouse, comp = 1:2, 
          ind.names = nutrimouse$genotype,
          group = nutrimouse$diet, rep.space = "XY-variate", 
          legend = TRUE, title = '(b) Nutrimouse, rCCA shrinkage XY-space')

method<-SparseCCA(X=as.matrix(X, nrow=40),
                  Y=as.matrix(Y, nrow=40),rank=2,
                  lambdaAseq=10^seq(-3,0.5, length.out = 10),
                  lambdaBseq= c(0, 0),
                  max.iter=100,conv=10^-2, selection.criterion=2, n.cv=5)



plotIndiv(shrink.rcc.nutrimouse, comp = 2:3, 
          ind.names = nutrimouse$genotype,
          group = nutrimouse$diet, rep.space = "XY-variate", 
          legend = TRUE, title = '(b) Nutrimouse, rCCA shrinkage XY-space')
network(shrink.rcc.nutrimouse, comp = 1:2, interactive = FALSE,
        lwd.edge = 2,
        cutoff = 0.5)
cim(shrink.rcc.nutrimouse, comp = 1:2, xlab = "genes", ylab = "lipids")

plotVar(shrink.rcc.nutrimouse, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.5,
        title = '(b) Nutrimouse, rCCA shrinkage comp 1 - 2')



plotArrow(shrink.rcc.nutrimouse, group = nutrimouse$diet, 
          col.per.group = color.mixo(1:5),
          title = '(b) Nutrimouse, shrinkage method')



Uhat_comp = as.data.frame(Uhat_comp)
colnames(Uhat_comp) = c("X1", "X2", "X3")
Uhat_comp["variable"]= colnames(X)


Vhat_comp = as.data.frame(Vhat_comp)
colnames(Vhat_comp) = c("X1", "X2", "X3")
Vhat_comp["variable"]= colnames(Y)


ggplot(data.frame(pivot_longer(rbind(Uhat_comp, Vhat_comp), cols = c("X1", "X2", "X3"))))+
  geom_tile(aes(x=variable, y=name, fill =value)) + 
  scale_fill_gradient2()


Uhat_SAR["method"] = "SAR"
Uhat["method"] = "RRR"
combined = rbind(Uhat_SAR, Uhat)
ggplot(data.frame(pivot_longer(combined, cols = c("X1", "X2", "X3"))))+
  geom_tile(aes(x=gene, y=name, fill =value)) + 
  scale_fill_gradient2() + 
  facet_grid(method~.)

ggplot(data.frame(pivot_longer(Vhat, cols = c("X1", "X2", "X3"))))+
  geom_tile(aes(x=lipid, y=name, fill =value)) + 
  scale_fill_gradient2()


Vhat_SAR["method"] = "SAR"
Vhat["method"] = "RRR"
combinedV = rbind(Vhat_SAR, Vhat)
ggplot(data.frame(pivot_longer(combinedV, cols = c("X1", "X2", "X3"))))+
  geom_tile(aes(x=lipid, y=name, fill =value)) + 
  scale_fill_gradient2() + 
  facet_grid(method~.)
t(res_alt$vfinal) %*% cov(as.matrix(Y)) %*% (res_alt$vfinal)

test1_v_norm <-test1$v   %*% sqrtm(t(test1$v) %*% cov(as.matrix(Y)) %*% (test1$v))$Binv
t(test1_v_norm) %*% cov(as.matrix(Y)) %*% (test1_v_norm)

ggplot(data.frame(pivot_longer(Vhat, cols = c("X1", "X2", "X3"))))+
  geom_tile(aes(x=lipid, y=name, fill =value)) + 
  scale_fill_gradient2()
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
  
