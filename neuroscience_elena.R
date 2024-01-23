library(tidyverse)
library(RNifti)
library(oro.nifti)
library(neurobase)

library(tidyverse)
library(igraph)
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
source("elena/graph_reduced_rank_regression.R")
#store all values above diagonal of connectome matrices in matrix c
dir = '../fMRI-data/data/'

activations <- read_csv2("~/Downloads/activation (1).csv")
behaviour <- read_csv("~/Downloads/behavior.csv")
group_assignment <- readxl::read_xlsx("~/Downloads/activation_groups.xlsx", col_names = FALSE)
colnames(group_assignment) <- "group.id"
index_groups = which(group_assignment$group.id!=0)
activations  = activations[,c(1, 1+index_groups)]
groups <- sapply(1:length(unique(group_assignment$group.id[index_groups])),
                 function(g){which(group_assignment$group.id[index_groups] == g)})

name_subject = intersect(activations$Row, behaviour$participant_id)


new_data = activations %>% inner_join(behaviour, join_by(Row == participant_id ))
X = new_data[, 2:ncol(activations)]
Y = new_data[, c("demo_age",  "bio_sex",  "bas_drive" ,
                 "bas_funseeking",  "bas_rewardrespons",
                 "bis_total", "masq_distress", "masq_anhedonia",   
                 "masq_anxarousal", "panas_positive","panas_negative")]

gender = Y$bio_sex
Y = Y[ , c("bas_drive" ,
        "bas_funseeking",  "bas_rewardrespons",
        "bis_total", "masq_distress", "masq_anhedonia",   
        "masq_anxarousal", "panas_positive","panas_negative")]
females = which(gender == "F")
males = which(gender == "M")
means_f = colMeans(X[females, ])
means_m = colMeans(X[males, ])
X[females, ] = X[females, ] + (means_m - means_f)

means_qf = colMeans(Y[females,  ], na.rm=TRUE)
means_qm = colMeans(Y[males, ], na.rm=TRUE)
Y[females, ] = Y[females, ] + (means_m - means_f)

# Calculate column means, excluding NAs
column_means <- colMeans(Y, na.rm = TRUE)

# Replace NA values with column means
for (i in 1:ncol(Y)){
  Y[is.na(Y[,i]), i] <- column_means[i]
}


STOP
GM_mask = readnii("/Users/cdonnat/Documents/group-CCA/experiments/real-data/fMRI-data/data/gm_mask020_bin.nii.gz")
GM_mask = readnii("/Users/cdonnat/Downloads/BN_Atlas_246_2mm.nii.gz")
GM.d = GM_mask@.Data
GM.d <- reshape2::melt(GM.d)
colnames(GM.d)<- c("x", "y", "z", "GM")
print(dim(GM.d))


res_alt <- CCA_group_rrr(X, Y, 
                         groups = groups[index_groups],
                         Sx = NULL, Sy=NULL,
                         lambda =lambda_chosen, Kx=NULL, 
                         r=r, 
                         solver="ADMM",
                         LW_Sy = LW_Sy, do.scale=TRUE,
                         rho=1, niter=2 * 1e4,
                         thresh=1e-6)


##### Grrrrr Need to adjust for gender

##### Let us approximate Sx


n = nrow(X)
p = ncol(X)
q = ncol(Y)
do.scale = T
if(q >n){
  X_temp <- X
  X <- Y
  Y <- X_temp
}
if (do.scale){
  X <- scale(X)
  Y <- scale(Y)
}


newX = t(X)
newX = data.frame(newX)
newX["group.id"] = group_assignment$group.id[index_groups]
#newX["ng"] = 1
newX_summ = newX %>% group_by(group.id) %>% 
  summarise_all(mean)
Sx = as.matrix(newX_summ[, -c(1)]) %*% t(as.matrix(newX_summ[, -c(1)]))/n


#test = pivot_longer(newX_summ, cols = -c("group.id", "ng"))
#test = test %>%
#  mutate(value = value/ sqrt(ng))
#test2 = pivot_wider(test, id_cols = c("group.id", "ng"),
#                    values_from = "value",
#                    names_from = "name")
#Sx = as.matrix(newX_summ[, -c(1, 2)]) %*% t(as.matrix(newX_summ[, -c(1, 2)]))/n

Sy = NULL
LW_Sy = FALSE
if (is.null(Sy)){
  ###
  Sy = t(Y) %*% Y /n
  if (LW_Sy){
    lw_cov <- corpcor::cov.shrink(Y)
    Sy <- as.matrix(lw_cov)
  }
}




svd_Sy = svd(Sy)
sqrt_inv_Sy = svd_Sy$u %*% diag(sapply(svd_Sy$d, function(x){ifelse(x > 1e-4, 
                                                                    1/sqrt(x), 0)}))  %*% t(svd_Sy$v)
tilde_Y = Y %*% sqrt_inv_Sy
Kx = NULL
if(!is.null(Kx)){
  Sx_tot = Sx + lambda_Kx * as.matrix(Kx) 
}else{
  Sx_tot = Sx
}

#lw_cov <- corpcor::cov.shrink(scale(t(as.matrix(newX_summ))))
#Sx <- as.matrix(lw_cov)

print("Using CVXR")
print("lambda is")
lambda = 0.1
print(lambda)

for (g in 1:length(groups)){
  print(g)
Sg = t(X[, groups[[g]]]) %*%X /n
write_csv(data.frame(Sg), paste0("~/Downloads/groups_Sg", g, ".csv" ))
}

### Use CVXR
if (solver=="CVXR"){
  B <- Variable(p, q)
  # Define the group lasso penalty
  group_lasso_penalty <- 0
  for (i in 1:length(groups)){
    # Assuming you have a way to define the rows in each group
    group_lasso_penalty <- group_lasso_penalty + norm(B[groups[[i]],], "2")
  }
  
  objective <- Minimize(norm(tilde_Y- X %*% B, 'F')^2 + 
                          lambda * group_lasso_penalty
  )
  
  
  problem <- Problem(objective)
  result <- solve(problem)
  B_opt <- result$getValue(B)
}else{
  # nb_groups = length(unique(groups$group.id[index_groups]))
  # groups
  # all_groups = sort(unique(groups$group.id[index_groups]))
  U = matrix(0, p, q)
  Z = matrix(0, p, q)
  B = matrix(0, p, q)
  
  prod_xy = t(X) %*% tilde_Y/n
  #invSx = solve(diag(sapply(1:g, function(gg){sqrt( length(groups[[gg]]))})) %*% Sx_tot %*% diag(sapply(1:g, function(gg){sqrt( length(groups[[gg]]))})) + rho *diag(rep(1, nrow(Sx_tot))))
  #### dispatch
 
  invSx0 = matrix(0, length(groups), p)
  rho = 10
  # for (gg in 1:length(groups)){
  #   for (ggg in 1:length(groups)){
  #     invSx0[gg, groups[[ggg]]] =  Sx[gg, ggg]/length(groups[[ggg]])
  #   }
  # }
  print(max(invSx0))
  niter = 100
  rho=10
  #Sx <- sapply(1:length(groups), function(g){ t(X[, groups[[g]]]) %*% X/n })
  for (i in 2:niter){
    Uold = U
    Zold = Z
    Bold = B
    for (g in 1:length(groups)){
      #print(paste0("group g = ", g))
      #eye = matrix(0, length(groups[[g]]), p)
      #eye[, groups[[g]]] =diag(1, length(groups[[g]]))
      #time1<-tic()
      #Sg <- read_csv2(paste0("~/Downloads/groups_Sg", g, ".csv" ))
      #time2<-toc()
      #print(time2)
      
      time1<-tic()
      Sg <- t(X[, groups[[g]]]) %*%X /n
      time2<-toc()
      print(time2)
      
      B[groups[[g]],] <- -1/rho^2 * invSx0[g,]  %*% (prod_xy  + rho* (Z - U))
      #B[groups[[g]],] <- -1/rho^2 *  (t(X[, groups[[g]]]) %*% X/n)   %*% (prod_xy  + rho* (Z - U))
      B[groups[[g]],] <- B[groups[[g]],]  + 1/rho * (prod_xy[groups[[g]],]   + rho* (Z - U)[groups[[g]],])
    }
    print(max(B))
    Z = B + U
    norm_col = sapply(1:length(groups), function(g){sqrt(sum(Z[groups[[g]],]^2))})
    for (g in 1:length(groups)){
      if(norm_col[g] < lambda * sqrt(length(groups[[g]]))){
        Z[groups[[g]],] = 0
        print(g)
      }else{
        print(c(g,  (1- (lambda * sqrt(length(groups[[g]])) /rho)/norm_col[g]) ))
        Z[groups[[g]],] =  (1- (lambda * sqrt(length(groups[[g]])) /rho)/norm_col[g]) * Z[groups[[g]],]
      }
      print(max( Z[groups[[g]],]))
    }
    U = U + B - Z
    print(c("ADMM iter", i, norm(Z - B)/sqrt(p), norm(Zold - Z)/sqrt(p), norm(Uold - U)/sqrt(p)))
    if (max(c(norm(Z - B), norm(Zold - Z)/sqrt(p))) <thresh){
      break
    }
  }
  B_opt = B
}

B_opt[which(abs(B_opt)<1e-5)] = 0
print(B_opt)

svd_Sx = svd(Sx)
sqrt_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, sqrt(x), 0)}))  %*% t(svd_Sx$u)
sqrt_inv_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sx$u)
sol = svd(sqrt_Sx %*% B_opt)
print(sqrt_Sx %*% B_opt)
#sol = svd(t(B_OLS[I, ]), nu = r, nv=r)
V = sqrt_inv_Sy %*% sol$v[, 1:r]
#U = matrix(0, p, r)
# B = U \tilde{V} 
inv_D = diag(sapply(1:r, FUN=function(x){ifelse(sol$d[x]<1e-4, 0, 1/sol$d[x])}))
U = B_opt %*% sol$v[, 1:r] %*% inv_D ### = U\lambda
print(t(U) %*% Sx %*% U)
print(t(V) %*% Sy %*% V)



loss = mean((Y %*% V - X %*% U)^2)


#atlas800 = readnii("/Users/cdonnat/Documents/group-CCA/experiments/real-data/fMRI-data/data/atlases/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz")
#atlas.d = atlas800@.Data
#atlas.d <- reshape2::melt(atlas.d)
#colnames(activations)
