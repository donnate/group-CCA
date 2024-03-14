library(tidyverse)
library(tidyverse)
library(igraph)
library(ggplot2)
library(vroom)
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

args <- commandArgs(trailingOnly=TRUE)
lambda <- as.numeric(args[1])
test_fold <- as.numeric(args[2])
val_fold <- ifelse(test_fold < 16, test_fold + 1, 1)

X <- (vroom("data/activations_X_preprocessed.csv"))
print("Done loading X")
behaviour <- read_csv("data/activations_Y_preprocessed.csv")
group_assignment <- readxl::read_xlsx("~/Downloads/activation_groups.xlsx", col_names = FALSE)
folds = t(read_csv("data/folds.csv"))

colnames(group_assignment) <- "group.id"
index_groups = which(group_assignment$group.id!=0)
groups <- sapply(1:length(unique(group_assignment$group.id[index_groups])),
                 function(g){which(group_assignment$group.id[index_groups] == g)})

Y = as.matrix(behaviour)
n = nrow(X)
p = ncol(X)
q = ncol(Y)
Sy = NULL
LW_Sy = FALSE
rho = 1
niter = 1000
r = 2
thresh = 1e-5
##### Split into different 
index_test = as.numeric(folds[test_fold,])
index_val = as.numeric(folds[val_fold,])
index_train = (1:n)[-c(index_test, index_val)]


if (is.null(Sy)){
  ###
  Sy = t(Y[index_train,]) %*% Y[index_train,] /n
  if (LW_Sy){
    lw_cov <- corpcor::cov.shrink(Y)
    Sy <- as.matrix(lw_cov)
  }
}



svd_Sy = svd(Sy)
sqrt_inv_Sy = svd_Sy$u %*% diag(sapply(svd_Sy$d, function(x){ifelse(x > 1e-4, 
                                                                    1/sqrt(x), 0)}))  %*% t(svd_Sy$v)
tilde_Y = Y %*% sqrt_inv_Sy




print("Using CVXR")
print("lambda is")
print(lambda)

# large_groups = c(19, 18, 15, 13, 14,11, 12)
# for (g in large_groups){
#  print(g)
#  Sg = t(X[index_train, groups[[g]]]) %*%X[index_train, ] /n
#  write_csv(data.frame(Sg), paste0("~/Downloads/groups_Sg", g, ".csv" ))
# }


#STOP
svd_X = svd(X[index_train, ]/sqrt(n))
print("Done with the SVD")
svd_left = (svd_X$v) %*% diag(1/ (svd_X$d^2 + rho))
print("Done with svd_left")
{

  U = matrix(0, p, q)
  Z = matrix(0, p, q)
  B = matrix(0, p, q)
  
  prod_xy = t(X[index_train,]) %*% tilde_Y[index_train,]/length(index_train)
  print("Starting iter")
  print(max(prod_xy))
  for (i in 1:niter){
    Uold = U
    Zold = Z
    Bold = B
    time1<-tic()
    B = svd_left  %*% (t(svd_X$v) %*% (prod_xy  + rho * (Z - U)))
    time2<-toc()
    print(paste0("iteration ", i, " time = ", time2))
    print(max(B))
    Z = B + U
    norm_col = sapply(1:length(groups), function(g){sqrt(sum(Z[groups[[g]],]^2))})
    for (g in 1:length(groups)){
      if(norm_col[g] < lambda * sqrt(length(groups[[g]]))/rho){
        Z[groups[[g]],] = 0
        print(g)
      }else{
        #print(c(g,  (1- (lambda * sqrt(length(groups[[g]])) /rho)/norm_col[g]) ))
        Z[groups[[g]],] =  (1- (lambda * sqrt(length(groups[[g]])) /rho)/norm_col[g]) * Z[groups[[g]],]
      }
      print(paste("Max of group ", g, " is ", max( Z[groups[[g]],])))
    }
    U = U + B - Z
    print(c("ADMM iter", i, norm(Z - B)/sqrt(p), norm(Zold - Z)/sqrt(p), norm(Uold - U)/sqrt(p)))
    if (max(c(norm(Z - B), norm(Zold - Z)/sqrt(p))) <thresh){
      break
    }
  }
  B_opt = B
}

### groups 3, 4, 7,, 8, 10 , 11, 15, 18,19 should be computed ahead
B_opt[which(abs(B_opt)<1e-5)] = 0
print(B_opt)

sol = svd(((svd_X$v) %*% diag(sapply(svd_X$d^2, function(x){ifelse(x > 1e-4, sqrt(x), 0)})))  %*% (t(svd_X$v)  %*% B_opt))
#sol = svd(t(B_OLS[I, ]), nu = r, nv=r)
V = sqrt_inv_Sy %*% sol$v[, 1:r]
inv_D = diag(sapply(1:r, FUN=function(x){ifelse(sol$d[x]<1e-4, 0, 1/sol$d[x])}))
U = B_opt %*% sol$v[, 1:r] %*% inv_D ### = U\lambda

correlation <-
  data.frame("CCA_graph_rrr",
             lambda,
              test_fold,
             "val_fold" = val_fold,
             diag(cov(as.matrix(X)[index_train,] %*% U,
                      as.matrix(Y)[index_train,] %*%  V)),
             apply(((as.matrix(X)[index_train,] %*% U) -
                      (as.matrix(Y)[index_train,] %*%  V))^2, 2, mean),
             diag(t(as.matrix(X)[index_test,] %*% U) %*%
                    as.matrix(Y)[index_test,] %*% V),
             apply(((as.matrix(X)[index_test,] %*% U) -
                      (as.matrix(Y)[index_test,] %*%  V))^2, 2, mean),
             diag(cor(as.matrix(X)[index_test,] %*% U,
                    as.matrix(Y)[index_test,] %*% V)),
             diag(t(as.matrix(X)[index_val,] %*% U) %*%
                    as.matrix(Y)[index_val,] %*% V),
             apply(((as.matrix(X)[index_val,] %*% U) -
                      (as.matrix(Y)[index_val,] %*%  V))^2, 2, mean),
             diag(cor(as.matrix(X)[index_val,] %*% U,
                  as.matrix(Y)[index_val,] %*% V))
  )

write_csv(correlation, file = paste0("data/results_l", lambda, "_test_fold", test_fold, ".csv"))


#### Follow up: load results:

file_list <- list.files(path = "~/Documents/group-CCA/data", 
                         pattern = "results*", full.names = TRUE)
results <- bind_rows(lapply(file_list, read.csv))
colnames(results) <- c("Method", "Lambda",
                       "Test_fold","Val_fold",
                       "train_cov",
                       "train_mse",
                       "test_prod",
                       "test_mse",
                        "test_cor",
                       "val_prod",
                       "val_mse",
                       "val_cor")
results2 = pivot_longer(results, cols = -c("Method", "Lambda",
                                            "Test_fold","Val_fold")
                       )
results2 = pivot_wider(results2, id_cols = c("Method", "Lambda",
                                            "Test_fold","Val_fold"),
                       names_from = "name",
                       values_from = "value")

summary_correlation = results %>% 
  group_by(Method, Lambda, Test_fold, Val_fold) %>%
  summarise(test_mse = sum(test_mse),
         train_mse = mean(train_mse),
         val_mse =  mean(val_mse),
         test_cor = mean(test_cor),
         val_cor = mean(val_cor)) %>%
  group_by(Method, Lambda) %>% summarise_all(mean, na.rm=TRUE)%>%
  arrange((test_mse)) %>% ungroup()

ordered_res = summary_correlation %>% 
  filter(lambda< 0.1)%>% ungroup()%>% arrange(test_mse)
lambda_opt = 0.0001


relevant_correlations = summary_correlation %>% 
  filter( Lambda == lambda_opt ) %>%
  dplyr::select(Method, test_mse, test_cor, val_mse, val_cor)


ggplot(summary_correlation%>% ungroup())+
  geom_line(aes(x = as.numeric(Lambda), y=test_mse))+
  geom_point(aes(x = as.numeric(Lambda), y=test_mse))+
  scale_x_log10()

Vhat_comp = V
rownames(Vhat_comp) = colnames(Y)
print(Vhat_comp)
heatmap(Vhat_comp)
df_V = data.frame(Vhat_comp)
colnames(df_V) = c("CD-1", "CD-2")
df_V["question"] = colnames(Y)
df_V = pivot_longer(df_V, cols=-c("question"))



# Example new labels
new_labels <- c("Drive" , 
                "Funseeking",
                "Reward Response" ,
                "Total",
                "Distress",
                "Anhedonia",
                "Anxarousal",
                "PANAS Positive Affect",
                "PANAS Negative Affect")

# Assuming your questions are in a column named 'question'
old_labels <- unique(df_V$question)
label_mapping <- setNames(new_labels, old_labels)
df_V <- df_V %>% mutate(question_new = label_mapping[question])
theme_set(theme_bw(base_size = 18))
ggplot(df_V, aes(x = value, y = reorder(question_new, value), fill=name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ name) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Question", x = "Loading Intensity") + 
  labs(fill = "Canonical\nDirection")
