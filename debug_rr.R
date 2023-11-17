library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
setwd("~/Documents/group-CCA/")

source("elena/generate_example_rrr.R")
source('experiments/sparse_CCA/experiment_functions.R')
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source("elena/missing/evaluation.R")
#source("elena/missing/original_CCA_impute.R")
source("elena/gradient_descent.r")
#source("elena/iterative_cca.R")
source("elena/reduced_rank_regression.R")


#Simulation for missing values in both X and Y

args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
print(seed)
name_exp <- args[2]
set.seed(seed)
#n <- as.numeric(args[3])
strength_theta <- args[3]
#p_val <- as.numeric(args[5])
overlaps <- c(0, 1)
#props <- c(0, 0.1, 0.2)
props <- c(0)
noise = 1
seeds = 1:100
normalize_diagonal = TRUE
result = c()
for(seed_n in seeds){
  for (n in c(100, 300,500, 1000, 10000)){
  set.seed(seed * 100 + seed_n)
  for(nnzeros in c(20, 10, 15, 50, 5)){
    #for(p in c(100,  200, 300,  500, 800, 80, 20)){
    for (p in c(20, 50, 80, 100, 200, 500, 1000)){
      for (q in c(10, 20, 30, 50, 80)){
      #for(nnzeros in c(5, 10, 15, 20, 50)){
      rs = ifelse( p <6, c(2,  5), c(2, 5, 10))
      for (r in rs){
        
        if ( strength_theta == "high"){
          thetas <- diag(seq(0.9, 0.75, length.out = r))
        }else{
          if ( strength_theta == "medium"){
            thetas <- diag(seq(0.7, 0.55, length.out = r))
          }
          else{
            thetas <- diag(seq(0.5, 0.35, length.out = r))
          }
        }
        for (r_pca in c(0, 5)){
          if (max(r_pca * nnzeros, r * nnzeros) < p) {
            for (overlapping_amount in overlaps){
              for(prop_missing in props){
                cat("seed:")
                cat(seed, " ")
                print(c(n, r, r_pca, strength_theta))
                #gen = generate(n, p, q, s, prop_missing)
                gen = generate_example_sparse_U(n, p, q,
                                                r_pca = r_pca,
                                                nnzeros = nnzeros,
                                                theta = thetas,
                                                lambda_pca = 1,
                                                r = r,
                                                overlapping_amount = overlapping_amount,
                                                normalize_diagonal = normalize_diagonal) 
                X = gen$X
                Y = gen$Y
                #Xna = gen$Xna
                #Yna = gen$Yna
                Sigma0_sqrt_inv = sqrtm(gen$Sigma)$Binv
                Sigma_hat_sqrt_inv = sqrtm(gen$S  + 1e-4 *diag(nrow=nrow(gen$S)))$Binv
                
                
                if (p < n){
                  
                  start_time_rrr <- system.time({
                    rrr <- CCA_rrr(X, Y, lambda =0, Kx=NULL, r, highdim=FALSE)
                  })
                  
                  result = rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, rrr$U, rrr$V, gen$u, 
                                                             gen$v,
                                                             Sigma_hat_sqrt_inv = Sigma_hat_sqrt_inv, 
                                                             Sigma0_sqrt_inv = Sigma0_sqrt_inv), 
                                                    "noise" = noise, "method" = "RRR", 
                                                    "prop_missing" = prop_missing,
                                                    "overlapping_amount" = overlapping_amount,
                                                    "nnzeros" = nnzeros,
                                                    "theta_strength" = strength_theta,
                                                    "n" = n,
                                                    "r_pca" = r_pca,
                                                    "exp" = seed * 100 + seed_n,
                                                    "normalize_diagonal" = normalize_diagonal,
                                                    "time" = start_time_rrr[[1]]))
                }
                
                for (lambda in c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 
                                 0.25, 0.5, 1)){
                  if (lambda == 0.0001){
                    init_coef = NULL
                  }
                  tryCatch({
                    #### if it's all zero then just stop
                    if (is.null(init_coef) || ((norm(init_coef$U, "F") > 1e-5) & (norm(init_coef$V, "F") > 1e-5))){
                      start_time_alt <- system.time({
                        alt <- CCA_rrr(X, Y, lambda =lambda, Kx=NULL, r, highdim=TRUE, penalty = "l21")
                      })
                      init_coef = list(U = alt$U, V = alt$V)
                      result = rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                                 alt$U, alt$V, gen$u, 
                                                                 gen$v,
                                                                 Sigma_hat_sqrt_inv = Sigma_hat_sqrt_inv, 
                                                                 Sigma0_sqrt_inv = Sigma0_sqrt_inv), 
                                                        "noise" = noise, "method" = paste0("RRR-", lambda),
                                                        "prop_missing" = prop_missing,
                                                        "overlapping_amount" = overlapping_amount,
                                                        "nnzeros" = nnzeros,
                                                        "theta_strength" = strength_theta,
                                                        "n" = n,
                                                        "r_pca" = r_pca,
                                                        "exp" = seed * 100 + seed_n,
                                                        "normalize_diagonal" = normalize_diagonal,
                                                        "time" = start_time_alt[[1]]))
                      
                    }
                  }, error = function(e) {
                    # Print the error message
                    cat("Error occurred in Alt", lambda, ":", conditionMessage(e), "\n")
                    # Skip to the next iteration
                  })
                }
                
                if ( p < n && prop_missing > 0){
                  Ximp = data.frame(gen$Xna) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE)))
                  Yimp = data.frame(gen$Yna) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE)))
                  
                  Ximp2 = data.frame(gen$Xna) %>% mutate_all(~replace_na(., median(., na.rm = TRUE)))
                  Yimp2 = data.frame(gen$Yna) %>% mutate_all(~replace_na(., median(., na.rm = TRUE)))
                  
                  start_time_cca <- system.time({
                    cca = CCA::cc(Ximp, Yimp)
                  })
                  
                  result = rbind(result, data.frame(evaluate(gen$newX, gen$newY, cca$xcoef[, 1:r], 
                                                             cca$ycoef[, 1:r], 
                                                             gen$u, gen$v,
                                                             Sigma_hat_sqrt_inv = Sigma_hat_sqrt_inv, 
                                                             Sigma0_sqrt_inv = Sigma0_sqrt_inv),
                                                    "noise" = noise,  method = "CCA-mean",  
                                                    "prop_missing" = prop_missing, 
                                                    "overlapping_amount" = overlapping_amount,
                                                    "nnzeros" = nnzeros,
                                                    "theta_strength" = strength_theta,
                                                    "n" = n,
                                                    "r_pca" = r_pca,
                                                    "exp" = seed * 100 + seed_n,
                                                    "normalize_diagonal" = normalize_diagonal,
                                                    "time" = start_time_cca[[1]]))
                  
                  start_time_cca2 <- system.time({
                    cca = CCA::cc(Ximp2, Yimp2)
                  })
                  result = rbind(result, data.frame(evaluate(gen$newX, gen$newY, cca$xcoef[, 1:r], 
                                                             cca$ycoef[, 1:r], 
                                                             gen$u, gen$v,
                                                             Sigma_hat_sqrt_inv = Sigma_hat_sqrt_inv, 
                                                             Sigma0_sqrt_inv = Sigma0_sqrt_inv),
                                                    "noise" = noise,  method = "CCA-median",  
                                                    "prop_missing" = prop_missing, 
                                                    "nnzeros" = nnzeros,
                                                    "theta_strength" = strength_theta,
                                                    "overlapping_amount" = overlapping_amount,
                                                    "r_pca" = r_pca,
                                                    "n" = n,
                                                    "exp" = seed * 100 + seed_n,
                                                    "normalize_diagonal" = normalize_diagonal,
                                                    "time" = start_time_cca2[[1]]))
                  
                  
                }
                
                #### Try out alternative approaches
                #### Oracle
                print("beginning oracle")
                set_u =  which(apply(gen$u,1, norm)>0)
                set_v =  which(apply(gen$v,1, norm)>0)
                t=CCA::cc(as.matrix(gen$X[,set_u]), as.matrix(gen$Y[, set_v]))
                Uhat = matrix(0, p, r)
                Vhat = matrix(0, q, r)
                Uhat[set_u, ] <-  t$xcoef[, 1:r]
                Vhat[set_v, ] <-  t$ycoef[, 1:r]
                result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, Uhat, 
                                                            Vhat, 
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt_inv = Sigma_hat_sqrt_inv, 
                                                            Sigma0_sqrt_inv = Sigma0_sqrt_inv),
                                                   "noise" = noise,  method = "Oracle",  
                                                   "prop_missing" = prop_missing, 
                                                   "nnzeros" = nnzeros,
                                                   "theta_strength" = strength_theta,
                                                   "overlapping_amount" = overlapping_amount,
                                                   "r_pca" = r_pca,
                                                   "n" = n,
                                                   "normalize_diagonal" = normalize_diagonal,
                                                   "exp" = seed * 100 + seed_n,
                                                   "time" = 0
                                                   
                )
                )
                
                start_time_alt2 <- system.time({
                  res_alt = CCA_rrr.CV(X, Y, 
                             r=r, Kx = NULL, lambda_Kx = 0,
                             param_lambda=10^seq(-4, 2, 0.25),
                             kfolds=10, penalty="l21")
                })
                
                print(paste0("Starting ", "Alt opt") )
                result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                            res_alt$ufinal[, 1:r], 
                                                            res_alt$vfinal[, 1:r], 
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt_inv = Sigma_hat_sqrt_inv, 
                                                            Sigma0_sqrt_inv = Sigma0_sqrt_inv),
                                                   "noise" = noise,  method = "Alt-opt",  
                                                   "prop_missing" = prop_missing, 
                                                   "nnzeros" = nnzeros,
                                                   "theta_strength" = strength_theta,
                                                   "overlapping_amount" = overlapping_amount,
                                                   "r_pca" = r_pca,
                                                   "n" = n,
                                                   "exp" = seed * 100 + seed_n,
                                                   "normalize_diagonal" = normalize_diagonal,
                                                   "time" = start_time_alt2[[1]]
                )
                )
                
                for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                                 "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                                 "SCCA_Parkhomenko", "Canonical Ridge-Author")){
                  
                  print(paste0("Starting ", method))
                  
                  
                  tryCatch({
                    start_time_additional_method <- system.time({
                      test1<-additional_checks(gen$X,
                                               gen$Y, S=NULL, 
                                               rank=r, kfolds=5, 
                                               method.type = method)
                    })
                    result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, test1$u[, 1:r], 
                                                                test1$v[, 1:r], 
                                                                gen$u, gen$v,
                                                                Sigma_hat_sqrt_inv = Sigma_hat_sqrt_inv, 
                                                                Sigma0_sqrt_inv = Sigma0_sqrt_inv),
                                                       "noise" = noise,  method = method,  
                                                       "prop_missing" = prop_missing, 
                                                       "overlapping_amount" = overlapping_amount,
                                                       "nnzeros" = nnzeros,
                                                       "theta_strength" = strength_theta,
                                                       "r_pca" = r_pca,
                                                       "n" = n,
                                                       "exp" = seed * 100 + seed_n,
                                                       "normalize_diagonal" = normalize_diagonal,
                                                       "time" = start_time_additional_method[[1]]
                    )
                    )
                  }, error = function(e) {
                    # Print the error message
                    cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
                    # Skip to the next iteration
                  })
                }
                
                write_csv(result, paste0("elena/missing/results/new_RRR_results", name_exp, ".csv"))
                #write.csv(result, "missing/simulation-RRR-results-sparse.csv", row.names = F)
              }
            }
          }
        }
      }
      }
    }
  }
}
}




