library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
library(rrpack)
library(corpcor)
#setwd("~/Documents/group-CCA/")

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
source("elena/group_reduced_rank_regression.R")

#Simulation for missing values in both X and Y

args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
print(seed)
name_exp <- args[2]
set.seed(seed)
n <- as.numeric(args[3])
strength_theta <- args[4]
# <- as.numeric(300, 500, 700)
rs <- c(as.numeric(args[6]))
ps <- as.numeric(args[5])
overlaps <- c(0)
#props <- c(0, 0.1, 0.2)
props <- c(0)
noise = 1
seeds = 1:100
normalize_diagonal = TRUE
LW_Sy = TRUE
nnzero_values = c(5)

result = c()
for(seed_n in seeds){
  #for (n in c(100, 300, 500, 1000, 10000)){
  set.seed(seed * 100 + seed_n)
  print("Start loop")
  for(nnzeros in nnzero_values){
    for(p in c(ps)){
    #for (p in c(20, 50, 80, 100, 200, 500, 1000)){
      for (q in c(10, 30, 50, 70)){
      #for(nnzeros in c(5, 10, 15, 20, 50)){
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
        for (r_pca in c(5)){
          if ( (max(r_pca, r, nnzeros) < p) ) {
            for (overlapping_amount in overlaps){
              for(prop_missing in props){
                cat("seed:")
                cat(seed, " ")
                print(c(n, r, r_pca, strength_theta))
                start_time_creation <- system.time({
                gen = generate_example_group(n, p, q, 
                                  r_pca = r_pca,
                                  nnzeros = nnzeros,
                                  theta = thetas,
                                  lambda_pca = 1,
                                  r = r, overlapping_amount = overlapping_amount,
                                  normalize_diagonal = normalize_diagonal,
                                  n_new = 5000) 
              })
                print("Generation time")
              print(start_time_creation[[1]])
                print("We're done generating")
                X = gen$X
                Y = gen$Y
                groups = gen$groups
                Sigma0_sqrt = sqrtm(gen$Sigma)$B
                Sigma_hat_sqrt = sqrtm(gen$S)$B
                
                
                # for (lambda in c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 
                #                  0.25, 0.5, 1, 5, 7.5, 10)){
                #   tryCatch({
                #     #### if it's all zero then just stop
                #       start_time_alt <- system.time({
                #         alt <- CCA_rrr(X, Y, Sx=NULL, Sy=NULL,
                #                        lambda =lambda, Kx=NULL, r=r,
                #                        solver="CVXR", LW_Sy =  LW_Sy)
                #       })
                #       #init_coef = list(U = alt$U, V = alt$V)
                #       alt$U[which(is.na(alt$U))] <- 0
                #       alt$V[which(is.na(alt$V))] <- 0
                #       
                #       
                #       result = rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                #                                                  alt$U, alt$V, gen$u, 
                #                                                  gen$v,
                #                                                  Sigma_hat_sqrt = Sigma_hat_sqrt, 
                #                                                  Sigma0_sqrt = Sigma0_sqrt), 
                #                                         "noise" = noise, "method" = paste0("RRR-", lambda),
                #                                         "prop_missing" = prop_missing,
                #                                         "overlapping_amount" = overlapping_amount,
                #                                         "nnzeros" = nnzeros,
                #                                         "theta_strength" = strength_theta,
                #                                         "n" = n,
                #                                         "r_pca" = r_pca,
                #                                         "exp" = seed * 100 + seed_n,
                #                                         "normalize_diagonal" = normalize_diagonal,
                #                                         "lambda_opt" = lambda,
                #                                         "time" = start_time_alt[[1]]))
                #       
                #     
                #   }, error = function(e) {
                #     # Print the error message
                #     cat("Error occurred in Alt", lambda, ":", conditionMessage(e), "\n")
                #     # Skip to the next iteration
                #   })
                # 
                #   tryCatch({
                #       start_time_alt <- system.time({
                #         alt <- CCA_group_rrr(X, Y, 
                #                        groups = gen$groups,
                #                        Sx=NULL, Sy=NULL,
                #                        lambda =lambda, Kx=NULL, r=r,
                #                         solver="ADMM", LW_Sy =  LW_Sy,
                #                         scale = TRUE,
                #                         rho=1,
                #                         niter=1e3)
                #       })
                #       #init_coef = list(U = alt$U, V = alt$V)
                #       alt$U[which(is.na(alt$U))] <- 0
                #       alt$V[which(is.na(alt$V))] <- 0
                #       
                #       
                #       result = rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                #                                                  alt$U, alt$V, gen$u, 
                #                                                  gen$v,
                #                                                  Sigma_hat_sqrt = Sigma_hat_sqrt, 
                #                                                  Sigma0_sqrt = Sigma0_sqrt), 
                #                                         "noise" = noise, "method" = paste0("group-RRR-CVX-", lambda),
                #                                         "prop_missing" = prop_missing,
                #                                         "overlapping_amount" = overlapping_amount,
                #                                         "nnzeros" = nnzeros,
                #                                         "theta_strength" = strength_theta,
                #                                         "n" = n,
                #                                         "r_pca" = r_pca,
                #                                         "exp" = seed * 100 + seed_n,
                #                                         "normalize_diagonal" = normalize_diagonal,
                #                                         "lambda_opt" = lambda,
                #                                         "time" = start_time_alt[[1]]))
                #       
                #     
                #   }, error = function(e) {
                #     # Print the error message
                #     cat("Error occurred in group", lambda, ":", conditionMessage(e), "\n")
                #     # Skip to the next iteration
                #   })
                # }
                
                
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
                                                            Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                            Sigma0_sqrt = Sigma0_sqrt),
                                                   "noise" = noise,  method = "Oracle",  
                                                   "prop_missing" = prop_missing, 
                                                   "nnzeros" = nnzeros,
                                                   "theta_strength" = strength_theta,
                                                   "overlapping_amount" = overlapping_amount,
                                                   "r_pca" = r_pca,
                                                   "n" = n,
                                                   "normalize_diagonal" = normalize_diagonal,
                                                   "exp" = seed * 100 + seed_n,
                                                   "lambda_opt" = 0,
                                                   "time" = 0
                                                   
                )
                )
                
                print(paste0("Starting ", "Alt opt") )
                
                tryCatch({
                start_time_alt2 <- system.time({
                  res_alt = CCA_rrr.CV(X, Y, 
                             r=r, Kx = NULL, lambda_Kx = 0,
                             do.scale = TRUE,
                             param_lambda=c(10^seq(-3, 1, length.out = 30)),
                             kfolds=3, solver="rrr", LW_Sy = LW_Sy)
                })
                res_alt$ufinal[which(is.na( res_alt$ufinal))] <- 0
                res_alt$vfinal[which(is.na( res_alt$vfinal))] <- 0
                Uhat <- res_alt$ufinal[, 1:r]
                Vhat <- res_alt$vfinal[, 1:r]
                if (sum(apply(res_alt$ufinal, 1, function(x){sum(x!=0)}) >0) <r){
                  #### Choose another lambda
                  lambda_chosen = max(res_alt$resultsx$lambda[which(res_alt$resultsx$rmse > 1.05 * min(res_alt$resultsx$rmse))])
                  start_time_alt <- system.time({
                    res_alt <- CCA_rrr(X, Y, Sx = NULL, Sy=NULL,
                                       lambda =lambda_chosen, Kx=NULL, r, highdim=TRUE,
                                       solver="rrr",
                                       LW_Sy = LW_Sy, do.scale=TRUE)
                    res_alt$U[which(is.na(res_alt$U))] <- 0
                    res_alt$V[which(is.na(res_alt$v))] <- 0
                    Uhat <- res_alt$U[, 1:r]
                    Vhat <- res_alt$V[, 1:r]
                    
                  })
                }
                print(res_alt$ufinal)
                result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                            Uhat, 
                                                            res_alt$vfinal[, 1:r], 
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                            Sigma0_sqrt = Sigma0_sqrt),
                                                   "noise" = noise,  method = "RRR-opt",  
                                                   "prop_missing" = prop_missing, 
                                                   "nnzeros" = nnzeros,
                                                   "theta_strength" = strength_theta,
                                                   "overlapping_amount" = overlapping_amount,
                                                   "r_pca" = r_pca,
                                                   "n" = n,
                                                   "exp" = seed * 100 + seed_n,
                                                   "normalize_diagonal" = normalize_diagonal,
                                                   "lambda_opt" = res_alt$lambda,
                                                   "time" = start_time_alt2[[1]]
                )
                )
                }, error = function(e) {
                  # Print the error message
                  cat("Error occurred in RRR CV:", conditionMessage(e), "\n")
                  # Skip to the next iteration
                })
              

                print(paste0("Starting ", "Alt opt2") )
                tryCatch({
                start_time_alt3 <- system.time({
                  res_alt = CCA_rrr.CV(X, Y,
                             r=r, Kx = NULL, lambda_Kx = 0,
                             param_lambda=c(10^seq(-3, 1, length.out = 30)),
                             kfolds=3, solver="ADMM", LW_Sy = LW_Sy, do.scale = TRUE,
                             rho=1, niter=2 * 1e4, thresh = 1e-6)
                })
                res_alt$ufinal[which(is.na( res_alt$ufinal))] <- 0
                res_alt$vfinal[which(is.na( res_alt$vfinal))] <- 0
                Uhat <- res_alt$ufinal[, 1:r]
                Vhat <- res_alt$vfinal[, 1:r]
                lambda_chosen = res_alt$lambda
                if (sum(apply(res_alt$ufinal, 1, function(x){sum(x!=0)}) >0) <r){
                  #### Choose another lambda
                  while(sum(apply(Uhat, 1, function(x){sum(x!=0)}) >0) <r){
                    lambda_chosen = lambda_chosen / 2
                    #start_time_alt <- system.time({
                      res_alt <- CCA_rrr(X, Y, Sx = NULL, Sy=NULL,
                                         lambda =lambda_chosen, Kx=NULL, 
                                         r=r, 
                                         highdim=TRUE,
                                         solver="ADMM",
                                         LW_Sy = LW_Sy, do.scale=TRUE,
                                         thresh = 1e-6)
                      res_alt$U[which(is.na(res_alt$U))] <- 0
                      res_alt$V[which(is.na(res_alt$v))] <- 0
                      Uhat <- res_alt$U[, 1:r]
                      Vhat <- res_alt$V[, 1:r]
                      
                    #})
                    
                  }
                }
                result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew,
                                                            Uhat,
                                                            Vhat[, 1:r],
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt = Sigma_hat_sqrt,
                                                            Sigma0_sqrt = Sigma0_sqrt),
                                                   "noise" = noise,
                                                   method = "RRR-ADMM-opt",
                                                   "prop_missing" = prop_missing,
                                                   "nnzeros" = nnzeros,
                                                   "theta_strength" = strength_theta,
                                                   "overlapping_amount" = overlapping_amount,
                                                   "r_pca" = r_pca,
                                                   "n" = n,
                                                   "exp" = seed * 100 + seed_n,
                                                   "normalize_diagonal" = normalize_diagonal,
                                                   "lambda_opt" = lambda_chosen,
                                                   "time" = start_time_alt3[[1]]
                )
                )
            }, error = function(e) {
              # Print the error message
              cat("Error occurred in CVXR CV:", conditionMessage(e), "\n")
              # Skip to the next iteration
            })
          
                
          print(paste0("Starting ", "CV opt-group") )

              tryCatch({
                start_time_alt4<- system.time({
                  res_alt = CCA_group_rrr.CV(X, Y, 
                             groups = gen$groups,
                             r=r, Kx = NULL, lambda_Kx = 0,
                             param_lambda=c(10^seq(-3, 1, length.out = 30)),
                             kfolds=3, solver="ADMM", LW_Sy = LW_Sy, do.scale=TRUE,
                             rho=1, niter=2 * 1e4,
                             thresh=1e-6)
                })
                res_alt$ufinal[which(is.na( res_alt$ufinal))] <- 0
                res_alt$vfinal[which(is.na( res_alt$vfinal))] <- 0
                Uhat <- res_alt$ufinal[, 1:r]
                Vhat <- res_alt$vfinal[, 1:r]
                lambda_chosen = res_alt$lambda
                
                if (sum(apply(res_alt$ufinal, 1, function(x){sum(x!=0)}) >0) <r){
                  #### Choose another lambda
                  while(sum(apply(Uhat, 1, function(x){sum(x!=0)}) >0) <r){
                    lambda_chosen = lambda_chosen / 2
                    #start_time_alt <- system.time({
                    res_alt <- CCA_group_rrr(X, Y, 
                                             groups = gen$groups,
                                             Sx = NULL, Sy=NULL,
                                       lambda =lambda_chosen, Kx=NULL, 
                                       r=r, 
                                       solver="ADMM",
                                       LW_Sy = LW_Sy, do.scale=TRUE,
                                       rho=1, niter=2 * 1e4,
                                       thresh=1e-6)
                    res_alt$U[which(is.na(res_alt$U))] <- 0
                    res_alt$V[which(is.na(res_alt$v))] <- 0
                    Uhat <- res_alt$U[, 1:r]
                    Vhat <- res_alt$V[, 1:r]
                    
                    #})
                    
                  }
                }

                print(res_alt$ufinal)
                result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                            Uhat, 
                                                            Vhat[, 1:r], 
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                            Sigma0_sqrt = Sigma0_sqrt),
                                                   "noise" = noise,  
                                                   method = "CVX-opt-group",  
                                                   "prop_missing" = prop_missing, 
                                                   "nnzeros" = nnzeros,
                                                   "theta_strength" = strength_theta,
                                                   "overlapping_amount" = overlapping_amount,
                                                   "r_pca" = r_pca,
                                                   "n" = n,
                                                   "exp" = seed * 100 + seed_n,
                                                   "normalize_diagonal" = normalize_diagonal,
                                                   "lambda_opt" = lambda_chosen,
                                                   "time" = start_time_alt4[[1]]
                )
                )
          }, error = function(e) {
            # Print the error message
            cat("Error occurred in group CV", ":", conditionMessage(e), "\n")
            # Skip to the next iteration
          })
        

          for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                                 "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                                 "SCCA_Parkhomenko")){
                  
                  print(paste0("Starting ", method))
                  
                  
                  tryCatch({
                    start_time_additional_method <- system.time({
                      test1<-additional_checks(gen$X,
                                               gen$Y, S=NULL, 
                                               rank=r, kfolds=3, 
                                               method.type = method,
                                               lambdax= 10^seq(-3,1, length.out = 30),
                                               lambday = c(0))
                    })
                    result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                                test1$u[, 1:r], 
                                                                test1$v[, 1:r], 
                                                                gen$u, gen$v,
                                                                Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                                Sigma0_sqrt = Sigma0_sqrt),
                                                       "noise" = noise,  method = method,  
                                                       "prop_missing" = prop_missing, 
                                                       "overlapping_amount" = overlapping_amount,
                                                       "nnzeros" = nnzeros,
                                                       "theta_strength" = strength_theta,
                                                       "r_pca" = r_pca,
                                                       "n" = n,
                                                       "exp" = seed * 100 + seed_n,
                                                       "normalize_diagonal" = normalize_diagonal,
                                                       "lambda_opt" = 0,
                                                       "time" = start_time_additional_method[[1]]
                    )
                    )
                  }, error = function(e) {
                    # Print the error message
                    cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
                    # Skip to the next iteration
                  })
                }
                write_csv(result, paste0("elena/missing/results/2024-group-newest_RRR_efficient_results", name_exp, ".csv"))
                print("Done loop")
              
                #write.csv(result, "missing/simulation-RRR-results-sparse.csv", row.names = F)
              }
              print("Done inner most loop")
            }
          }
        }
      }
    }
  }

  }
}

# name_exp = "NA"
# result <- rbind(read_csv(paste0("elena/missing/results/2024-group-newest_RRR_efficient_resultsexp_new_high.csv")),
#                 read_csv(paste0("elena/missing/results/2024-group-newest_RRR_efficient_resultsexp_new_high_q30.csv")),
#                 read_csv(paste0("elena/missing/results/2024-group-newest_RRR_efficient_resultsexp_new_medium.csv")),
#                 read_csv(paste0("elena/missing/results/2024-group-newest_RRR_efficient_resultsexp_new_medium_q30.csv")),
#                 read_csv(paste0("elena/missing/results/2024-group-newest_RRR_efficient_resultsexp_new_low.csv")),
#                 read_csv(paste0("elena/missing/results/2024-group-newest_RRR_efficient_resultsexp_new_low_q30.csv")))
# unique(result$method)
# res = result %>% 
#   group_by(n, p1, p2, nnzeros, r, r_pca, method) %>%
#   summarise_if(is.numeric,mean)
# 
# unique(res$method)
# legend_order <- c("Oracle",  "FIT_SAR_CV", 
#                   "FIT_SAR_BIC", "Witten_Perm", "Witten.CV",
#                   "SCCA_Parkhomenko", "Waaijenborg-CV", "Waaijenborg-Author",
#                   #"RRR-0.5" ,"RRR-7.5","RRR-10","RRR-12.5",  "RRR-20",   
#                   "RRR-opt",    "RRR-ADMM-opt", "CVX-opt-group"   )
# my_colors <- c( "black", "red", "indianred4",
#                 "orange", "yellow", "chartreuse2",
#                 "burlywood2", "burlywood4",
#                 # "lightblue", "lightblue3","cyan", "dodgerblue", "dodgerblue4", 
#                 "navyblue", "cyan", "dodgerblue")
# 
# labels_n <-    c("Oracle",  "SAR CV (Wilms et al)", 
#                  "SAR BIC (Wilms et al)", 
#                  "Sparse CCA, permuted\n(Witten et al)", 
#                  "Sparse CCA with CV\n(Witten et al)",
#                  "SCCA (Parkhomenko et al)", "Sparse CCA with CV\n(Waaijenborg et al)",
#                  "Sparse CCA(Waaijenborg et al)",
#                  # "RRR-0.5" ,"RRR-7.5","RRR-10","RRR-12.5",  "RRR-20",   
#                  "RRR-CCA (this paper)",   "RRR-CCA 2(this paper)",
#                  "RRR-CCA-group (this paper)")
# theme_set(theme_bw(base_size = 14))
# colnames(res)
# 
# result %>% filter(method == "Oracle")
# ggplot(res,
#        aes(x=p1, 
#            y =distance_tot, 
#            colour =method)) +
#   geom_point()+
#   geom_line()+
#   scale_y_log10()+
#   scale_color_manual(values = my_colors, breaks = legend_order,
#                      labels = labels_n) 
# 
# ggplot(result,
# aes(x=p1, 
#     y = distance_tot, 
#     colour =method)) +
#   geom_point()+
#   geom_line()+
#   scale_color_manual(values = my_colors, breaks = legend_order,
#                      labels = labels_n) +
#   facet_grid(theta_strength~p2)
