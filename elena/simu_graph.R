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
source("elena/graph_reduced_rank_regression.R")

#Simulation for missing values in both X and Y

args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
print(seed)
name_exp <- args[2]
set.seed(seed)
n <- as.numeric(args[3])
strength_theta <- args[4]
#p <- as.numeric(300, 500, 700)
rs <- c(as.numeric(args[5]))
p_val <- as.numeric(args[6])
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
    for(p in p_val){
      #for (p in c(20, 50, 80, 100, 200, 500, 1000)){
      for (q in c(10, 30, 50, 80)){
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
                    gen = generate_example_graph(n,
                                                 p1=p,
                                                 type_graph="2d-grid",
                                                 p2=q, order = 3,
                                                 r_pca = r_pca,
                                                 nnzeros = nnzeros,
                                                 do_plot = FALSE,
                                                 theta = thetas,
                                                 lambda_pca = 1,
                                                 nnzeros_pca = 20,
                                                 r = r, 
                                                 overlapping_amount = overlapping_amount,
                                                 normalize_diagonal = normalize_diagonal,
                                                 gen.using.gamma = TRUE,
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
                  Gamma_dagger = pinv(gen$Gamma)
                  
                  
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
                  tryCatch({
                  set_u =  which(apply(gen$u,1, norm)>0)
                  set_v =  which(apply(gen$v,1, norm)>0)
                  t=CCA::cc(as.matrix(gen$X[,set_u]), as.matrix(gen$Y[, set_v]))
                  Uhat = matrix(0, dim(gen$u)[1], r)
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
                }, error = function(e) {
                  # Print the error message
                  cat("Error occurred in CVXR CV:", conditionMessage(e), "\n")
                  # Skip to the next iteration
                })
                  
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
                                                       "time" = start_time_alt2[[4]]
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
                                           rho=1, niter=3 * 1e4, thresh = 1e-6)
                    })
                    res_alt$ufinal[which(is.na( res_alt$ufinal))] <- 0
                    res_alt$vfinal[which(is.na( res_alt$vfinal))] <- 0
                    Uhat <- res_alt$ufinal[, 1:r]
                    Vhat <- res_alt$vfinal[, 1:r]
                    
                    
                    print(res_alt$ufinal)
                    result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew,
                                                                Uhat,
                                                                res_alt$vfinal[, 1:r],
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
                                                       "lambda_opt" = res_alt$lambda,
                                                       "time" = start_time_alt3[[4]]
                    )
                    )
                  }, error = function(e) {
                    # Print the error message
                    cat("Error occurred in CVXR CV:", conditionMessage(e), "\n")
                    # Skip to the next iteration
                  })
                  
                  
                  print(paste0("Starting ", "CV opt-graph") )
                  
                  
                 
                  tryCatch({
                    start_time_alt4<- system.time({
                      res_alt = CCA_graph_rrr.CV(X, Y, 
                                                 Gamma = gen$Gamma,
                                                 r=r, Kx = NULL, lambda_Kx = 0,
                                                 param_lambda=c(10^seq(-3, 1, length.out = 30)),
                                                 kfolds=5, 
                                                 LW_Sy = LW_Sy, do.scale=TRUE,
                                                 rho=1, niter=2 * 1e4,
                                                 thresh=1e-4,
                                                 Gamma_dagger = Gamma_dagger)
                    })
                    res_alt$ufinal[which(is.na( res_alt$ufinal))] <- 0
                    res_alt$vfinal[which(is.na( res_alt$vfinal))] <- 0
                    Uhat <- res_alt$ufinal[, 1:r]
                    Vhat <- res_alt$vfinal[, 1:r]
                    
                    print(res_alt$ufinal)
                    result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                                Uhat, 
                                                                res_alt$vfinal[, 1:r], 
                                                                gen$u, gen$v,
                                                                Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                                Sigma0_sqrt = Sigma0_sqrt),
                                                       "noise" = noise,  
                                                       method = "CVX-opt-graph",  
                                                       "prop_missing" = prop_missing, 
                                                       "nnzeros" = nnzeros,
                                                       "theta_strength" = strength_theta,
                                                       "overlapping_amount" = overlapping_amount,
                                                       "r_pca" = r_pca,
                                                       "n" = n,
                                                       "exp" = seed * 100 + seed_n,
                                                       "normalize_diagonal" = normalize_diagonal,
                                                       "lambda_opt" = res_alt$lambda,
                                                       "time" = start_time_alt4[[4]]
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
                                                         "time" = start_time_additional_method[[4]]
                      )
                      )
                    }, error = function(e) {
                      # Print the error message
                      cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
                      # Skip to the next iteration
                    })
                  }
                  write_csv(result, paste0("elena/missing/results/2024-graph-newest_RRR_efficient_results", name_exp, ".csv"))
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



# df <- read_csv("elena/missing/results/2024-graph-newest_RRR_efficient_resultstest.csv")
# summ = df %>% group_by(n, p1, p2, r, r_pca,
#                             nnzeros, 
#                             overlapping_amount, noise, 
#                             #lambda_opt,
#                             method,
#                             theta_strength,
#                             normalize_diagonal,
#                             prop_missing) %>% 
#   summarise(distance_tot_mean = mean(distance_tot),
#             distance_U_mean = mean(distance_U),
#             distance_V_mean = mean(distance_V),
#             distance_tot_q50 = quantile(distance_tot, 0.5, na.rm=TRUE),
#             distance_tot_q75 = quantile(distance_tot, 0.75, na.rm=TRUE),
#             distance_tot_q25 = quantile(distance_tot, 0.25, na.rm=TRUE),
#             distance_tot_q975 = quantile(distance_tot, 0.975, na.rm=TRUE),
#             distance_tot_q025 = quantile(distance_tot, 0.025, na.rm=TRUE),
#             distance_tot_q90 = quantile(distance_tot, 0.9, na.rm=TRUE),
#             distance_tot_q10 = quantile(distance_tot, 0.1, na.rm=TRUE),
#             prediction_tot_mean= mean(prediction_tot),
#             prediction_tot_q50 = quantile(prediction_tot, 0.5, na.rm=TRUE),
#             prediction_tot_q25 = quantile(prediction_tot, 0.75, na.rm=TRUE),
#             prediction_tot_q75 = quantile(prediction_tot, 0.25, na.rm=TRUE),
#             distance_U_q50 = quantile(distance_U, 0.5, na.rm=TRUE),
#             distance_U_q25 = quantile(distance_U, 0.75, na.rm=TRUE),
#             distance_U_q75 = quantile(distance_U, 0.25, na.rm=TRUE),
#             prediction_U_mean= mean(distance_U),
#             prediction_U_q50 = quantile(distance_U, 0.5, na.rm=TRUE),
#             prediction_U_q25 = quantile(distance_U, 0.75, na.rm=TRUE),
#             prediction_U_q75 = quantile(distance_U, 0.25, na.rm=TRUE),
#             distance_V_q50 = quantile(distance_V, 0.5, na.rm=TRUE),
#             distance_V_q25 = quantile(distance_V, 0.75, na.rm=TRUE),
#             distance_V_q75 = quantile(distance_V, 0.25, na.rm=TRUE),
#             prediction_V_mean= mean(distance_V),
#             prediction_V_q50 = quantile(distance_V, 0.5, na.rm=TRUE),
#             prediction_V_q25 = quantile(distance_V, 0.75, na.rm=TRUE),
#             prediction_V_q75 = quantile(distance_V, 0.25, na.rm=TRUE),
#             TPR_q50 = quantile(TPR, 0.5, na.rm=TRUE),
#             TPR_q25 = quantile(TPR, 0.75, na.rm=TRUE),
#             TPR_q75 = quantile(TPR, 0.25, na.rm=TRUE),
#             FPR_mean = mean(FPR, na.rm=TRUE),
#             FPR_q50 = quantile(FPR, 0.5, na.rm=TRUE),
#             FPR_q25 = quantile(FPR, 0.75, na.rm=TRUE),
#             FPR_q75 = quantile(FPR, 0.25, na.rm=TRUE),
#             FNR_mean = mean(FNR, na.rm=TRUE),
#             FNR_q50 = quantile(FNR, 0.5, na.rm=TRUE),
#             FNR_q25 = quantile(FNR, 0.75, na.rm=TRUE),
#             FNR_q75 = quantile(FNR, 0.25, na.rm=TRUE),
#             time_med = quantile(time, 0.5, na.rm=TRUE),
#             time_mean = mean(time),
#             counts = n()
#             
#   ) %>%
#   ungroup() 
# 
# unique(summ$method)
# legend_order <- c("Oracle",  "FIT_SAR_CV",
#                   "FIT_SAR_BIC", "Witten_Perm", "Witten.CV",
#                   "SCCA_Parkhomenko", "Waaijenborg-CV", "Waaijenborg-Author",
#                   #"RRR-0.5" ,"RRR-7.5","RRR-10","RRR-12.5",  "RRR-20",
#                   #"RRR-opt",    
#                   "RRR-opt" ,
#                   "RRR-ADMM-opt",
#                   "CVX-opt-graph" )
# my_colors <- c( "black", "red", "indianred4",
#                 "orange", "yellow", "chartreuse2",
#                 "burlywood2", "burlywood4",
#                 # "lightblue", "lightblue3","cyan", "dodgerblue", "dodgerblue4",
#                 "navyblue", 
#                 "cyan", 
#                 "dodgerblue")
# 
# labels_n <-    c("Oracle",  "SAR CV (Wilms et al)",
#                  "SAR BIC (Wilms et al)",
#                  "Sparse CCA, permuted\n(Witten et al)",
#                  "Sparse CCA with CV\n(Witten et al)",
#                  "SCCA (Parkhomenko et al)", "Sparse CCA with CV\n(Waaijenborg et al)",
#                  "Sparse CCA(Waaijenborg et al)",
#                  # "RRR-0.5" ,"RRR-7.5","RRR-10","RRR-12.5",  "RRR-20",
#                  "RRR-CCA (this paper, using rrr solver)",
#                  "RRR-CCA (this paper)",
#                  "graph-RRR-CCA (this paper)")
# theme_set(theme_bw(base_size = 18))
# 
# 
# unique(summ$method)
# colnames(results)
# colnames(summ)
# unique(summ$nnzeros)
# unique(summ$r)
# unique(summ$r_pca)
# unique(summ$p1)
# unique(summ$p2)
# unique(summ$counts)
# 
# summ$theta_strength <- factor(summ$theta_strength, levels = c("high", "medium", "low"))
# 
# 
# ggplot(summ %>% filter( r_pca == 5, r==2,
#                         #nnzeros==10, 
#                         n==500,
#                         method %in% legend_order
# ),
# aes(x=p1, 
#     y = distance_tot_mean, 
#     colour =method)) +
#   geom_point(size=2)+
#   geom_line(linewidth=0.8)+
#   #geom_errorbar(aes(ymin=distance_tot_q975, ymax=distance_tot_q025,
#   #                  colour =method), width=0.05, alpha=0.5)+
#   scale_color_manual(values = my_colors, breaks = legend_order,
#                      labels = labels_n) +
#   facet_grid(theta_strength~ p2, scales = "free",labeller = as_labeller(c(`10` = "q = 10",
#                                                                           `30` = "q = 30",
#                                                                           `50` = "q = 50",
#                                                                           `80` = "q = 80",
#                                                                           `100` = "n = 100",
#                                                                           `200` = "n = 200",
#                                                                           `300` = "n = 300",
#                                                                           `500` = "n = 500",
#                                                                           `high` = "High",
#                                                                           `1000` = "n = 1,000",
#                                                                           `2000` = "n = 2,000",
#                                                                           `10000` = "n = 10,000",
#                                                                           `medium` = "Medium",
#                                                                           `low` = "Low"
#                                                                           
#   ))) +
#   xlab("p") + 
#   ylab(expression("Subspace Distance")) +
#   labs(colour="Method") + 
#   scale_y_log10()+
#   scale_x_log10()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))



