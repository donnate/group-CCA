library(ggplot2)
library(tidyverse)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
#setwd("~/Documents/group-CCA/")
source("elena/missing/helper.R")
source("elena/missing/evaluation.R")
source("elena/missing/original_CCA_impute.R")
source("elena/iterative_cca.R")

#Simulation for missing values in both X and Y

seeds = 1:100
ps = c(5, 10, 20, 50, 80, 100, 200)

props = seq(0, 0.5, 0.05)
args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
print(seed)
name_exp <- args[2]
set.seed(seed)
strength_theta <- args[3]



overlaps = c(1)
result <- c()
for(seed_n in seeds){
  set.seed(seed + seed_n)
  for (n in c(100, 200, 500, 1000)){
    for(p in ps){
      if (p < n){
        for(s in c(1)){
          rs = ifelse( p < 6, c(2, 3, 5), c(2, 3, 5, 10))
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
            
            q = p
            for (r_pca in c(0, p)){
              for (overlapping_amount in overlaps){
                for(prop_missing in props){
                  cat("\n\ns:", s, "p:", p, "props:", prop_missing, "\n")
                  cat("seed:")
                  cat(seed, " ")
                  gen = generate_example_non_trivial_pca(n, p, q,
                                                    r_pca = r_pca,
                                                    nnzeros = p,
                                                    noise = s,
                                                    theta = thetas,
                                                    lambda_pca = 1,
                                                    r = r,
                                                    overlapping_amount = overlapping_amount,
                                                    normalize_diagonal = FALSE,
                                                    prop_missing = prop_missing) 
                  X = gen$X
                  Y = gen$Y
                  Xna = gen$Xna
                  Yna = gen$Yna
                  Sigma0_sqrt_inv = sqrtm(gen$Sigma)$Binv
                  Sigma_hat_sqrt_inv = sqrtm(gen$S)$Binv
                  Ximp = data.frame(gen$Xna) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE)))
                  Yimp = data.frame(Yna) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE)))
                  
                  Ximp2 = data.frame(gen$Xna) %>% mutate_all(~replace_na(., median(., na.rm = TRUE)))
                  Yimp2 = data.frame(gen$Yna) %>% mutate_all(~replace_na(., median(., na.rm = TRUE)))
                  
                  rrr = CCAimpute(gen$Xna, gen$Yna, k = r, eps = 1e-4, verbose = F)
                  result = rbind(result, data.frame(evaluate(gen$newX, gen$newY, rrr$U, rrr$V, gen$u, 
                                                            gen$v,
                                                            Sigma_hat_sqrt_inv = Sigma_hat_sqrt_inv, 
                                                            Sigma0_sqrt_inv = Sigma0_sqrt_inv), 
                                                    "noise" = s, "method" = "RRR",  "prop_missing" = prop_missing,
                                                    "overlapping_amount" = overlapping_amount,
                                                    "n" = n,
                                                    "strength_theta" = strength_theta,
                                                    "r_pca" = r_pca,
                                                    "exp" = seed + seed_n))
                  
                  
                  alt = alternating_cca(gen$Xna, gen$Yna, r)
                  result = rbind(result, data.frame(evaluate(gen$newX, gen$newY, alt$U, alt$V, gen$u, 
                                                            gen$v,
                                                            Sigma_hat_sqrt_inv = Sigma_hat_sqrt_inv, 
                                                            Sigma0_sqrt_inv = Sigma0_sqrt_inv), 
                                                    "noise" = s, "method" = "Alt",  "prop_missing" = prop_missing,
                                                    "overlapping_amount" = overlapping_amount,
                                                    "n" = n,
                                                     "strength_theta" = strength_theta,
                                                    "r_pca" = r_pca,
                                                    "exp" = seed + seed_n))
                  cca = CCA::cc(Ximp, Yimp)
                  result = rbind(result, data.frame(evaluate(gen$newX, gen$newY, cca$xcoef[, 1:r], 
                                                            cca$ycoef[, 1:r], 
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt_inv = Sigma_hat_sqrt_inv, 
                                                            Sigma0_sqrt_inv = Sigma0_sqrt_inv),
                                                    "noise" = s,  method = "CCA-mean",  
                                                    "prop_missing" = prop_missing, 
                                                    "overlapping_amount" = overlapping_amount,
                                                    "n" = n,
                                                    "r_pca" = r_pca,
                                                    "strength_theta" = strength_theta,
                                                    "exp" = seed + seed_n))
                  
                  cca = CCA::cc(Ximp2, Yimp2)
                  result = rbind(result, data.frame(evaluate(gen$newX, gen$newY, cca$xcoef[, 1:r], 
                                                            cca$ycoef[, 1:r], 
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt_inv = Sigma_hat_sqrt_inv, 
                                                            Sigma0_sqrt_inv = Sigma0_sqrt_inv),
                                                    "noise" = s,  method = "CCA-median",  
                                                    "prop_missing" = prop_missing, 
                                                    "overlapping_amount" = overlapping_amount,
                                                    "r_pca" = r_pca,
                                                    "strength_theta" = strength_theta,
                                                    "n" = n,
                                                    "exp" = seed + seed_n))
                  write_csv(result, paste0("elena/missing/simulation-RRR-results", name_exp, ".csv"))
            }
          }
        }
        }
      }
    }
  }
 }
}
