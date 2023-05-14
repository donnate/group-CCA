wd = getwd()
#wd = "~/Documents/group-CCA/"
print(wd)
setwd(wd)
source("experiments/sparse_CCA/experiment_localized_functions.R")
source("experiments/sparse_CCA/experiment_functions.R")
results <- c()


args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
print(seed)
name_exp <- args[2]
N <- as.integer(as.numeric(args[3]))
type_graph <- args[4]
power <- as.numeric(args[5])
set.seed(seed)
it = seed
for (n in c(N)){
  for (psize in c(0.25 *n, 0.5*n, 0.75*n, n, 1.5 *n , 2 *n, 5 *n)){
    for (nnz in c(5, 10, 20)){
    p1=as.integer(psize); p2=as.integer(psize)
    #nnz = ceil(sparsity * n)
    print(c(n, p1, p2))
    p = p1+ p2
    example <- generate_localized_example(n=n, p1=p1,  
                                nnzeros = min(nnz, min(p1,p2)-1),
                                theta = diag( c(0.9,  0.8)),
                                r=2, type_graph=type_graph, power=power, probs=PROBA,
                                threshold_limit=0.8)
    
      
    final_nnz = sum(apply(example$a^2,1,sum)>1e-4)
      
    print("here")
    param1 = 10^(seq(-5, 3, by = 0.2))
    max1 = 20 * sqrt(log(p1)/n)
    min1 = 0.05 * sqrt(log(p1)/n) 
    param1 = param1[which(param1 < max1 & param1 > min1)]
     # c(5, 10, 20, 30, 50, 80, 100, 200, 300,  500, 700, 1000)
    maxk = 0.5 * p
    mink = 0.01 * p 
    param2 = ceiling(seq(max(ceiling(mink),5), ceiling(maxk), length.out = 8))
      
    transformed_data =  as.matrix(example$Data %*% example$daggerD)
    transformed_sigma0hat = as.matrix(t(example$daggerD) %*% example$sigma0hat %*% (example$daggerD))
    transformed_Sigmax =as.matrix( t(example$daggerDx) %*% example$Sigmax %*% (example$daggerDx))
    transformed_Sigmay =as.matrix( t(example$daggerDy) %*% example$Sigmay %*% (example$daggerDy))
    for (adaptive in c(TRUE, FALSE)){
        for (create_folds in c(TRUE, FALSE)){
          name_method = ifelse(adaptive, "adaptive_regularised_lasso", "regularised_lasso")
          name_method = paste0(name_method, ifelse(create_folds, "_with_folds", ""))
          
          
          res = pipeline_adaptive_lasso(transformed_data, example$Mask_dual, 
                                        transformed_sigma0hat, r=2, 
                                        nu=1, transformed_Sigmax, 
                                        transformed_Sigmay, maxiter=100, lambdax=NULL,
                                        adaptive=adaptive, kfolds=5,  param1=param1,
                                        create_folds=create_folds)
          Uhat = rbind(example$daggerDx %*% res$Uhat, example$daggerDy %*% res$Vhat)
          temp <- data.frame("method" = name_method,
                             "exp" = it,
                             "n" = n,
                             "nnz" = final_nnz,
                             "p1" = p1,
                             "p2" = p2,
                             "nb_discoveries" = sum(apply(Uhat^2, 1, sum)>1e-4),
                             "param1" = res$lambdax,
                             "param2" = res$lambday,
                             "distance" = subdistance(Uhat, example$a),
                             "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                             "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                             "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                             "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum)))
          if (length(results)==0){
            results=temp
          }else{
            results <- rbind(results, temp )
            
          }
          
        } 
      }
      
      print("Done witht hte ")
      for (adaptive in c(TRUE, FALSE)){
        for (create_folds in c(FALSE)){
          name_method = ifelse(adaptive, "adaptive_lasso", "lasso")
          name_method = paste0(name_method, ifelse(create_folds, "_with_folds", ""))
          
          res = pipeline_adaptive_lasso(example$Data, example$Mask, example$sigma0hat, r=2, 
                                        nu=1, example$Sigmax, 
                                        example$Sigmay, maxiter=100, lambdax=NULL,
                                        adaptive=adaptive, kfolds=5,  param1=param1,
                                        create_folds=create_folds)
          Uhat = rbind(res$Uhat, res$Vhat)
          temp <- data.frame("method" = name_method,
                             "exp" = it,
                             "n" = n,
                             "nnz" = nnz,
                             "p1" = p1,
                             "p2" = p2,
                             "nb_discoveries" = sum(apply(Uhat^2, 1, sum)>0),
                             "param1" = res$lambdax,
                             "param2" = res$lambday,
                             "distance" = subdistance(Uhat, example$a),
                             "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                             "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                             "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                             "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum)))
          if (length(results)==0){
            results=temp
          }else{
            results <- rbind(results, temp )
            
          }
          
        } 
      }

      print("there")
      res_tg <- pipeline_thresholded_gradient(example$Data, example$Mask, example$sigma0hat, 
                                              r=2, nu=1,Sigmax=example$Sigmax, 
                                              Sigmay=example$Sigmay, maxiter.init=100, lambda=NULL,k=NULL,
                                              kfolds=5, maxiter=2000, convergence=1e-3, eta=1e-3,
                                              param1=param1,
                                              param2=param2, normalize=normalize)
      Uhat = rbind(res_tg$ufinal, res_tg$vfinal)
      temp <- data.frame("method" = "TG",
                         "exp" = it,
                         "n" = n,
                         "nnz" = nnz,
                         "p1" = p1,
                         "p2" = p2,
                         "nb_discoveries" = sum(apply(Uhat^2, 1, sum)>0),
                         "param1" = res_tg$lambda,
                         "param2" = res_tg$k,
                         "distance" = subdistance(Uhat, example$a),
                         "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum)))
      
      results <- rbind(results, temp)
      
      for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                       "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                       "SCCA_Parkhomenko", "Canonical Ridge-Author"
      )){
        test1<-additional_checks(example$Data[,1:p1],
                                 example$Data[,(p1+1):(p2+p1)], S=NULL, 
                                 rank=2, kfolds=5, method.type = method)
        Uhat = rbind(test1$u, test1$v)
        temp <- data.frame("method" = method,
                           "exp" = it,
                           "n" = n,
                           "nnz" = nnz,
                           "p1" = p1,
                           "p2" = p2,
                           "nb_discoveries" = sum(apply(Uhat^2, 1, sum)>0),
                           "param1" = NA,
                           "param2" = NA,
                           "distance" = subdistance(Uhat, example$a),
                           "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                           "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                           "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                           "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum)))
        if (length(results)==0){
          results=temp
        }else{
           results <- rbind(results, temp )

        }
        
      }
      write_excel_csv(results, paste0("experiments/sparse_CCA/results/results_exp_localized_cca_", name_exp, ".csv"))
      
    }
  }
}


