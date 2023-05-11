wd = getwd()
print(wd)
setwd(wd)
source("experiments/sparse_CCA/experiment_functions.R")
results <- c()


args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
print(seed)
name_exp <- args[2]
N <- as.numeric(args[3])
set.seed(seed)
it = seed
for (n in c(N)){
  for (psize in c(0.25 * n, 0.5 * n, 0.75 *n, n, 1.5 *n , 2 *n, 5 *n)){
    for (sparsity in c(0.1,  0.2, 0.3)){
    p1=as.integer(psize); p2=as.integer(psize)
    nnz = ceil(sparsity * p1)
    example <- generate_example(n=n, p1=p1, p2=p2,   nnzeros = nnz,
                                  theta = diag( c(0.9,  0.8)),
                                  a = 0, r=2)
      print("here")
      res = pipeline_adaptive_lasso(example$Data, example$Mask, example$sigma0hat, r=2, 
                                    nu=1, example$Sigmax, 
                                    example$Sigmay, maxiter=100, lambdax=NULL,
                                    adaptive=TRUE, kfolds=5,  param1=10^(seq(-5, 1, by = 0.5)),
                                    create_folds=FALSE)
      Uhat = rbind(res$Uhat, res$Vhat)
      temp <- data.frame("method" = "adaptive_lasso",
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
      print("there")
      
      res = pipeline_adaptive_lasso(example$Data, example$Mask, example$sigma0hat, r=2, 
                                    nu=1, example$Sigmax, 
                                    example$Sigmay, maxiter=100, lambdax=NULL,
                                    adaptive=TRUE, kfolds=5,  param1=10^(seq(-5, 1, by = 0.5)),
                                    create_folds=TRUE)
      Uhat = rbind(res$Uhat, res$Vhat)
      #plot(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum) )
      temp <- data.frame("method" = "adaptive_lasso_with_folds",
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
      
      results <-rbind(results, temp )
      
      res_tg <- pipeline_thresholded_gradient(example$Data, example$Mask, example$sigma0hat, 
                                              r=2, nu=1,Sigmax=example$Sigmax, 
                                              Sigmay=example$Sigmay, maxiter.init=100, lambda=NULL,k=NULL,
                                              kfolds=5,maxiter=2000, convergence=1e-3, eta=1e-3,
                                              param1=10^(seq(-5, 1, by = 0.5)),
                                              param2=c(10, 20, 30, 50, 100, 250, 500, 700, 1000))
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
      write_excel_csv(results, paste0("experiments/sparse_CCA/results/results_exp_sparse_cca_", name_exp, ".csv"))
      
    }
  }
}


