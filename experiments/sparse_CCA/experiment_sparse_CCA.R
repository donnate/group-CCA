wd = getwd()
#wd = "~/Documents/group-CCA/"
print(wd)
setwd(wd)
source("experiments/sparse_CCA/experiment_functions.R")
results <- c()


args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
print(seed)
name_exp <- args[2]
N <- as.integer(as.numeric(args[3]))
set.seed(seed)
it = seed
for (n in c(N)){
<<<<<<< HEAD
  for (psize in c(n, 1.5 *n , 2 *n, 5 *n)){
    for (nnz in c(5, 10, 20, 30, 50)){
      psize= 1.5 * n
      nnz= ceiling(0.1*  n/(2 * 2) )
=======
#  for (psize in c(0.25 *n, 0.5 *n , 0.75 *n, n, 1.25 *n,  1.5 *n , 2 *n, 3*n,  5 *n)){
for (psize in c( 1.25 *n,  1.5 *n , 2 *n, 3*n,  5 *n)){   
 for (sparsity in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5)){
>>>>>>> 8026c28060dea1210e3b0bbe13b8e8da60aebc97
    p1=as.integer(psize); p2=as.integer(psize)
   
    nnz = ceil(sparsity * p1)
    print(c(n, p1, p2))
    p = p1+ p2
    example <- generate_example(n=n, p1=p1, p2=p2,   
                                nnzeros = min(nnz, min(p1,p2)-1),
                                theta = diag( c(0.9,  0.8)),
                                a = 0, r=2)
    silly_benchmark = subdistance(matrix(0, p,2), example$a)
      print("here")
      max1 = 100 * sqrt(log(p1)/n)
      min1 = 0.01 * sqrt(log(p1)/n) 
      param1 = seq(min1, max1, length.out=10)
       # c(5, 10, 20, 30, 50, 80, 100, 200, 300,  500, 700, 1000)
      maxk = 0.25 * p
      mink = 0.01 * p 
      param2 = ceiling(seq(max(ceiling(mink),5), ceiling(maxk), length.out = 8))

      for (adaptive in c(TRUE, FALSE)){
        for (create_folds in c(TRUE, FALSE)){
          result <- tryCatch({
        name_method = ifelse(adaptive, "adaptive_lasso", "lasso")
        name_method = paste0(name_method, ifelse(create_folds, "_with_folds", ""))

        res = pipeline_adaptive_lasso(example$Data, example$Mask, example$sigma0hat, 
                                    r=2, 
                                    nu=1, example$Sigmax, 
                                    example$Sigmay, maxiter=100, lambdax=NULL,
                                    adaptive=adaptive, kfolds=5,  param1=param1,
                                    create_folds=create_folds, normalize=FALSE)
         Uhat = rbind(res$Uhat, res$Vhat)

         temp <- data.frame("method" = name_method,
                         "exp" = it,
                         "n" = n,
                         "nnz" = nnz,
                         "p1" = p1,
                         "p2" = p2,
                         "sparsity" = sparsity,
                         "zero_benchmark" = silly_benchmark,
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
       }, error = function(e) {
    # Error handling code goes here
    
    # Print the error message
    cat("Error occurred in method", name_method, ":", conditionMessage(e), "\n")
    
    # Skip to the next iterationt
  })
        } 

      
     
      
      }

       print("there")
      
     result <- tryCatch({
      normalize=FALSE
      name_method = paste0("TG", ifelse(normalize, "_normalized", ""))
      res_tg <- pipeline_thresholded_gradient(example$Data, example$Mask, example$sigma0hat, 
                                              r=2, nu=1,Sigmax=example$Sigmax, 
                                              Sigmay=example$Sigmay, maxiter.init=100, lambda=NULL,k=NULL,
                                              kfolds=5, maxiter=2000, convergence=1e-3, eta=1e-3,
                                              param1=param1,
                                              param2=param2, normalize=normalize)
      Uhat = rbind(res_tg$ufinal, res_tg$vfinal)
      temp <- data.frame("method" = name_method,
                         "exp" = it,
                         "n" = n,
                         "nnz" = nnz,
                         "p1" = p1,
                         "p2" = p2,
                         "sparsity" = sparsity,
                          "zero_benchmark" = silly_benchmark,
                         "nb_discoveries" = sum(apply(Uhat^2, 1, sum)>0),
                         "param1" = res_tg$lambda,
                         "param2" = res_tg$k,
                         "distance" = subdistance(Uhat, example$a),
                         "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum)))
      
      results <- rbind(results, temp)


      Uhat = rbind(res_tg$initu, res_tg$initv)
      temp <- data.frame("method" = "Fantope",
                         "exp" = it,
                         "n" = n,
                         "nnz" = nnz,
                         "p1" = p1,
                         "p2" = p2,
                         "sparsity" = sparsity,
                          "zero_benchmark" = silly_benchmark,
                         "nb_discoveries" = sum(apply(Uhat^2, 1, sum)>0),
                         "param1" = res_tg$lambda,
                         "param2" = res_tg$k,
                         "distance" = subdistance(Uhat, example$a),
                         "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum)))
      results <- rbind(results, temp)
     selected_rows = which(apply(res_tg$initu^2, 1, sum)>0)
     selected_rows.v = which(apply(res_tg$initv^2, 1, sum)>0)
      print("Selected rows.v")
     print(selected_rows.v)

     print("Selected rows")
     print(selected_rows)
     if (min(length(selected_rows), length(selected_rows.v))<n & min(length(selected_rows), length(selected_rows.v))>0){
          r=2
          print("here CCA")
          t=cancor(as.matrix(example$Data[,selected_rows]), as.matrix(example$Data[, (selected_rows.v + p1)]))
          Uhat = matrix(0, p, r)
          Uhat[selected_rows,]  = t$xcoef[,1:r]
           Uhat[selected_rows.v + p1,]  = t$ycoef[,1:r]
           print(paste0("Done with Uhat: dimension:"))
          print(dim(Uhat))
            print(paste0("GT: dimension:"))
           print(dim(example$a))
           temp <- data.frame("method" = "thresholded-lasso",
                         "exp" = it,
                         "n" = n,
                         "nnz" = nnz,
                         "p1" = p1,
                         "p2" = p2,
                         "sparsity" = sparsity,
                          "zero_benchmark" = silly_benchmark,
                         "nb_discoveries" = sum(apply(Uhat^2, 1, sum)>0),
                         "param1" = res_tg$lambda,
                         "param2" = res_tg$k,
                         "distance" = subdistance(Uhat, example$a),
                         "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum)),
                         "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum)))
       print(dim(temp))
       print(dim(results))
       print("done")
       }else{
     print("Couldnt threshold")
         temp <- data.frame("method" = "thresholded-lasso",
                             "exp" = it,
                             "n" = n,
                             "nnz" = nnz,
                             "p1" = p1,
                             "p2" = p2,
                               "sparsity" = sparsity,
                             "zero_benchmark" = silly_benchmark,
                             "nb_discoveries" = NA,
                             "param1" = NA,
                             "param2" = NA,
                             "distance" = NA,
                             "TPR" = NA,
                             "TNR" = NA,
                             "FPR" = NA,
                             "FNR" = NA)
   }
      results <- rbind(results, temp )
     

       }, error = function(e) {
    # Error handling code goes here
    
    # Print the error message
    cat("Error occurred in method", name_method, ":", conditionMessage(e), "\n")
    
    # Skip to the next iteration
   
  })

      for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                       "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                       "SCCA_Parkhomenko", "Canonical Ridge-Author"
      )){
      

        result <- tryCatch({
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
                           "sparsity" = sparsity,
                            "zero_benchmark" = silly_benchmark,
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
        }, error = function(e) {
    # Error handling code goes here
    
    # Print the error message
    cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
    
    # Skip to the next iteration
      })
        
      }
      write_excel_csv(results, paste0("experiments/sparse_CCA/results/extended_results_exp_sparse_cca_", name_exp, ".csv"))
      
    }
  }
}


