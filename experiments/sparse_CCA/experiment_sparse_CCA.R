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
criterion <- args[4]
set.seed(seed)
it = seed
THRES = 1e-3
r = 2 

for (n in c(N)){
for (psize in c( 0.25 *n, 0.5 *n , 0.75 *n, n, 1.25 *n,  1.5 *n , 2 *n, 2.5*n, 3*n, 4*n, 5 *n)){   
 for (sparsity in c(0.05, 0.1, 0.2, 0.3, 0.5)){
for (sparsity in c(0.1)){
    p1=as.integer(psize); p2=as.integer(psize)
    nnz = ceil(sparsity * p1)
    print(c(n, p1, p2))
    p = p1+ p2
    example <- generate_example(n=n, p1=p1, p2=p2,   
                                nnzeros = min(nnz, min(p1,p2)-1),
                                theta = diag( c(0.9,  0.8)),
                                a = 0.3, r=5)
    silly_benchmark = subdistance(matrix(0, p1,2), example$u)
    print("here")
    max1 = 500 * sqrt(log(p1)/n)
    min1 = 0.001 * sqrt(log(p1)/n) 
    param1 = exp(seq(log(min1), log(max1), length.out=20))
       # c(5, 10, 20, 30, 50, 80, 100, 200, 300,  500, 700, 1000)
    maxk = 0.25 * p
    mink = 0.01 * p 
    param2 = ceiling(seq(max(ceiling(mink),5), ceiling(maxk), length.out = 10))
    fantope_solution = NULL

      for (adaptive in c(TRUE, FALSE)){
        for (initialize in c("Fantope", "Selection")){
          result <- tryCatch({
          name_method = ifelse(adaptive, "adaptive_lasso", "lasso")
          name_method = paste0(name_method, "_",  initialize)
        
          res = pipeline_adaptive_lasso(example$Data, example$Mask, example$sigma0hat, 
                                    r=r, 
                                    nu=1, example$Sigmax, 
                                    example$Sigmay, maxiter=100, lambdax=NULL,
                                    adaptive=adaptive, kfolds=5,  param1=param1,
                                    create_folds=FALSE, normalize=FALSE,
                                    init=initialize, alpha = 0.75,
                                    criterion=criterion, 
                                    fantope_solution = fantope_solution )
          if (initialize == 'Fantope'){
            fantope_solution = rbind(res$initu, res$initv)
          }
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
                         "nb_real_discoveries" = sum(apply(Uhat^2, 1, sum)>THRES),
                         "param1" = res$lambdax,
                         "param2" = res$lambday,
                         "distance" = subdistance(Uhat, example$a),
                          "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum), tol=THRES))
         
      if (length(results)==0){
          results=temp
        }else{
           results <- rbind(results, temp )
           if(name_method == "lasso_Selection"){
             Uhat = rbind(res$initu, res$initv)
             temp <- data.frame("method" = "Selection",
                                "exp" = it,
                                "n" = n,
                                "nnz" = nnz,
                                "p1" = p1,
                                "p2" = p2,
                                "sparsity" = sparsity,
                                "zero_benchmark" = silly_benchmark,
                                "nb_discoveries" = sum(apply(Uhat^2, 1, sum)>0),
                                "nb_real_discoveries" = sum(apply(Uhat^2, 1, sum)>THRES),
                                "param1" = res_tg$lambda,
                                "param2" = res_tg$k,
                                "distance" = subdistance(Uhat, example$a),
                                "principal_angles" = principal_angles(Uhat, example$a),
                                "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                                "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                                "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                                "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum), tol=THRES))
             results <- rbind(results, temp)
           }

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
      
      
       test = pipeline_alternating_lasso(example$Data, example$Mask, example$sigma0hat, 
                                  r=2, 
                                  nu=1, example$Sigmax, 
                                  example$Sigmay, maxiter=100, lambdax=NULL,
                                  adaptive=adaptive, kfolds=5,  param1=param1,
                                  create_folds=FALSE, normalize=FALSE,
                                  init=initialize, alpha = 0.75,
                                  criterion=criterion, 
                                  fantope_solution = fantope_solution )
       
       Uhat = rbind(test$ufinal$Uhat, test$vfinal$Uhat)
       temp <- data.frame("method" = name_method,
                          "exp" = it,
                          "n" = n,
                          "nnz" = nnz,
                          "p1" = p1,
                          "p2" = p2,
                          "sparsity" = sparsity,
                          "zero_benchmark" = silly_benchmark,
                          "nb_discoveries" = sum(apply(Uhat^2, 1, sum)>0),
                          "nb_real_discoveries" = sum(apply(Uhat^2, 1, sum)>THRES),
                          "param1" = res$lambdax,
                          "param2" = res$lambday,
                          "distance" = subdistance(Uhat, example$a),
                          "principal_angles" = principal_angles(Uhat, example$a),
                          "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                          "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                          "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                          "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum), tol=THRES))
       
       
     result <- tryCatch({
      normalize=FALSE
      name_method = paste0("TG", ifelse(normalize, "_normalized", ""))
      res_tg <- pipeline_thresholded_gradient(example$Data, example$Mask, example$sigma0hat, 
                                              r=2, nu=1,Sigmax=example$Sigmax, 
                                              Sigmay=example$Sigmay, maxiter.init=100, lambda=NULL,k=NULL,
                                              kfolds=5, maxiter=2000, convergence=1e-3, eta=1e-3,
                                              param1=param1,
                                              param2=param2, normalize=normalize,
                                              criterion=criterion,
                                              fantope_solution=fantope_solution)
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
                         "nb_real_discoveries" = sum(apply(Uhat^2, 1, sum)>THRES),
                         "param1" = res_tg$lambda,
                         "param2" = res_tg$k,
                         "distance" = subdistance(Uhat, example$a),
                         "principal_angles" = principal_angles(Uhat, example$a),
                         "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum), tol=THRES))
      
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
                         "nb_real_discoveries" = sum(apply(Uhat^2, 1, sum)>THRES),
                         "param1" = res_tg$lambda,
                         "param2" = res_tg$k,
                         "distance" = subdistance(Uhat, example$a),
                         "principal_angles" = principal_angles(Uhat, example$a),
                         "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum), tol=THRES))
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
                         "nb_real_discoveries" = sum(apply(Uhat^2, 1, sum)>THRES),
                         "param1" = res_tg$lambda,
                         "param2" = res_tg$k,
                         "distance" = subdistance(Uhat, example$a),
                         "principal_angles" = principal_angles(Uhat, example$a),
                         "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                         "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum), tol=THRES))
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
                            "nb_real_discoveries" = NA,
                             "param1" = NA,
                             "param2" = NA,
                             "distance" = NA,
                            "principal_angles" = NA,
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
     
     
  #### Oracle
  set_u =  which(apply(example$u,1, norm)>0)
  set_v = nrow(example$u) + which(apply(example$v,1, norm)>0)
  t=CCA::cc(as.matrix(example$Data[,set_u]), as.matrix(example$Data[, set_v]))
  Uhat = matrix(0, p, r)
  Uhat[set_u, ] =  t$xcoef[,1:r]
  Uhat[set_v, ] =  t$ycoef[,1:r]
  temp <- data.frame("method" = "Oracle",
                     "exp" = it,
                     "n" = n,
                     "nnz" = nnz,
                     "p1" = p1,
                     "p2" = p2,
                     "sparsity" = sparsity,
                     "zero_benchmark" = silly_benchmark,
                     "nb_discoveries" = sum(apply(Uhat^2, 1, sum)>0),
                     "nb_real_discoveries" = sum(apply(Uhat^2, 1, sum)>THRES),
                     "param1" = NA,
                     "param2" = NA,
                     "distance" = subdistance(Uhat, example$a),
                     "principal_angles" = principal_angles(Uhat, example$a),
                     "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                     "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                     "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                     "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum), tol=THRES))
  results <- rbind(results, temp )

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
                           "nb_real_discoveries" = sum(apply(Uhat^2, 1, sum)>THRES),
                           "param1" = NA,
                           "param2" = NA,
                           "distance" = subdistance(Uhat, example$a),
                           "principal_angles" = principal_angles(Uhat, example$a),
                           "TPR" =TPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                           "TNR" = TNR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                           "FPR" = FPR(apply(Uhat^2, 1, sum), apply(example$a^2, 1, sum), tol=THRES),
                           "FNR" = FPR(apply(example$a^2, 1, sum),apply(Uhat^2, 1, sum), tol=THRES))
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
      write_excel_csv(results, paste0("experiments/sparse_CCA/results/extended_results_exp_sparse_cca_", name_exp, "_", criterion, ".csv"))
    }
  }
}


