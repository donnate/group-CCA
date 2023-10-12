wd = getwd()
#wd = "~/Documents/group-CCA/"
print(wd)
setwd(wd)
source("experiments/sparse_CCA/experiment_functions.R")
source("experiments/sparse_CCA/generate_examples.R")
source("src/ssvd/my_ssvd.R")
results <- c()

args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
print(seed)
name_exp <- args[2]
n <- as.integer(as.numeric(args[3]))
r <- ceiling(as.numeric(args[4]))
r_pca <- ceiling(as.numeric(args[5]))
criterion <- args[6]
normalize_diagonal <- as.logical(as.numeric(args[7]))
set.seed(seed)
it = seed
THRES = 1e-3


for (psize in c(0.5 *n , 0.75 *n, n, 
                1.25 *n,  1.5 *n , 2 *n, 2.5*n, 3*n, 4*n, 5 *n)){   
    for (sparsity in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4)){
      nnz = ceil(sparsity * psize)
      if (nnz > 2){
         for (overlapping_amount in seq(0,1 , 0.1)){
          p1 <- as.integer(psize); 
          p2 <- as.integer(psize)
          print(c(n, p1, p2, nnz))
          p <-  p1 + p2
          example <- generate_example_none_trivial_pca(n, p1, p2, 
                                              r_pca = 3, 
                                              nnzeros = nnz,
                                              theta = diag(seq(from = 0.9, to=0.6, 
                                                               length.out = r)),
                                              lambda_pca = 1,
                                              r=r, 
                                              overlapping_amount = overlapping_amount,
                                              normalize_diagonal=normalize_diagonal)
          print("here")
          max1 = 500 * sqrt(log(p1)/n)
          min1 = 0.001 * sqrt(log(p1)/n) 
          param1 = exp(seq(log(min1), log(max1), length.out=20))
          maxk = 0.25 * p
          mink = 0.01 * p 
          param2 = ceiling(seq(max(ceiling(mink),5), ceiling(maxk), length.out = 10))
          fantope_solution = NULL
          
          ssvd_results <- tryCatch({
            test1<-my.ssvd(example$S[1:p1, (p1+1):p],
                           Sigma_u = example$S[1:p1, 1:p1],
                           Sigma_v=example$S[(p1+1):p, (p1+1):p],
                           r=r, method = "theory")
            Uhat = rbind(test1$u, test1$v)
            temp <- evaluate_results(Uhat= test1$u, 
                                     Vhat = test1$v, 
                                     example = example, 
                                     name_method="SSVD-theory", 
                                     overlapping_amount=overlapping_amount,
                                     lambdax= NA,
                                     lambday = NA, 
                                     thres = THRES,
                                     it=it,
                                     normalize_diagonal=normalize_diagonal,
                                     criterion=criterion)
            if (length(results)==0){
              results=temp
            }else{
              results <- rbind(results, temp )
              
            }
          }, error = function(e) {
            # Print the error message
            cat("Error occurred in method sparse SSVD (theory)", ":", conditionMessage(e), "\n")
            # Skip to the next iteration
          })
          
          ssvd_results <- tryCatch({
            test1<-my.ssvd(example$S[1:p1, (p1+1):p],
                           Sigma_u = example$S[1:p1, 1:p1],
                           Sigma_v=example$S[(p1+1):p, (p1+1):p],
                           r=r, 
                           method = "method")
            Uhat = rbind(test1$u, test1$v)
            temp <- evaluate_results(Uhat= test1$u, 
                                     Vhat = test1$v, 
                                     example = example, 
                                     name_method="SSVD-method", 
                                     overlapping_amount=overlapping_amount,
                                     lambdax= NA,
                                     lambday = NA, 
                                     thres = THRES,
                                     it=it,
                                     normalize_diagonal=normalize_diagonal,
                                     criterion=criterion)
            if (length(results)==0){
              results=temp
            }else{
              results <- rbind(results, temp )
              
            }
          }, error = function(e) {
            # Print the error message
            cat("Error occurred in method sparse SSVD (theory)", ":", conditionMessage(e), "\n")
            # Skip to the next iteration
          })
          
          
          for (adaptive in c(TRUE, FALSE)){
              for (initialize in c("Fantope", "Selection")){
                result <- tryCatch({
                                    name_method = ifelse(adaptive, "adaptive_lasso", 
                                                        "lasso")
                                    name_method = paste0(name_method, "_",  
                                                         initialize)
              
                                    res = pipeline_adaptive_lasso(example$Data, 
                                                                  example$Mask, 
                                                                  example$sigma0hat,
                                                                  r=r, nu=1, 
                                                                  example$Sigmax, 
                                                                  example$Sigmay, 
                                                                  maxiter=100, 
                                                                  lambdax=NULL, 
                                                                  adaptive=adaptive, 
                                                                  kfolds=5,  
                                                                  param1=param1,
                                                                  create_folds=FALSE, 
                                                                  normalize=FALSE,
                                                                  init=initialize, 
                                                                  alpha = 0.75,
                                                                  criterion=criterion, 
                                                                  fantope_solution = fantope_solution )
                                    if (initialize == 'Fantope'){
                                      fantope_solution = rbind(res$initu, res$initv)
                                    }
                                    temp <- evaluate_results(Uhat= res$Uhat, 
                                                             Vhat = res$Vhat, 
                                                             example = example, 
                                                             name_method=name_method, 
                                                             overlapping_amount=overlapping_amount,
                                                             lambdax= res$lambdax,
                                                             lambday = res$lambday, 
                                                             thres = THRES,
                                                             it=it,
                                                             normalize_diagonal=normalize_diagonal,
                                                             criterion=criterion)
                                    if (length(results)==0){
                                        results=temp
                                      }else{
                                        results <- rbind(results, temp )
                                        if(name_method == "lasso_Selection"){
                                          temp <- evaluate_results(Uhat= res$initu, 
                                                                   Vhat = res$initv, 
                                                                   example = example, 
                                                                   name_method=paste0(name_method, "_init"), 
                                                                   overlapping_amount=overlapping_amount,
                                                                   lambdax= res$lambdax,
                                                                   lambday = res$lambday, 
                                                                   thres = THRES,
                                                                   it=it,
                                                                   normalize_diagonal=normalize_diagonal,
                                                                   criterion=criterion)
                                          results <- rbind(results, temp)
                                        }
                    
                                      }
                                    }, error = function(e) {
                                                             # Print the error message
                                      cat("Error occurred in method", name_method, ":", conditionMessage(e), "\n")
                                      # Skip to the next iteration
                               })
            } 
          }

          print("Done with the adaptive stuff")
          
          result <- tryCatch({
            normalize=FALSE
            name_method = paste0("TG", ifelse(normalize, "_normalized", ""))
            res_tg <- pipeline_thresholded_gradient(example$Data, example$Mask, 
                                                    example$sigma0hat, 
                                                    r=r, nu=1,
                                                    Sigmax=example$Sigmax, 
                                                    Sigmay=example$Sigmay, 
                                                    maxiter.init=100, 
                                                    lambda=NULL,k=NULL,
                                                    kfolds=5, maxiter=2000, 
                                                    convergence=1e-3, eta=1e-3,
                                                    param1=param1,
                                                    param2=param2, normalize=normalize,
                                                    criterion=criterion,
                                                    fantope_solution=fantope_solution)
            Uhat = rbind(res_tg$ufinal, res_tg$vfinal)
            temp <- evaluate_results(Uhat= res_tg$ufinal, 
                                     Vhat = res_tg$vfinal, 
                                     example = example, 
                                     name_method=name_method, 
                                     overlapping_amount=overlapping_amount,
                                     lambdax=  res_tg$lambda,
                                     lambday =   res_tg$k, 
                                     thres = THRES,
                                     it=it,
                                     normalize_diagonal=normalize_diagonal,
                                     criterion=criterion)
            
          results <- rbind(results, temp)

          #Uhat = rbind(res_tg$initu, res_tg$initv)
          temp <- evaluate_results(Uhat= res_tg$initu, 
                                   Vhat = res_tg$initv, 
                                   example = example, 
                                   name_method="Fantope", 
                                   overlapping_amount=overlapping_amount,
                                   lambdax=  res_tg$lambda,
                                   lambday =   res_tg$k, 
                                   thres = THRES,
                                   it=it,
                                   normalize_diagonal=normalize_diagonal,
                                   criterion=criterion)
          results <- rbind(results, temp)
          selected_rows = which(apply(res_tg$initu^2, 1, sum)>0)
          selected_rows.v = which(apply(res_tg$initv^2, 1, sum)>0)
          print("Selected rows.v")
          print(selected_rows.v)

          print("Selected rows")
          print(selected_rows)
          
          
          
          if (min(length(selected_rows), length(selected_rows.v))<n & min(length(selected_rows), length(selected_rows.v))>0){
                print("here CCA")
                t = cancor(as.matrix(example$Data[,selected_rows]), as.matrix(example$Data[, (selected_rows.v + p1)]))
                Uhat = matrix(0, p, r)
                Uhat[selected_rows,]  = t$xcoef[,1:r]
                Uhat[selected_rows.v + p1,]  = t$ycoef[,1:r]
                print(paste0("Done with Uhat: dimension:"))
                
                temp <- evaluate_results(Uhat= Uhat[1:p1,],
                                         Vhat = Uhat[(p1+1):p,], 
                                         example = example, 
                                         name_method="thresholded-lasso", 
                                         overlapping_amount=overlapping_amount,
                                         lambdax= res_tg$lambda,
                                         lambday =  res_tg$k, 
                                         thres = THRES,
                                         it=it,
                                         normalize_diagonal=normalize_diagonal,
                                         criterion=criterion)
              
            print(dim(temp))
            print(dim(results))
            print("done")
            }else{
                print("Couldnt threshold")
          
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
        temp <- evaluate_results(Uhat= Uhat[1:p1,],
                                 Vhat = Uhat[(p1+1):p,], 
                                 example = example, 
                                 name_method="Oracle", 
                                 overlapping_amount=overlapping_amount,
                                 lambdax= NA,
                                 lambday =  NA, 
                                 thres = THRES,
                                 it=it,
                                 normalize_diagonal=normalize_diagonal,
                                 criterion=criterion)
        results <- rbind(results, temp )

        for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                            "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                            "SCCA_Parkhomenko", "Canonical Ridge-Author")){
            

              result <- tryCatch({
                            test1<-additional_checks(example$Data[,1:p1],
                                                    example$Data[,(p1+1):(p2+p1)], S=NULL, 
                                                    rank=r, kfolds=5, 
                                                    method.type = method)
                            Uhat = rbind(test1$u, test1$v)
                            temp <- evaluate_results(Uhat= test1$u,
                                                     Vhat = test1$v, 
                                                     example = example, 
                                                     name_method=method, 
                                                     overlapping_amount=overlapping_amount,
                                                     lambdax= NA,
                                                     lambday =  NA, 
                                                     thres = THRES,
                                                     it=it,
                                                     normalize_diagonal=normalize_diagonal,
                                                     criterion=criterion)
                            
                            if (length(results)==0){
                              results=temp
                            }else{
                              results <- rbind(results, temp )
              
                            }
                            }, error = function(e) {
                              # Print the error message
                              cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
                              # Skip to the next iteration
                    })
      }
      write_excel_csv(results, paste0("experiments/sparse_CCA/results/extended_results_exp_sparse_cca_", name_exp, "_", criterion, ".csv"))
      }
    }
  }
}



