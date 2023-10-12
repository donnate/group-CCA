source("src/cd_solver.R")
source("src/penalties.R")
source("src/losses.R")

library(cvwrapr)

genglm <- function(x, y, D, family = c("gaussian", "binomial", "poisson"),
                   lambda1,
                   lambda2=0,
                   standardize = TRUE,
                   solver = c("ECOS", "CGD", "GLMNET" ),
                   intercept = TRUE,
                   eps=1e-3,
                   max.it = 200){


  n = nrow(x)
  p = ncol(x)
  mu  <- apply(X, 2, mean)
  if (standardize){
     X <-  as.matrix(data.frame(X) %>% mutate_all(~(scale(.) %>% as.vector)))
     if (intercept){
        X <- X + mu
     }
  }
  if (intercept){
    p <- p + 1
    X = cbind(matrix(rep(1, n), nrow=n), X)
    D = cbind(matrix(rep(0, p), nrow=n), D)
  }

  if (family == "gaussian"){
    if (solver == "ECOS"){
      beta1 = Variable(p)
      prob1 = Problem(Minimize(l2_loss(x, y, beta1) +
                                 gen_penalty(beta1,
                                             lambda1,
                                             lambda2,
                                             D)))
      CVXR_solution <- solve(prob1)
      hat_beta  = CVXR_solution$getValue(beta1)
      return(list("xcoef"=hat_beta, "lambda1"=lambda1,
                  "lambda2"=lambda2))
    }else{
      if (solver == "CGD"){
        hat_beta  = cgd_solver(x ,y, D, lambda1, lambda2,
                               eps = eps,  max.it)
        return(list("xcoef"=hat_beta,"lambda1"=lambda1,
                    "lambda2"=lambda2))
      }else{
            if (solver == "XCGD"){
        hat_beta  = xcgd_solver(x ,y, D, lambda1, lambda2,
                               eps = eps,  max.it)
        return(list("xcoef"=hat_beta,"lambda1"=lambda1,
                    "lambda2"=lambda2))
            }else{
              if (solver =="GLMNET"){
                X_tilde = x %*% pinv(D)
                glm.model = glmnet(X_tilde, y,
                              lambda = lambda1,
                              intercept = FALSE)
                hat_beta <- pinv(D) %*% as.matrix(glm.model$beta)
                return(list("xcoef"=hat_beta, "lambda1"=lambda1,
                            "lambda2"=lambda2))
                }else{
                  print("Not implemented yet")
                  return(NA)
                }
            }
      }
    }
  }
  if (family == "binomial"){
    if (solver == "ECOS"){
      print("Not implemented yet")
      return(NA)
    }
  }

  if (family == "poisson"){
    if (solver == "ECOS"){
      print("Not implemented yet")
      return(NA)
    }
  }


}
