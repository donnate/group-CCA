library(tidyverse)
library(igraph)
library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
setwd("~/Documents/group-CCA/")

source('experiments/sparse_CCA/experiment_functions.R')
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source("elena/missing/helper.R")
source("elena/missing/evaluation.R")
source("elena/missing/original_CCA_impute.R")
source("elena/gradient_descent.r")
source("elena/iterative_cca.R")
source("elena/reduced_rank_regression.R")
source("elena/graph_reduced_rank_regression.R")

get_edge_incidence <- function(g, weight = 1){
  n_nodes = vcount(g)
  d_max = max(degree(g))
  #d_max = 1
  edges = data.frame(as_edgelist(g)) %>%
    arrange(X1, X2)
  Gamma = matrix(0, nrow(edges), n_nodes)
  
  # Make beta_v into a matrix
  names_st = unique(c(edges$X1, edges$X2))
  for (e in 1:nrow(edges)){
    ind1 = which( edges$X1[e] == names_st)
    ind2 = which( edges$X2[e] == names_st)
    Gamma[e, ind1] = weight
    Gamma[e,ind2] = - weight
  }
  return(Gamma)
}

data_movies <- readxl::read_xlsx("~/Downloads/movies.xlsx")
data_movies <- data_movies %>% 
  filter(Year > 2022)
data_movies <- data_movies %>%
  mutate(`Total Gross` = str_remove_all(`Total Gross`, "\\$|,")) %>%
  mutate(`Total Gross` = as.numeric(`Total Gross`))

data_movies = data_movies[, -c(1, 3, 6, 7, 8, 9)]
states <- read_csv("~/Downloads/movie_states.csv")
movie_titles<- colnames(states)
state_names <- states$Region
states = t(as.matrix(states[, 2:22]))
states = data.frame(states)
colnames(states) = state_names
states["Title"] = movie_titles[2:22]


##### Analysis of politicians scores vs opinions on certain questions

data = merge(data_movies,
             states,
             #candidate_data,
             by = c("Title"))


ind_hawaii = which(state_names == "Hawaii")
ind_alaska = which(state_names == "Alaska")
X = data[, state_names[-c(ind_hawaii, ind_alaska)]]
numeric_columns = colnames(data_movies)[2:6]
Y = data[, numeric_columns]

### keep only small number
#numeric_columns <- numeric_columns[which(apply(Y, 2, function(x){mean(x>0)}) > 0.3)]
### party affiliation as a label for analysis
#Y = Y[, numeric_columns]
#### Find adacency map for
#minX = min(X[X>0])
#minY = min(Y[Y>0])


# Apply the logit function with smoothing to all columns
#X_transformed <- X %>%
#  mutate(across(everything(), log))
X_transformed <- X %>% 
  mutate(across(everything(), scale))
Y_transformed <- Y %>% 
  mutate(across(everything(), log)) %>%
  mutate(across(everything(), scale))


#### Download the adjacency matrix
A = read_csv("~/Downloads/state_adjacency.csv")
A = A[,2:ncol(A)]
### erase Hawaii
ind_hawaii = which(colnames(A) == "HI")
A = A[-c(ind_hawaii), -c(ind_hawaii)]

library(igraph)
# Create a graph from the adjacency matrix
g <- graph_from_adjacency_matrix(as.matrix(A), mode = "undirected")
plot(g,
     vertex.label = V(g)$name,  # Node labels
     vertex.size = 10,  # Adjust node size
     vertex.color = "skyblue",  # Node color
     edge.color = "grey",  # Edge color
     vertex.label.color = "black",  # Node label color
     edge.label.color = "red",  # Edge label color
     vertex.label.dist = 1.5,  # Distance of labels from nodes
     edge.label.cex = 0.8,  # Edge label size
     vertex.label.cex = 1.2  # Node label size
)

Gamma <- get_edge_incidence(g, weight = 1)
#### Now apply the algorithm
r = 3

p = ncol(X_transformed)
Y_transformed = scale(Y)

folds = createFolds(1:nrow(X_transformed),5)
correlation <- c()
lambda = 0.05
correlation <-c()
order = sample(1:length(folds), length(folds))
for (i in  1:length(folds)){
  index = order[ifelse(i < 18, i + 1, (i+1)%%18)]
  index2 =order[ifelse(i < 17, i + 2, (i+2)%%18)]
  print(c(i, index, index2))
  for (lambda in 10^seq(from=-3, 1, length.out=30)){
    final = CCA_graph_rrr(as.matrix(X_transformed)[-c(folds[[index]],
                                                      folds[[index2]]),], 
                          as.matrix(Y_transformed)[-c(folds[[index]],
                                                      folds[[index2]]),],  
                          Gamma, 
                          Sx=NULL, Sy=NULL, Sxy = NULL,
                          lambda = lambda, 
                          Kx=NULL, r=r,
                          rho=1, niter=2 * 1e4,
                          do.scale = FALSE, lambda_Kx=0,
                          thresh=1e-6,
                          LW_Sy = FALSE)
    
    correlation <- rbind(
      correlation,
      c("CCA_graph_rrr",
        lambda,
        i,
        diag(cov(as.matrix(X_transformed)[-c(folds[[index]],
                                             folds[[index2]]),] %*% final$U,
                 as.matrix(Y_transformed)[-c(folds[[index]],
                                             folds[[index2]]),] %*%  final$V)),
        apply(((as.matrix(X_transformed)[-c(folds[[index]],
                                            folds[[index2]]),] %*% final$U) -
                 (as.matrix(Y_transformed)[-c(folds[[index]],
                                              folds[[index2]]),] %*%  final$V))^2, 2, mean),
        diag(t(as.matrix(X_transformed)[folds[[index]][1],] %*% final$U) %*%
               (as.matrix(Y_transformed)[folds[[index]][1],] %*%  final$V)),
        ((as.matrix(X_transformed)[folds[[index]][1],] %*% final$U) -
           (as.matrix(Y_transformed)[folds[[index]][1],] %*%  final$V))^2,
        diag(t(as.matrix(X_transformed)[folds[[index2]][1],] %*% final$U) %*%
               (as.matrix(Y_transformed)[folds[[index2]][1],] %*%  final$V)),
        ((as.matrix(X_transformed)[folds[[index2]][1],] %*% final$U) -
           (as.matrix(Y_transformed)[folds[[index2]][1],] %*%  final$V))^2
      ))
    
  }
  for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                   "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                   "SCCA_Parkhomenko")){
    
    print(paste0("Starting ", method))
    tryCatch({
      test1<-additional_checks(as.matrix(X_transformed)[-c(
        folds[[index2]]),], as.matrix(Y_transformed)[-c(folds[[index2]]),],
        S=NULL, 
        rank=r, kfolds=5, 
        method.type = method,
        lambdax= 10^seq(-3,1, length.out = 30),
        lambday =  10^seq(-3,1, length.out = 30))
      
      testX = diag(t(as.matrix(X_transformed)[folds[[index]][1],] %*% test1$u) %*%  (as.matrix(X_transformed)[folds[[index]][1],] %*% test1$u))^(-0.5)
      testY = diag(t(as.matrix(Y_transformed)[folds[[index]][1],] %*% test1$v) %*%  (as.matrix(Y_transformed)[folds[[index]][1],] %*% test1$v))^(-0.5)
      valX = diag(t(as.matrix(X_transformed)[folds[[index2]][1],] %*% test1$u) %*%  (as.matrix(X_transformed)[folds[[index2]][1],] %*% test1$u))^(-0.5)
      valY = diag(t(as.matrix(Y_transformed)[folds[[index2]][1],] %*% test1$v) %*%  (as.matrix(Y_transformed)[folds[[index2]][1],] %*% test1$v))^(-0.5)
      
      correlation <- rbind(
        correlation,
        c(method,
          lambda,
          i,
          diag(cov(as.matrix(X_transformed)[-c(folds[[index]],
                                               folds[[index2]]),] %*% test1$u,
                   as.matrix(Y_transformed)[-c(folds[[index]],
                                               folds[[index2]]),] %*%  test1$v)),
          apply(((as.matrix(X_transformed)[-c(folds[[index]],
                                              folds[[index2]]),] %*% test1$u) -
                   (as.matrix(Y_transformed)[-c(folds[[index]],
                                                folds[[index2]]),] %*%  test1$v))^2, 2, mean),
          diag(t(as.matrix(X_transformed)[folds[[index]][1],] %*% test1$u) %*%
                 (as.matrix(Y_transformed)[folds[[index]][1],] %*%  test1$v)),
          ((as.matrix(X_transformed)[folds[[index]][1],] %*% test1$u) -
             (as.matrix(Y_transformed)[folds[[index]][1],] %*%  test1$v))^2,
          diag(t(as.matrix(X_transformed)[folds[[index2]][1],] %*% test1$u) %*%
                 (as.matrix(Y_transformed)[folds[[index2]][1],] %*%  test1$v)),
          ((as.matrix(X_transformed)[folds[[index2]][1],] %*% test1$u) -
             (as.matrix(Y_transformed)[folds[[index2]][1],] %*%  test1$v))^2
        ))
    }, error = function(e) {
      # Print the error message
      cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
      # Skip to the next iteration
    })
    
  }
}
STOP

final = CCA_graph_rrr.CV(as.matrix(X_transformed), as.matrix(Y_transformed),  
                         Gamma, 
                         kfolds = 2,
                         param_lambda = 10^seq(from=-3, 1, length.out=30), 
                         Kx=NULL, r=3,
                         rho=1, 
                         niter=2 * 1e4,
                         do.scale = FALSE, lambda_Kx=0,
                         thresh=1e-6,
                         LW_Sy = FALSE,
                         Gamma_dagger =  pinv(Gamma))

foldVector <- seq(from = 1, to = nrow(X_transformed), by = 3)
folds = split(sample(1:nrow(X_transformed), nrow(X_transformed)), foldVector)
correlation <-c()
r=2
order = 1:length(folds)
correlation <- c()
for (i in  1:length(folds)){
  index = order[ifelse(i < 6, i + 1, (i+1)%%6)]
  index2 =order[ifelse(i < 5, i + 2, (i+2)%%6)]
  print(c(i, index, index2))
  for (lambda in 10^seq(from=-3, 0, by =0.25)){
    final = CCA_graph_rrr(as.matrix(X_transformed)[-c(folds[[index]],
                                                      folds[[index2]]),], 
                          as.matrix(Y_transformed)[-c(folds[[index]],
                                                      folds[[index2]]),],  
                          Gamma, 
                          Sx=NULL, Sy=NULL, Sxy = NULL,
                          lambda = lambda, 
                          Kx=NULL, r=r,
                          rho=1, niter=2 * 1e4,
                          do.scale = FALSE, lambda_Kx=0,
                          thresh=1e-6,
                          LW_Sy = FALSE)
    
    correlation <- rbind(
      correlation,
      c("CCA_graph_rrr",
        lambda,
        i,
        diag(cov(as.matrix(X_transformed)[-c(folds[[index]],
                                             folds[[index2]]),] %*% final$U,
                 as.matrix(Y_transformed)[-c(folds[[index]],
                                             folds[[index2]]),] %*%  final$V)),
        apply(((as.matrix(X_transformed)[-c(folds[[index]],
                                            folds[[index2]]),] %*% final$U) -
                 (as.matrix(Y_transformed)[-c(folds[[index]],
                                              folds[[index2]]),] %*%  final$V))^2, 2, mean),
        diag(t(as.matrix(X_transformed)[folds[[index]],] %*% final$U) %*%
               (as.matrix(Y_transformed)[folds[[index]],] %*%  final$V)),
        diag(cor(as.matrix(X_transformed)[folds[[index]],] %*% final$U, (as.matrix(Y_transformed)[folds[[index]],] %*%  final$V))),
        apply(((as.matrix(X_transformed)[folds[[index]],] %*% final$U) -
                 (as.matrix(Y_transformed)[folds[[index]],] %*%  final$V))^2,2,mean),
        diag(t(as.matrix(X_transformed)[folds[[index2]],] %*% final$U) %*%
               (as.matrix(Y_transformed)[folds[[index2]],] %*%  final$V)),
        diag(cor(as.matrix(X_transformed)[folds[[index2]],] %*% final$U, (as.matrix(Y_transformed)[folds[[index2]],] %*%  final$V))),
        apply(((as.matrix(X_transformed)[folds[[index2]],] %*% final$U) -
                 (as.matrix(Y_transformed)[folds[[index2]],] %*%  final$V))^2,2,mean)
      ))
    
  }
  for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                   "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                   "SCCA_Parkhomenko")){
    
    print(paste0("Starting ", method))
    tryCatch({
      test1<-additional_checks(as.matrix(X_transformed)[-c(
        folds[[index2]]),], as.matrix(Y_transformed)[-c(
          folds[[index2]]),],
        S=NULL, 
        rank=r, kfolds=5, 
        method.type = method,
        lambdax= 10^seq(from=-3, 1, by =0.25),
        lambday = 10^seq(-4,-3, length.out = 5))
      
      testX = diag(t(as.matrix(X_transformed)[folds[[index]][1],] %*% test1$u) %*%  (as.matrix(X_transformed)[folds[[index]][1],] %*% test1$u))^(-0.5)
      testY = diag(t(as.matrix(Y_transformed)[folds[[index]][1],] %*% test1$v) %*%  (as.matrix(Y_transformed)[folds[[index]][1],] %*% test1$v))^(-0.5)
      valX = diag(t(as.matrix(X_transformed)[folds[[index2]][1],] %*% test1$u) %*%  (as.matrix(X_transformed)[folds[[index2]][1],] %*% test1$u))^(-0.5)
      valY = diag(t(as.matrix(Y_transformed)[folds[[index2]][1],] %*% test1$v) %*%  (as.matrix(Y_transformed)[folds[[index2]][1],] %*% test1$v))^(-0.5)
      
      correlation <- rbind(
        correlation,
        (c(method,
           lambda,
           i,
           diag(cov(as.matrix(X_transformed)[-c(folds[[index]],
                                                folds[[index2]]),] %*% test1$u,
                    as.matrix(Y_transformed)[-c(folds[[index]],
                                                folds[[index2]]),] %*%  test1$v)),
           apply(((as.matrix(X_transformed)[-c(folds[[index]],
                                               folds[[index2]]),] %*% test1$u) -
                    (as.matrix(Y_transformed)[-c(folds[[index]],
                                                 folds[[index2]]),] %*%  test1$v))^2, 2, mean),
           diag(t(as.matrix(X_transformed)[folds[[index]],] %*% test1$u) %*%
                  (as.matrix(Y_transformed)[folds[[index]],] %*%  test1$v)),
           diag(cor(as.matrix(X_transformed)[folds[[index]],] %*% test1$u, (as.matrix(Y_transformed)[folds[[index]],] %*%  test1$v))),
           apply(((as.matrix(X_transformed)[folds[[index]],] %*% test1$u) -
                    (as.matrix(Y_transformed)[folds[[index]],] %*%  test1$v))^2, 2, sum),
           diag(t(as.matrix(X_transformed)[folds[[index2]],] %*% test1$u) %*%
                  (as.matrix(Y_transformed)[folds[[index2]],] %*%  test1$v)),
           diag(cor(as.matrix(X_transformed)[folds[[index2]],] %*% test1$u, (as.matrix(Y_transformed)[folds[[index2]],] %*%  test1$v))),
           apply(((as.matrix(X_transformed)[folds[[index2]],] %*% test1$u) -
                    (as.matrix(Y_transformed)[folds[[index2]],] %*%  test1$v))^2, 2, sum)
        )))
    }, error = function(e) {
      # Print the error message
      cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
      # Skip to the next iteration
    })
    
  }
}
#folds = createFolds(1:nrow(X_transformed),6) ### 3 people per fold
### shuffle:

old_correlation = correlation_df

correlation_df = data.frame(correlation)
colnames(correlation_df) = c("method", "lambda", "fold",
                             "train_cov1",  "train_cov2",  
                             "train_mse1",  "train_mse2",  
                             "test_cov1",  "test_cov2",  
                             "test_cor1",  "test_cor2",  
                             "test_mse1", "test_mse2",
                             "val_cov1",  "val_cov2",  
                             "val_cor1",  "val_cor2",  
                             "val_mse1",  "val_mse2")
for (i in 3:19){
  correlation_df[,i] = as.numeric(correlation_df[,i])
}

correlation_df = rbind(correlation_df,
                       old_correlation %>% filter(method =="Witten.CV" | method == "Witten_Perm"))

summary_correlation = correlation_df %>% 
  mutate(test_mse = (test_mse1 + test_mse2 )/2,
         train_mse = train_mse1 + train_mse2,
         val_mse =  (val_mse1 + val_mse2 )/2,
         test_cov = (test_cov1 + test_cov2)/2,
         test_cor = (test_cor1 + test_cor2)/2,
         train_cov = (train_cov1 + train_cov2)/2,
         val_cov = (val_cov1 + val_cov2)/2,
         val_cor = (val_cor1 + val_cor2)/2) %>%
  group_by(method, lambda) %>% summarise_all(mean, na.rm=TRUE)%>%
  arrange((test_mse)) %>% ungroup()

ggplot(summary_correlation %>% 
         filter(method == "CCA_graph_rrr")%>% ungroup())+
  geom_line(aes(x = as.numeric(lambda), y=test_mse))+
  geom_point(aes(x = as.numeric(lambda), y=test_mse))+
  #scale_x_log10()
  
  summary_correlation$lambda = as.numeric(summary_correlation$lambda)
lambda_opt = summary_correlation$lambda[which.min(summary_correlation$train_mse[which(summary_correlation$method ==  "CCA_graph_rrr")])]

lambda_opt = 0.056234133
relevant_correlations = summary_correlation %>% 
  filter( (method == "CCA_graph_rrr" & lambda < 0.1 & lambda>0.05 ) | 
            method!= "CCA_graph_rrr" ) %>%
  dplyr::select(method, test_mse, test_cor, val_mse, val_cor)

install.packages("knitr")
install.packages("kableExtra")
library(knitr)
library(kableExtra)
relevant_correlations
latex_table <- kable(relevant_correlations, format = "latex", booktabs = TRUE)


r = 2
test1<-additional_checks(as.matrix(X_transformed), as.matrix(Y_transformed),  
                         S=NULL, 
                         rank=r, kfolds=3, 
                         method.type = "FIT_SAR_CV",
                         lambdax= 10^seq(-3,1, length.out = 30),
                         lambday = c(0))

df = data.frame(1/sqrt(18) * as.matrix(X_transformed)%*% test1$u) #)
#df = data.frame(1/sqrt(18) * as.matrix(Y_transformed)%*%  test1$v) #)
pol = sapply(data$candidate_simplified,
             function(u){candidate_data$Cand_Party_Affiliation[which(candidate_data$candidate_simplified == u)[1]]})
pol = as.character(pol)
df["pol"] = pol
df["name"] = data$candidate_simplified

library(ellipse)
legend_order <- c("REF", "LIB", "REP",
                  "GRE", "IND", "DEM", "CON")
my_colors <- c(  "purple","gold", "red", 
                 "chartreuse2", "orchid1", "dodgerblue",
                 "orchid4"
                 # "lightblue", "lightblue3","cyan", "dodgerblue", "dodgerblue4",
                 #"navyblue", #"cyan", 
                 #"dodgerblue"
)

labels_n <-    c("Reform", "Libertarian",
                 "Republican", "Green", "Independent",
                 "Democratic", "Constitution"
)

ellipse.level =0.8
theme_set(theme_bw(base_size = 24))
ggplot(df, aes(x=X1, y=X2))+
  geom_point(aes( colour=pol), size = 3)+
  geom_text(aes(label = name,  colour=pol), vjust = -0.4, show.legend = FALSE)+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(pol=="LIB")], 
                                                df$X2[which(pol=="LIB")]))), 
                                    centre=colMeans(t(rbind(df$X1[which(pol=="LIB")], 
                                                            df$X2[which(pol=="LIB")]))),
                                    level = ellipse.level
  )),
  aes(x=x, y=y, colour="LIB")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(pol=="REP")], 
                                                df$X2[which(pol=="REP")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(pol=="REP")], 
                          df$X2[which(pol=="REP")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="REP")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(pol=="DEM")], 
                                                df$X2[which(pol=="DEM")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(pol=="DEM")], 
                          df$X2[which(pol=="DEM")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="DEM")) +
  geom_path(data=data.frame(ellipse(cov(t(rbind(df$X1[which(pol=="GRE")], 
                                                df$X2[which(pol=="GRE")])
  )), 
  centre=colMeans(t(rbind(df$X1[which(pol=="GRE")], 
                          df$X2[which(pol=="GRE")]))),
  level = ellipse.level
  )),
  aes(x=x, y=y, colour="GRE")) +
  guides(colour = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid", 
                                                                "solid", "solid",
                                                                "solid", "solid")))) +
  xlab("CD-1") + 
  ylab("CD-2")+
  labs(colour = "Political\nAffiliation")+
  guides(colour = guide_legend(override.aes = list(size = 2)))








ggplot(df, aes(x=X1, y=X2, colour=pol))+
  geom_point()+
  stat_ellipse(geom = "polygon", alpha = 0.) +
  geom_text(aes(label = name), vjust = -1)+
  theme_bw() +
  xlab("CD-1") + 
  ylab("CD-2")+
  labs(colour = "Political\nAffiliation")

Uhat_comp = data.frame(test1$u)
Uhat_comp["state"] = str_to_lower(colnames(X))
map_data <- map_data("state")

merged_data <- merge(map_data, Uhat_comp , by.x = "region", by.y = "state")
# Plot
ggplot() +
  geom_polygon(data = merged_data, aes(x = long, y = lat, group = group, 
                                       fill = X1), colour = "grey") +
  coord_fixed(1.3) +
  labs(fill = "Value") +
  scale_fill_gradient2(mid = "white", high = "yellow", low = "green",
                       limits = c(-0.7, 0.7))

ggplot() +
  geom_polygon(data = merged_data, aes(x = long, y = lat, group = group, 
                                       fill = X2), colour = "grey") +
  coord_fixed(1.3) +
  labs(fill = "Value") +
  scale_fill_gradient2(mid = "white", high = "yellow", low = "green",
                       limits = c(-0.7, 0.8))


Vhat_comp = test1$v
rownames(Vhat_comp) = colnames(Y)
print(Vhat_comp)
heatmap(Vhat_comp)
df_V = data.frame(Vhat_comp)
colnames(df_V) = c("CD-1", "CD-2")
df_V["question"] = colnames(Y)
df_V = pivot_longer(df_V, cols=-c("question"))



# Example new labels
new_labels <- c("Should Abortion Remain  Legal?" , 
                "Should the Death Penalty\nBe Allowed?",
                "Should Former Felons Be\nAllowed to Vote?",
                "Should Federal Taxes Be Increased?" ,
                "Should the US Expand Its \n Nuclear Power?",
                "Are More Regulations On Guns Needed?",
                "Should the US Build a Fence\nAlong the US-Mexico Border?",
                "Is Obamacare Good for America?",
                "Are humans substantially responsible\nfor global climate change?",
                "Should the US tax carbon emissions?")

# Assuming your questions are in a column named 'question'
old_labels <- unique(df_V$question)
label_mapping <- setNames(new_labels, old_labels)
df_V <- df_V %>% mutate(question_new = label_mapping[question])
theme_set(theme_bw(base_size = 18))
ggplot(df_V, aes(x = value, y = reorder(question_new, value), fill=name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ name) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(y = "Question", x = "Loading Value") + 
  labs(fill = "Canonical\nDirection")

plot(log10(final$resultsx$lambda[which(final$resultsx$rmse<1e5)]),
     final$resultsx$rmse[which(final$resultsx$rmse<1e5)])

selected_Y <- c(1:10)#[-c(2, 3,5,6,7,9,10)]

index = 1:nrow(X)#which(apply(X, 1, sum) > 1)



################## CCA
lambda_opt = 1
final = CCA_graph_rrr(as.matrix(X_transformed), as.matrix(Y_transformed),  
                      Gamma, 
                      Sx=NULL, Sy=NULL, Sxy = NULL,
                      lambda = lambda_opt, 
                      Kx=NULL, r=2,
                      rho=1, niter=2 * 1e4,
                      do.scale = FALSE, lambda_Kx=0,
                      thresh=1e-6,
                      LW_Sy = TRUE)


Uhat_comp = data.frame(final$U)
#### Im not sure how to interpret the results
Uhat_comp["state"] = str_to_lower(colnames(X))
Vhat_comp = final$V
rownames(Vhat_comp) = colnames(Y)
print(Vhat_comp)
heatmap(Vhat_comp)

rownames(Vhat_comp) = colnames(Y)
print(Vhat_comp)
heatmap(Vhat_comp)
df_V = data.frame(Vhat_comp)
colnames(df_V) = c("CD-1", "CD-2")
df_V["question"] = colnames(Y)
df_V = pivot_longer(df_V, cols=-c("question"))

df = data.frame(1/sqrt(18) * as.matrix(X_transformed)%*% final$U) #)
df = data.frame(1/sqrt(18) * as.matrix(Y_transformed)%*%  final$V) #)
pol = sapply(data$candidate_simplified,
             function(u){candidate_data$Cand_Party_Affiliation[which(candidate_data$candidate_simplified == u)[1]]})
pol = as.character(pol)
#df["pol"] = pol
df["name"] = data$Title

library(ellipse)
legend_order <- c("REF", "LIB", "REP",
                  "GRE", "IND", "DEM", "CON")
my_colors <- c(  "purple","gold", "red", 
                 "chartreuse2", "orchid1", "dodgerblue",
                 "orchid4"
                 # "lightblue", "lightblue3","cyan", "dodgerblue", "dodgerblue4",
                 #"navyblue", #"cyan", 
                 #"dodgerblue"
)

labels_n <-    c("Reform", "Libertarian",
                 "Republican", "Green", "Independent",
                 "Democratic", "Constitution"
)

ellipse.level =0.8
ggplot(df, aes(x=X1, y=X2))+
  geom_point( size = 3)+
  geom_text(aes(label = name), vjust = -0.4, show.legend = FALSE)+
  #scale_color_manual(values = my_colors, breaks = legend_order,
  #                   labels = labels_n) +
  xlab("CD-1") + 
  ylab("CD-2")+
  labs(colour = "Political\nAffiliation")+
  guides(colour = guide_legend(override.aes = list(size = 2)))


# Example new labels
new_labels <- c("Should Abortion Remain  Legal?" , 
                "Should the Death Penalty\nBe Allowed?",
                "Should Former Felons Be\nAllowed to Vote?",
                "Should Federal Taxes Be Increased?" ,
                "Should the US Expand Its \n Nuclear Power?",
                "Are More Regulations On Guns Needed?",
                "Should the US Build a Fence\nAlong the US-Mexico Border?",
                "Is Obamacare Good for America?",
                "Are humans substantially responsible\nfor global climate change?",
                "Should the US tax carbon emissions?")

# Assuming your questions are in a column named 'question'
old_labels <- unique(df_V$question)
label_mapping <- setNames(new_labels, old_labels)
df_V <- df_V %>% mutate(question_new = label_mapping[question])

ggplot(df_V, aes(x = value, y = reorder(question, value), fill=name)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ name) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  #theme_bw() +
  labs(y = "Question", x = "Response Intensity",
       fill = "Canonical\nDirection") 

relggplot(df_V)+
  geom_tile(aes(x = question, y=name, fill = 100 * value)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Canonical Component Loadings") + xlab('Question') +
  labs(fill = "Loading\nMass (%)")

map_data <- map_data("state")


test = data.frame(score= colMeans(X[which(pol == "REP"),]),
                  state = str_to_lower(colnames(X)))
merged_data <- merge(map_data, Uhat_comp , by.x = "region", by.y = "state")
# Plot
ggplot() +
  geom_polygon(data = merged_data, aes(x = long, y = lat, group = group, 
                                       fill = X1)) +
  coord_fixed(1.3) +
  labs(fill = "Value") +
  scale_fill_gradient2(mid = "white", high = "blue", low = "red")

ggplot() +
  geom_polygon(data = merged_data, aes(x = long, y = lat, group = group, 
                                       fill = X2)) +
  coord_fixed(1.3) +
  labs(fill = "Value") +
  scale_fill_gradient2(mid = "white", high = "blue", low = "red",
                       limits = c(-1, 1.3))

ggplot() +
  geom_polygon(data = merged_data, aes(x = long, y = lat, group = group, 
                                       fill = X2)) +
  coord_fixed(1.3) +
  labs(fill = "Your Value") +
  scale_fill_gradient(low = "blue", high = "red", 
                      limits = c(-0.5, 1))

ggplot() +
  geom_polygon(data = merged_data, aes(x = long, y = lat, group = group, 
                                       fill = X3)) +
  coord_fixed(1.3) +
  labs(fill = "Your Value") +
  scale_fill_gradient(low = "blue", high = "red", 
                      limits = c(-2, 2))

selected_Y
index = 1:18
shrink.rcc<- rcc( as.matrix(X_transformed)[index,],
                  as.matrix(Y)[index,selected_Y], ncomp=r, method = 'shrinkage') 
shrink.rcc$loadings$X = final$U #Uhat_comp
shrink.rcc$variates$X = 1/sqrt(18) * as.matrix(X_transformed)[index,] %*% final$U#XU_comp
shrink.rcc$variates$Y = 1/sqrt(18) *  as.matrix(Y)[index,selected_Y] %*% final$V##YV_comp
#Vhat_comp = data.frame(Vhat_comp)
#rownames(Vhat_comp) = colnames(Y_transformed)
shrink.rcc$loadings$Y =  final$V

pol = sapply(data$candidate_simplified,
             function(u){candidate_data$Cand_Party_Affiliation[which(candidate_data$candidate_simplified == u)[1]]})
pol = as.character(pol)

df = data.frame(1/sqrt(18) * as.matrix(X_transformed)%*% final$U) #)
df = data.frame(1/sqrt(18) * as.matrix(Y_transformed)%*% final$V) #)
df["pol"] = pol
df["name"] = data$candidate_simplified
ggplot(df, aes(x=X1, y=X2, colour=pol))+
  geom_point()+
  stat_ellipse(geom = "polygon", alpha = 0.) +
  geom_text(aes(label = name), vjust = -1)+
  theme_bw() +
  xlab("CC-1") + 
  ylab("CC-2")+
  labs(colour = "Political\nAffiliation")



plotIndiv(shrink.rcc, comp = c(1,2), 
          ind.names = data$candidate_simplified[index],
          ellipse = TRUE,  # plot using the ellipsesp
          rep.space = "X-variate", 
          group = pol[index],
          legend = TRUE, title = '(a) Our method  in XY-space',
          point.lwd = 1)

heatmap(final$V
)

selected_Y <- c(1:10)[-c(2, 3,5, 6,9,10)]

test1<-additional_checks(as.matrix(X_transformed), as.matrix(Y_transformed)[, selected_Y],  
                         S=NULL, 
                         rank=r, kfolds=3, 
                         method.type = "FIT_SAR_CV",
                         lambdax= 10^seq(-3,1, length.out = 30),
                         lambday = 10^seq(-3,1, length.out = 30))

Uhat_comp = data.frame(test1$u)
Uhat_comp["state"] = str_to_lower(colnames(X))
Vhat_comp = test1$v
rownames(Vhat_comp) = colnames(Y)[selected_Y]
print(Vhat_comp)
heatmap(Vhat_comp)

map_data <- map_data("state")

test = data.frame(score= colMeans(X[which(pol == "DEM"),]),
                  state = str_to_lower(colnames(X)))
merged_data <- merge(map_data,  , by.x = "region", by.y = "state")
# Plot
ggplot() +
  geom_polygon(data = merged_data, aes(x = long, y = lat, group = group, 
                                       fill = score)) +
  coord_fixed(1.3) +
  labs(fill = "Your Value") +
  scale_fill_gradient(low = "blue", high = "red")



df = data.frame(as.matrix(Y)[index,] %*% test1$v) #)
df = data.frame(as.matrix(X_transformed)[index,] %*% test1$u) #)
df["pol"] = pol
df["name"] = data$candidate_simplified[index]
ggplot(df, aes(x=X2, y=X3, colour=pol))+
  geom_point()+
  geom_text(aes(label = name), vjust = -1)+
  theme_bw()

shrink.rcc<- rcc( as.matrix(X_transformed)[index,],
                  as.matrix(Y)[index,selected_Y], ncomp=r, method = 'shrinkage') 
shrink.rcc$loadings$X =  test1$u #Uhat_comp
shrink.rcc$variates$X = as.matrix(X_transformed)[index,] %*% test1$u #XU_comp
shrink.rcc$variates$Y =  as.matrix(Y)[index,selected_Y] %*%  test1$v##YV_comp
#Vhat_comp = data.frame(Vhat_comp)
#rownames(Vhat_comp) = colnames(Y_transformed)
shrink.rcc$loadings$Y =   test1$v
plotIndiv(shrink.rcc, comp = c(1,2), 
          ind.names = data$candidate_simplified[index],
          ellipse = TRUE,  # plot using the ellipses
          rep.space = "XY-variate", 
          group = pol[index],
          legend = TRUE, title = '(a) Our method  in XY-space',
          point.lwd = 1)

df = data.frame(as.matrix(Y)[index,] %*%  test1$v) #)
df["pol"] = pol
df["name"] = data$candidate_simplified[index]
ggplot(df, aes(x=X1, y=X2, colour=pol))+
  geom_point()+
  geom_text(aes(label = name), vjust = -1)



Y_tilde = as.matrix(Y)[, selected_Y]
representation = data.frame(Y_tilde %*% test1$v)

pol = sapply(data$candidate_simplified,
             function(u){candidate_data$Cand_Party_Affiliation[which(candidate_data$candidate_simplified == u)[1]]})
pol = as.character(pol)
representation['color'] = pol
representation['candidate'] = data$candidate_simplified
ggplot(representation) +
  geom_point(aes(x= X1,
                 y=X2,
                 colour=color))

Vhat_comp = test1$v

rownames(Vhat_comp) = colnames(Y)[selected_Y]
Vhat_comp
plot(log(final$rmse[which(final$rmse<1e6)]))

source('experiments/sparse_CCA/experiment_functions.R')
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')

final =alternative_method(as.matrix(X_transformed), as.matrix(Y),  
                          Gamma, 
                          Sx=NULL, Sy=NULL, Sxy = NULL,
                          lambda = 0.1, Kx=NULL, r,
                          rho=1, niter=1e4,
                          scale = FALSE, lambda_Kx=0, 
                          LW_Sy = TRUE)

final$rmse
Uhat_comp = data.frame(final$u)
Vhat_comp = final$vfinal


map_data <- map_data("state")

merged_data <- merge(map_data, Uhat_comp , by.x = "region", by.y = "state")
# Plot
ggplot() +
  geom_polygon(data = merged_data, aes(x = long, y = lat, group = group, 
                                       fill = X1)) +
  coord_fixed(1.3) +
  labs(fill = "Your Value") +
  scale_fill_gradient(low = "blue", high = "red", 
                      limits = c(-0.5, 2))


final = CCA_graph_rrr(as.matrix(X_transformed), as.matrix(Y),  
                      Gamma, 
                      Sx=NULL, Sy=NULL, Sxy = NULL,
                      lambda = 0.001, 
                      Kx=NULL, r=r,
                      rho=1, niter=1e4,
                      scale = TRUE, lambda_Kx=0, 
                      LW_Sy = TRUE,
                      thresh=1e-5)


index = which(apply(X, 1, sum)>10)
final = CCA_graph_rrr(as.matrix(X_transformed), as.matrix(Y_transformed),  
                      Gamma, 
                      Sx=NULL, Sy=NULL, Sxy = NULL,
                      lambda = 0.01, 
                      Kx=NULL, r=2,
                      rho=10, niter=1e5,
                      lambda_Kx=0, 
                      LW_Sy = FALSE,
                      thresh=1e-5)


Uhat_comp = data.frame(final$U)
Uhat_comp["state"] = str_to_lower(colnames(X))
Vhat_comp = final$V
rownames(Vhat_comp) = colnames(Y)[-c(10)]
heatmap(Vhat_comp)
i = 1
# Get map data
library(maps)
map_data <- map_data("state")

merged_data <- merge(map_data, Uhat_comp , by.x = "region", by.y = "state")
# Plot
ggplot() +
  geom_polygon(data = merged_data, aes(x = long, y = lat, group = group, 
                                       fill = X2)) +
  coord_fixed(1.3) +
  labs(fill = "Your Value")

color_palette <- colorRampPalette(c("blue", "red"))(49)
# Assign colors based on the continuous variable
node_colors <- color_palette[cut(Uhat_comp[,i], breaks = length(color_palette), labels = FALSE)]
V(g)$color = node_colors


plot.igraph(g,
            vertex.label = V(g)$name,  # Node labels
            vertex.size = 10,  # Adjust node size
            #vertex.color =,  # Node color
            #palette = "viridis",
            edge.color = "grey",  # Edge color
            vertex.label.color = "black",  # Node label color
            edge.label.color = "red",  # Edge label color
            vertex.label.dist = 1.5,  # Distance of labels from nodes
            edge.label.cex = 0.8,  # Edge label size
            vertex.label.cex = 1.2  # Node label size
)


library(mixOmics)
shrink.rcc<- rcc( as.matrix(X_transformed),
                  as.matrix(Y), ncomp=r, method = 'shrinkage') 
map_data <- map_data("state")

merged_data <- merge(map_data, shrink.rcc$loadings$X, by.x = "region", by.y = "state")
# Plot
ggplot() +
  geom_polygon(data = merged_data, aes(x = long, y = lat, group = group, 
                                       fill = X1)) +
  coord_fixed(1.3) +
  labs(fill = "Your Value")

i=1
node_colors <- color_palette[cut(shrink.rcc$loadings$X[,i], breaks = length(color_palette), labels = FALSE)]
V(g)$color = node_colors
plot.igraph(g,
            vertex.label = V(g)$name,  # Node labels
            vertex.size = 10,  # Adjust node size
            #vertex.color =,  # Node color
            #palette = "viridis",
            edge.color = "grey",  # Edge color
            vertex.label.color = "black",  # Node label color
            edge.label.color = "red",  # Edge label color
            vertex.label.dist = 1.5,  # Distance of labels from nodes
            edge.label.cex = 0.8,  # Edge label size
            vertex.label.cex = 1.2  # Node label size
)


shrink.rcc<- rcc( as.matrix(X_transformed),
                  as.matrix(Y), ncomp=r, method = 'shrinkage') 
shrink.rcc$loadings$X = Uhat_comp #Uhat_comp
shrink.rcc$variates$X = as.matrix(X_transformed) %*% Uhat_comp#XU_comp
shrink.rcc$Y = as.matrix(Y_transformed) %*% Vhat_comp##YV_comp
#Vhat_comp = data.frame(Vhat_comp)
#rownames(Vhat_comp) = colnames(Y_transformed)
shrink.rcc$loadings$Y =  Vhat_comp


pol = sapply(data$candidate_simplified,
             function(u){candidate_data$Cand_Party_Affiliation[which(candidate_data$candidate_simplified == u)[1]]})
pol = as.character(pol)
plotIndiv(shrink.rcc, comp = c(1,2), 
          ind.names = data$candidate_simplified,
          ellipse = TRUE,  # plot using the ellipses
          rep.space = "X-variate", 
          group = pol,
          legend = TRUE, title = '(a) Our method  in XY-space',
          point.lwd = 1)

######################## 
########################
########################
#### Analysis finances

##### Analysis of politicians scores vs opinions on certain questions

data = merge(X,
             candidate_data,
             by = c("year", "candidate_simplified"))


state_names = colnames(X)[3:(ncol(X)-1)]
X = data[, state_names]

colnames(candidate_data)
numeric_columns <- colnames(candidate_data)[sapply(candidate_data, is.numeric)]
numeric_columns <- numeric_columns[1:(length(numeric_columns)-1)] ### remove "year"
#numeric_columns = colnames(Y2)[3:ncol(Y2)]
Y = data[, numeric_columns]
numeric_columns <- numeric_columns[which(apply(Y, 2, function(x){mean(x>0)}) > 0.1)]
numeric_columns <- c(23, 9, 15, 2, 1, 8, 14, 11, 4, 7, 12)
Y <- replace(Y, is.na(Y), 0)
pca_result <- prcomp(Y[, numeric_columns], center = TRUE, scale. = TRUE)


### keep only small number

### party affiliation as a label for analysis
#Y = Y[, numeric_columns]

log_signed <- function(p) {
  # Adjust values exactly equal to 0 or 1
  sign(p) * log(abs(p))
}

#index = which(apply(X,1, sum) * 100 > 1)
X_filter = X[index, ]
Y_filter = Y[index, ]

index_col_Z = which(colnames(X) %in% c("ALASKA" ,    "HAWAII"   ))
X = X[, -index_col_Z]
# Apply the logit function with smoothing to all columns
X_transformed <- X %>%
  mutate(across(everything(), logit_smooth))
X_transformed <- X_transformed %>% 
  mutate(across(everything(), scale))



Y_transformed <- Y_filter %>%
  mutate(across(everything(), log_signed))
Y_transformed[is.na(Y_transformed)]=0
Y_transformed <- Y %>% 
  mutate(across(everything(), scale))
hist(Y_transformed$Total_Disbursement)
# View the scaled dataframe
write_csv(as.data.frame(as.matrix(X_transformed)), file="~/Downloads/politicians_scores.csv")
write_csv(as.data.frame(as.matrix(Y_transformed)), "~/Downloads/politicians_data.csv")

X_transformed = as.data.frame(as.matrix(X_transformed))
Y_transformed = as.data.frame(as.matrix(Y_transformed))
#### Download the adjacency matrix
A = read_csv("~/Downloads/state_adjacency.csv")
A = A[,2:ncol(A)]
### erase Hawaii
ind_hawaii = which(colnames(A) == "HI")
A = A[-c(ind_hawaii), -c(ind_hawaii)]

library(igraph)
# Create a graph from the adjacency matrix
g <- graph_from_adjacency_matrix(as.matrix(A), mode = "undirected")

plot(g,
     vertex.label = V(g)$name,  # Node labels
     vertex.size = 10,  # Adjust node size
     vertex.color = "skyblue",  # Node color
     edge.color = "grey",  # Edge color
     vertex.label.color = "black",  # Node label color
     edge.label.color = "red",  # Edge label color
     vertex.label.dist = 1.5,  # Distance of labels from nodes
     edge.label.cex = 0.8,  # Edge label size
     vertex.label.cex = 1.2  # Node label size
)

Gamma <- get_edge_incidence(g, weight = 1)
#### Now apply the algorithm)


p = ncol(X_transformed)
q = ncol(Y_transformed)
r = 5
test1<- additional_checks(X_transformed,
                          Y_transformed, S=NULL, 
                          rank=r, kfolds=10, 
                          method.type = "FIT_SAR_CV",
                          lambdax= 10^seq(-3,0.5, length.out = 10),
                          lambday = c(0, 0))
test1$u = test1$u %*% sqrtm(t(test1$u) %*% cov(X) %*%test1$u )$Binv
test1$v = test1$v %*% sqrtm(t(test1$v) %*% cov(Y) %*%test1$v )$Binv
Uhat_comp = matrix(0, p, r)
index_zeros = which(apply(abs(test1$u), 1, sum) >0)
t = (varimax(test1$u[index_zeros, ], normalize=FALSE))$loadings
for (i in 1:r){
  print(i)
  Uhat_comp[index_zeros,i] = t[,i]
}
XU_comp = as.matrix(X_transformed) %*% Uhat_comp
Vhat_comp = matrix(0, q, r)
index_zeros = which(apply(abs(test1$v), 1, sum) >0)
t = (varimax(test1$v[index_zeros, ], normalize=FALSE))$loadings
for (i in 1:r){
  print(i)
  Vhat_comp[index_zeros,i] = t[,i]
}
YV_comp =  as.matrix(Y_transformed) %*% Vhat_comp
Uhat_comp = data.frame(Uhat_comp)
rownames(Uhat_comp) = colnames(X)

shrink.rcc<- rcc( as.matrix(X_transformed),
                  as.matrix(Y_transformed), ncomp=r, method = 'shrinkage') 

shrink.rcc$loadings$X = test1$u  #Uhat_comp
shrink.rcc$variates$X = as.matrix(X_transformed) %*% test1$u #XU_comp
shrink.rcc$Y = as.matrix(Y_transformed) %*% test1$v##YV_comp
#Vhat_comp = data.frame(Vhat_comp)
#rownames(Vhat_comp) = colnames(Y_transformed)
shrink.rcc$loadings$Y =  test1$v


pol = sapply(data$candidate_simplified,
             function(u){candidate_data$Cand_Party_Affiliation[which(candidate_data$candidate_simplified == u)[1]]})
pol = as.character(pol)
plotIndiv(shrink.rcc, comp = c(1,2), 
          ind.names = data$candidate_simplified,
          ellipse = TRUE,  # plot using the ellipses
          rep.space = "XY-variate", 
          group = pol,
          legend = TRUE, title = '(a) Our method  in XY-space',
          point.lwd = 1)




final =CCA_graph_rrr(X, Y,  Gamma, 
                     Sx=NULL, Sy=NULL,
                     lambda =0.1, Kx, r,
                     scale = FALSE, lambda_Kx=0, 
                     LW_Sy = FALSE)
library(mixOmics)
colbar <- rainbow(50)


rank = sort(Uhat_comp[,1], index.return=T)$ix
rank_x = 1:49
true_rank = rep(0,49)
true_rank[rank] = rank_x
V(g)$color = true_rank


Uhat_comp = matrix(0, p, r)
Vhat_comp = matrix(0, q, r)
t = (varimax(final$U, normalize=TRUE))
Uhat_comp = final$U %*% t$rotmat
t = (varimax(final$V, normalize=TRUE))
Vhat_comp = final$V %*% t$rotmat

mypalette<-brewer.pal(n = 50, name = "Greens")
i = 1
color_palette <- colorRampPalette(c("blue", "red"))(49)
# Assign colors based on the continuous variable
node_colors <- color_palette[cut(Uhat_comp[,i], breaks = length(color_palette), labels = FALSE)]
V(g)$color = node_colors
plot.igraph(g,
            vertex.label = V(g)$name,  # Node labels
            vertex.size = 10,  # Adjust node size
            #vertex.color =,  # Node color
            #palette = "viridis",
            edge.color = "grey",  # Edge color
            vertex.label.color = "black",  # Node label color
            edge.label.color = "red",  # Edge label color
            vertex.label.dist = 1.5,  # Distance of labels from nodes
            edge.label.cex = 0.8,  # Edge label size
            vertex.label.cex = 1.2  # Node label size
)


shrink.rcc.nutrimouse <- rcc(X_transformed,
                             Y_transformed, ncomp=2, method = 'shrinkage') 
# examine the optimal lambda values after shrinkage 
shrink.rcc.nutrimouse$lambda 

shrink.rcc<- rcc( as.matrix(X_transformed),
                  as.matrix(Y_transformed), ncomp=r, method = 'shrinkage') 

shrink.rcc$loadings$X = U #Uhat_comp
shrink.rcc$variates$X = X %*% U#XU_comp
shrink.rcc$Y = Y%*% V#YV_comp
#Vhat_comp = data.frame(Vhat_comp)
#rownames(Vhat_comp) = colnames(Y_transformed)
shrink.rcc$loadings$Y =  V #Vhat_comp


df = data.frame(U)
df["state"] = rownames(df)
df["names"] = colnames(X)
ggplot(df, aes(x = X2, y=X3)) + 
  geom_point() + 
  geom_text(aes(label = names), vjust = -1)

plot(Uhat_comp[,1],Uhat_comp[,2])
ggplot()
pol = sapply(data$candidate_simplified,
             function(u){candidate_data$Cand_Party_Affiliation[which(candidate_data$candidate_simplified == u)[1]]})
pol = as.character(pol)
plotIndiv(shrink.rcc, comp = c(1,2), 
          ind.names = data$candidate_simplified,
          ellipse = TRUE,  # plot using the ellipses
          rep.space = "XY-variate", 
          group = pol,
          legend = TRUE, title = '(a) Our method  in XY-space',
          point.lwd = 1)


#### Now let's try our method


data =read_csv("~/Downloads/binarized_new_results_algo_semi_synthetictest.csv")
test = data %>% 
  group_by(Method, Lambda, graph_type) %>%
  summarise_all(mean)
ggplot(test, aes(x=Lambda,
                 y=Accuracy_true_p_pos,
                 colour = Method))+
  geom_line() +
  scale_x_log10() +
  geom_hline(data=test %>% filter(Method == "SSNAL-opt"),
              aes(yintercept=Accuracy_true_p_pos,
                  colour = Method))


data =read_csv("~/Downloads/sim_all_cv_mean_interpolation_5.csv")
test = data %>% 
  group_by(N) %>%
  summarise(time_splsi_mean =  mean(time_splsi),
            time_splsi_sd = sd(time_splsi),
            time_slda_mean =  mean(time_slda),
            time_slda_sd = sd(time_slda),
            time_vanilla_mean =  mean(time_v),
            time_vanilla_sd = sd(time_v),
            vanilla_acc_mean =  mean(vanilla_acc),
            vanilla_acc_sd = sd(vanilla_acc),
            vanilla_err_mean =  mean(vanilla_err),
            vanilla_err_sd = sd(vanilla_err),
            splsi_acc_mean =  mean(splsi_acc),
            splsi_acc_sd = sd(splsi_acc),
            splsi_err_mean =  mean(splsi_err),
            splsi_err_sd = sd(splsi_err),
            slda_acc_mean =  mean(slda_acc),
            slda_acc_sd = sd(slda_acc),
            slda_err_mean =  mean(slda_err),
            slda_err_sd = sd(slda_err)
            )

ggplot(test, aes(x=N,
                 y=vanilla_err_mean,
                 colour = "Standard PLSI"))+
  geom_line(aes(colour = "Standard PLSI"), linewidth=1.2) +
  geom_point(aes(colour = "Standard PLSI"), size=3) +
  geom_ribbon(aes(x =N, 
                  fill = "Standard PLSI",
                  ymin = vanilla_err_mean - vanilla_err_sd,
                  ymax = vanilla_err_mean + vanilla_err_sd), alpha = 0.3,
              show.legend = FALSE)+
  geom_line(aes(x=N,
                y=splsi_err_mean,
                colour = "Spatial PLSI"),linewidth=1.2) +
  geom_point(aes(x=N,
                 y=splsi_err_mean,
                 colour = "Spatial PLSI"),
             size=3) +
  geom_ribbon(aes(x =N, 
                  fill = "Spatial PLSI",
                  colour = "Spatial PLSI",
                  ymin = splsi_err_mean - splsi_err_sd,
                  ymax = splsi_err_mean + splsi_err_sd), alpha = 0.3,
              show.legend = FALSE)+
  geom_line(aes(x=N,
                y=slda_err_mean,
                colour = "Spatial LDA"),linewidth=1.2) +
  geom_point(aes(x=N,
                 y=slda_err_mean,
                 colour = "Spatial LDA"),
             size=3) +
  geom_ribbon(aes(x =N, 
                  fill = "Spatial LDA",
                  colour = "Spatial LDA",
                  ymin = slda_err_mean - slda_err_sd,
                  ymax = slda_err_mean + slda_err_sd), alpha = 0.3,
              show.legend = FALSE)+
  scale_x_log10() +
  ylab(" Error")



ggplot(test, aes(x=N,
                 y=time_vanilla_mean,
                 colour = "Standard PLSI"))+
  geom_line(aes(colour = "Standard PLSI"), linewidth=1.2) +
  geom_point(aes(colour = "Standard PLSI"), size=3) +
  geom_ribbon(aes(x =N, 
                  fill = "Standard PLSI",
                  ymin = time_vanilla_mean - time_vanilla_sd,
                  ymax = time_vanilla_mean + time_vanilla_sd), alpha = 0.3,
              show.legend = FALSE)+
  geom_line(aes(x=N,
                y=time_splsi_mean,
                colour = "Spatial PLSI"),linewidth=1.2) +
  geom_point(aes(x=N,
                 y=time_splsi_mean,
                 colour = "Spatial PLSI"),
             size=3) +
  geom_ribbon(aes(x =N, 
                  fill = "Spatial PLSI",
                  colour = "Spatial PLSI",
                  ymin = time_splsi_mean - time_splsi_sd,
                  ymax = time_splsi_mean + time_splsi_sd), alpha = 0.3,
              show.legend = FALSE)+
  geom_line(aes(x=N,
                y=time_slda_mean,
                colour = "Spatial LDA"),linewidth=1.2) +
  geom_point(aes(x=N,
                 y=time_slda_mean,
                 colour = "Spatial LDA"),
             size=3) +
  geom_ribbon(aes(x =N, 
                  fill = "Spatial LDA",
                  colour = "Spatial LDA",
                  ymin = time_slda_mean - time_slda_sd,
                  ymax = time_slda_mean + time_slda_sd), alpha = 0.3,
              show.legend = FALSE)+
  scale_x_log10() +
  ylab(" Time (in s)")
