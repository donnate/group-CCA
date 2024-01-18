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

data_presidents <- read_csv("~/Downloads/1976-2020-president.csv")
data_presidents <- data_presidents %>% 
  filter(year > 2007, is.na(candidate) == FALSE) %>%
  mutate(percentage_votes = candidatevotes/ totalvotes)

X = pivot_wider(data_presidents %>% 
                  dplyr::select(year, state, candidate, percentage_votes) %>%
                  group_by(year, state, candidate) %>%
                  summarise(percentage_votes = sum(percentage_votes)) %>%
                  ungroup(), 
                id_cols = c("year", "candidate"),
                names_from = state,
                values_from =  "percentage_votes",
                values_fill=0)
all_candidates = unique(X$candidate)
### Maybe select candidates represented in at least 20 states
X_bis <- mean(apply(X[, 3:ncol(X)] >0, 1, mean) > 0.1)

candidate_data20 <- read_csv("~/Downloads/candidate_summary_2020.csv") %>% 
                        filter(Cand_Office == "P")%>%
  mutate(year = 2020)
candidate_data16 <- read_csv("~/Downloads/candidate_summary_2016.csv") %>% 
  filter(Cand_Office == "P") %>%
  mutate(year = 2016)
candidate_data12 <- read_csv("~/Downloads/candidate_summary_2012.csv") %>% 
  filter(Cand_Office == "P")%>%
  mutate(year = 2012)
candidate_data08 <- read_csv("~/Downloads/candidate_summary_2008.csv") %>% 
  filter(Cand_Office == "P")%>%
  mutate(year = 2008)
candidate_data = rbind(candidate_data20,
                       candidate_data16,
                       candidate_data12,
                       candidate_data08)


extract_format <- function(name) {
  parts <- strsplit(name, ", ")[[1]]
  last_name <- parts[1]
  first_initial <- substr(parts[2], 1, 2)
  return(paste(last_name, first_initial, sep = ", "))
}

X$candidate_simplified = sapply(X$candidate, extract_format)
candidate_data$candidate <- sub("([A-Z]+, [A-Z]+).*", "\\1", candidate_data$Cand_Name)
candidate_data$candidate_simplified <- sapply(candidate_data$Cand_Name, extract_format)
index = which(candidate_data$Cand_Name %in% c("WEST, KANYE DEEZ NUTZ", "TRUMP, DON'T VOTE FOR"))
candidate_data = candidate_data[-index,] #remove  "WEST, KANYE DEEZ NUTZ", "TRUMP, DON'T VOTE FOR"
#### Checks

index = which(X$candidate %in% c( "WHITE, JEROME \"\"JERRY\"\"" ))
X = X[-c(index),]
sort(setdiff(unique(X$candidate_simplified),
             unique(candidate_data$candidate_simplified)
             )
     )

Y2 <- read_csv("~/Downloads/politicians_positions.csv")


##### Analysis of politicians scores vs opinions on certain questions

data = merge(X,
             Y2,
             #candidate_data,
             by = c("year", "candidate_simplified"))


state_names = colnames(X)[3:(ncol(X)-1)]
X = data[, state_names]

colnames(candidate_data)
#numeric_columns <- colnames(candidate_data)[sapply(candidate_data, is.numeric)]
#numeric_columns <- numeric_columns[1:(length(numeric_columns)-1)] ### remove "year"
numeric_columns = colnames(Y2)[3:ncol(Y2)]
Y = data[, numeric_columns]
Y <- replace(Y, is.na(Y), 0)
### keep only small number
#numeric_columns <- numeric_columns[which(apply(Y, 2, function(x){mean(x>0)}) > 0.3)]
### party affiliation as a label for analysis
#Y = Y[, numeric_columns]
#### Find adacency map for
#minX = min(X[X>0])
#minY = min(Y[Y>0])
hist(apply(X, 1, function(x){mean(x>0)}))
dim(X)
dim(Y)

logit_smooth <- function(p) {
  # Adjust values exactly equal to 0 or 1
  minX = 1e-8
  p[p == 0] <- minX
  log(p / (1 - p))
}
#X_filter = X[index, ]
#Y_filter = Y[index, ]

index_col_Z = which(colnames(X) %in% c("ALASKA" ,    "HAWAII"   ))
X = X[, -index_col_Z]
# Apply the logit function with smoothing to all columns
X_transformed <- X %>%
  mutate(across(everything(), logit_smooth))
X_transformed <- X_transformed %>% 
  mutate(across(everything(), scale))


###
Sxy = matrix(0, p, q)
for (i in 1:p){
  for (j in 1:q){
    Sxy[i,j] = polyserial(X_transformed[,i], Y[,j], ML = FALSE, control = list(), 
                          std.err = FALSE, maxcor=.9999, bins=4, start, thresholds=FALSE)
  }
}
Sy = diag(rep(1, q))
for (i in 1:(q-1)){
  for (j in (i+1):q){
    Sy[i,j] = polyserial(Y[,i], Y[,j], ML = FALSE, control = list(), 
                         std.err = FALSE, maxcor=.9999, bins=4, start, thresholds=FALSE)
    Sy[j,i] = Sy[i,j]
  }
}


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
r = 2
p = ncol(X_transformed)
q = ncol(Y_transformed)
final =CCA_graph_rrr(as.matrix(X_transformed), as.matrix(Y_transformed),  
                     Gamma, 
                     Sx=NULL, Sy=Sy, Sxy = Sxy,
                     lambda = 0.05, Kx=NULL, r,
                     scale = FALSE, lambda_Kx=0, 
                     LW_Sy = TRUE)

Uhat_comp = final$U
Vhat_comp = final$V

mypalette<-brewer.pal(n = 11, name = "Spectral")
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


shrink.rcc<- rcc( as.matrix(X_transformed),
                  as.matrix(Y_transformed), ncomp=r, method = 'shrinkage') 

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




  