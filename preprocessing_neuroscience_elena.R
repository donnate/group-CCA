library(tidyverse)
library(RNifti)
library(oro.nifti)
library(neurobase)

library(tidyverse)
library(igraph)
library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)
setwd("/scratch/cdonnat/group-CCA/")
#store all values above diagonal of connectome matrices in matrix c

activations <- read_csv("data/activation_neuroscience.csv")
behaviour <- read_csv("data/behavior.csv")
group_assignment <- readxl::read_xlsx("data/activation_groups.xlsx", col_names = FALSE)
colnames(group_assignment) <- "group.id"
index_groups = which(group_assignment$group.id!=0)
activations  = activations[,c(1, 1+index_groups)]
groups <- sapply(1:length(unique(group_assignment$group.id[index_groups])),
                 function(g){which(group_assignment$group.id[index_groups] == g)})


new_data = activations %>% inner_join(behaviour, join_by(Row == participant_id ))
X = new_data[, 2:ncol(activations)]
Y = new_data[, c("demo_age",  "bio_sex",  "bas_drive" ,
                 "bas_funseeking",  "bas_rewardrespons",
                 "bis_total", "masq_distress", "masq_anhedonia",   
                 "masq_anxarousal", "panas_positive","panas_negative")]

gender = Y$bio_sex
Y = Y[ , c("bas_drive" ,
           "bas_funseeking",  "bas_rewardrespons",
           "bis_total", "masq_distress", "masq_anhedonia",   
           "masq_anxarousal", "panas_positive","panas_negative")]
females = which(gender == "F")
males = which(gender == "M")
means_f = colMeans(X[females, ])
means_m = colMeans(X[males, ])
X[females, ] = X[females, ] + (means_m - means_f)

means_qf = colMeans(Y[females,  ], na.rm=TRUE)
means_qm = colMeans(Y[males, ], na.rm=TRUE)
Y[females, ] = Y[females, ] + (means_m - means_f)

# Calculate column means, excluding NAs
column_means <- colMeans(Y, na.rm = TRUE)

# Replace NA values with column means
for (i in 1:ncol(Y)){
  Y[is.na(Y[,i]), i] <- column_means[i]
}


##### Split into different 

n = nrow(X)
p = ncol(X)
q = ncol(Y)
do.scale = T
if(q >n){
  X_temp <- X
  X <- Y
  Y <- X_temp
}
if (do.scale){
  X <- scale(X)
  Y <- scale(Y)
}
write_csv(X, "data/activations_X_preprocessed.csv")
write_csv(Y, "data/activations_Y_preprocessed.csv")

folds = createFolds(1:nrow(X), k=16)
write_csv(folds, "data/folds.csv")
