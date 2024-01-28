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
# Compute initial means
means_qf = colMeans(X[females,], na.rm=TRUE)
means_qm = colMeans(X[males,], na.rm=TRUE)

# Calculate the difference in means
mean_diff = means_qm - means_qf

# Adjust the female data
X[females,] = X[females,] + matrix(rep(mean_diff, nrow(X[females,])), nrow=nrow(X[females,]), byrow=TRUE)

# Recompute means to check
new_means_qf = colMeans(X[females,], na.rm=TRUE)
new_means_qm = colMeans(X[males,], na.rm=TRUE)

# Check if the means are now equal
all.equal(new_means_qf, new_means_qm)


# Compute initial means
means_qf = colMeans(Y[females,], na.rm=TRUE)
means_qm = colMeans(Y[males,], na.rm=TRUE)

# Calculate the difference in means
mean_diff = means_qm - means_qf

# Adjust the female data
Y[females,] = Y[females,] + matrix(rep(mean_diff, nrow(Y[females,])), nrow=nrow(Y[females,]), byrow=TRUE)

# Recompute means to check
new_means_qf = colMeans(Y[females,], na.rm=TRUE)
new_means_qm = colMeans(Y[males,], na.rm=TRUE)

# Check if the means are now equal
all.equal(new_means_qf, new_means_qm)

# Replace NA values with column means
column_means <- colMeans(Y, na.rm = TRUE)

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
write_csv(data.frame(X), "data/activations_X_preprocessed.csv")
write_csv(data.frame(Y), "data/activations_Y_preprocessed.csv")

folds = createFolds(1:nrow(X), k=16)
foldVector <- seq(from = 1, to = nrow(X), by = 10)
folds = split(sample(1:nrow(X), nrow(X)), foldVector)
write_csv(t(data.frame(folds)), "data/folds.csv")


