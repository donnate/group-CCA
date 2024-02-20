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
#store all values above diagonal of connectome matrices in matrix c
dir = '../fMRI-data/data/'
files = list.files("/Users/cdonnat/Documents/group-CCA/experiments/real-data/fMRI-data/data/wm_ntgtvsbl/")
ids_length = 6
c = c()
c_ids = c()

Y = read_csv("/Users/cdonnat/Documents/group-CCA/experiments/real-data/fMRI-data/data/cognitive_data.csv")
#### demean 
prop_missing = apply(Y, 1, function(x){mean(is.na(x))})
selected_index = which(prop_missing < 0.1)
prop_col_missing = apply(Y[selected_index, 1:74], 2, function(x){mean(is.na(x))})
newY = Y[selected_index, 1:74]

Y[is.na(Y)] = 0 #### need to replace it by the mean of the column rather
GM_mask = readnii("/Users/cdonnat/Documents/group-CCA/experiments/real-data/fMRI-data/data/gm_mask020_bin.nii.gz")
GM.d = GM_mask@.Data
GM.d <- reshape2::melt(GM.d)
colnames(GM.d)<- c("x", "y", "z", "GM")


atlas800 = readnii("/Users/cdonnat/Documents/group-CCA/experiments/real-data/fMRI-data/data/atlases/Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz")
atlas.d = atlas800@.Data
atlas.d <- reshape2::melt(atlas.d)
colnames(atlas.d)<- c("x", "y", "z", "Schaffer_group")
dataset <- c()
for(file in files){
  print(file)
  t1 = readnii(paste0("/Users/cdonnat/Documents/group-CCA/experiments/real-data/fMRI-data/data/wm_ntgtvsbl/", file))
  test = t1@.Data
  test[is.na(test)] = 0
  test <- reshape2::melt(test)
  colnames(test) = c('x', 'y', 'z', 'Response')
  test <- test %>%
    left_join(atlas.d, by = c('x', 'y', 'z'))
  test = cbind(test, atlas.d$Schaffer_group)
  #test["voxel"] = 1:nrow(test)
  
  t = test %>%
    dplyr::select(Response, Schaffer_group, x, y, z) %>%
    group_by(Schaffer_group) %>% summarise_all(mean)
  t["id"] = file
  
  dataset = rbind(dataset, t)

}

#### find  6 nearest neighbors.
#### This is probably the thing to do, but is currently too slow
# positions = atlas.d %>% 
#   filter(Schaffer_group > 0) 
# #### Find all contacts
# contacts <- c()
# contacts_s <-c()
# for (i in 1:nrow(positions)){
#   print(i)
#   ### Find all 6 nearest neighbor groups
#   ind_x = which((positions$x == (positions$x[i]-1)) & (positions$y == (positions$y[i])) &(positions$z == (positions$z[i])))
#   if (length(ind_x) >0){
#     contacts <- rbind(contacts, c(positions$Schaffer_group[i], positions$Schaffer_group[ind_x]))
#   }
#   ind_x = which((positions$x == (positions$x[i]+1)) & (positions$y == (positions$y[i])) &(positions$z == (positions$z[i])))
#   if (length(ind_x) >0){
#     contacts <- rbind(contacts, c(positions$Schaffer_group[i], positions$Schaffer_group[ind_x]))
#   }
#   
#   
#   ind_y = which((positions$x == (positions$x[i])) & (positions$y == (positions$y[i] -1)) &(positions$z == (positions$z[i])))
#   if (length(ind_y) >0){
#     contacts <- rbind(contacts, c(positions$Schaffer_group[i], positions$Schaffer_group[ind_y]))
#   }
#   ind_y = which((positions$x == (positions$x[i])) & (positions$y == (positions$y[i] + 1)) &(positions$z == (positions$z[i])))
#   if (length(ind_y) >0){
#     contacts <- rbind(contacts, c(positions$Schaffer_group[i], positions$Schaffer_group[ind_y]))
#   }
#   
#   ind_z = which((positions$x == (positions$x[i])) & (positions$y == (positions$y[i] )) &(positions$z == (positions$z[i]-1)))
#   if (length(ind_z) >0){
#     contacts <- rbind(contacts, c(positions$Schaffer_group[i], positions$Schaffer_group[ind_y]))
#   }
#   ind_z = which((positions$x == (positions$x[i])) & (positions$y == (positions$y[i])) &(positions$z == (positions$z[i]-1)))
#   if (length(ind_z) >0){
#     contacts <- rbind(contacts, c(positions$Schaffer_group[i], positions$Schaffer_group[ind_y]))
#   }
#   if (i %% 1000 == 0){
#     contacts_temp = data.frame(contacts)
#     colnames(contacts_temp) = c('C1', 'C2')
#     contacts_temp = contacts_temp %>% group_by(C1, C2) %>% summarize(counts=n())
#     contacts_s = rbind(contacts_s, contacts_temp)
#     contacts = c()
#     print(contacts_s)
#   }
# }
# 
# contacts_temp = data.frame(contacts)
# colnames(contacts_temp) = c('C1', 'C2')
# contacts_temp = contacts_temp %>% group_by(C1, C2) %>% summarize(counts=n())
# contacts_s = rbind(contacts_s, contacts_temp)
# contacts = c()
# print(contacts_s)
# 
# contacts_s2 = contacts_s %>% 
#   group_by(C1, C2) %>% 
#   summarize(counts=n())
# #### need to save
# write_csv(contacts_s2, "~/Downloads/contact_brain_atlas_800.csv")

positions = atlas.d %>% 
  filter(Schaffer_group > 0) %>%
  group_by(Schaffer_group) %>%
  summarise_all(mean)

# Load the igraph package
library(igraph)

# Assuming coordMatrix is your matrix of coordinates
# coordMatrix <- matrix(c(...), ncol=3) # Example format

# Compute the distance matrix
distanceMatrix <- as.matrix(dist(positions))

# Find the 6 nearest neighbors for each point
numNeighbors <- 6
edges <- c()
for (i in 1:nrow(distanceMatrix)){
  x = distanceMatrix[i,] 
  edges = rbind(edges,
                t(rbind(rep(i, numNeighbors), order(x)[1:numNeighbors+1])))
}
edges = data.frame(edges)
colnames(edges) <- c("Source", "Target")
# Create a graph from these edges
g <- graph_from_edgelist(as.matrix(edges), directed = FALSE)

# Optionally, plot the graph
plot(g)
# Find the connected components
comp <- components(g)
# Get the number of connected components
num_components <- comp$no
print(num_components)

#threed_coordinates = read.table("/Users/cdonnat/Documents/group-CCA/experiments/real-data/fMRI-data/data/atlases/Schaefer2018_800Parcels_7Networks_order.txt", sep="\t")
#colnames(threed_coordinates) <- c("id", "name", "x", "y", "z")





# dataset = pivot_wider(dataset, id_cols = c("id"),
#                       names_from = "Schaffer_group",
#                       values_from = "Response" )
#### Apply our method
record = Y$record_id[selected_index]
dataset = dataset %>%
  mutate("subject" = sub("_.*", "", id))

remaining_records = intersect(unique(dataset$subject), record)
X = pivot_wider(
  dataset %>% 
    filter(Schaffer_group >0, subject %in% record) %>%
    dplyr::select(Schaffer_group, Response, subject),
  id_cols = "subject",
  names_from = c("Schaffer_group"),
  values_from =  "Response") %>%
arrange(desc(subject)) 

X =X %>%
  dplyr::select(-c(subject))

X_transformed <- scale(X) 


newY <- newY %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))

Y_transformed <- scale(newY %>%
  filter(record_id %in% remaining_records) %>%
  arrange(desc(record_id)) %>%
  dplyr::select(-c(record_id, redcap_event_name))) 

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

Gamma <- get_edge_incidence(g, weight = 1)
Gamma_dagger = pinv(Gamma)
Pi = diag(rep(1, p)) - Gamma_dagger %*% Gamma
r = 2
p = ncol(X_transformed)
q = ncol(Y_transformed)
n = nrow(X_transformed)
final =CCA_graph_rrr(as.matrix(X_transformed), as.matrix(Y_transformed),  
                     Gamma, 
                     Sx=NULL, Sy=NULL, Sxy = NULL,
                     lambda = 0.031622777, Kx=NULL, r,
                     scale = FALSE, lambda_Kx=0, 
                     LW_Sy = TRUE)

final =CCA_graph_rrr.CV(as.matrix(X_transformed), as.matrix(Y_transformed),  
  Gamma,
  r=4, Kx = NULL, lambda_Kx = 0,
  param_lambda=10^seq(-3, 0, length.out = 10),
  kfolds=10, 
  parallelize = FALSE,
  scale=FALSE,
  LW_Sy = TRUE)
plot(final$rmse)
library(plotly)
# Create a graph from these edges

p <- plot_ly(x = positions$x, y = positions$y, 
             z = positions$z, type = 'scatter3d', mode = 'markers')

edges_df <- do.call(rbind, lapply(E(g), function(e) {
  start <- ends(g, e)[1]
  end <- ends(g, e)[2]
  rbind(cbind(positions[start,c("x", "y", "z")], type="start"), cbind(positions[end,c("x", "y", "z")], type="end"))
}))
edges_df$color = NA
i=4
positions$color = final$ufinal[,i]
# Calculate the thresholds
top_10_percent <- quantile(positions$color, 0.90)
bottom_10_percent <- quantile(positions$color, 0.10)


positions$labels = node_name$label_short

# Filter nodes for labeling
label_nodes <- positions %>%
  filter(color>= top_10_percent | color <= bottom_10_percent)
label_nodes$labels <- gsub("^LH_", "", label_nodes$labels)
label_nodes$labels <- gsub("^RH_", "", label_nodes$labels)

p <- plot_ly(data = positions, x = ~x, y = ~y, z = ~z, 
             color= ~color,
             text = ~labels,  # Add node labels
             type = 'scatter3d', mode = 'markers')

# Add edges
p <- p %>% add_trace(
  data = edges_df, 
  x = ~x, y = ~y, z = ~z,
  type = 'scatter3d', mode = 'lines',
  line = list(color = 'grey', width = 2)
)

# Add labels as a separate trace
p <- p %>% add_trace(
  data = label_nodes,
  x = ~x, y = ~y, z = ~z,
  type = 'scatter3d', mode = 'text',
  text = ~labels, # labels column in your positions dataframe
  textposition = 'top center',
  hoverinfo = 'none'
)


final$vfinal

node_name = read_table("/Users/cdonnat/Dropbox/My Mac (Claire’s MacBook Pro)/Downloads/data 3/atlases/Schaefer2018_200Parcels_7Networks_order.txt",
                      col_names = c("nb", "label", "x", "y", "z"))
node_name$label_short <- gsub("7Networks_", "", node_name$label)


node_colors <- color_palette[cut(final$U[,2], breaks = length(color_palette), labels = FALSE)]
V(g)$color = node_colors
plot.igraph(g,
            vertex.label = V(g)$name,  # Node labels
            vertex.size = 10,  # Adjust node size
            #vertex.color =,  # Node color
            #palette = "viridis",
            edge.color = "black",  # Edge color
            vertex.label.color = "black",  # Node label color
            edge.label.color = "red",  # Edge label color
            vertex.label.dist = 0,  # Distance of labels from nodes
            edge.label.cex = 0.8,  # Edge label size
            vertex.label.cex = 1.,  # Node label size
            layout = as.matrix(positions[,c("x", "y", "z")])
)

plot(final$U[,2], final$U[,1])
plot(1:p, final$ufinal[,2])

library(mixOmics)
shrink.rcc<- rcc( as.matrix(X_transformed),
                  as.matrix(Y_transformed), ncomp=r, method = 'ridge') 

i=1
node_colors <- color_palette[cut(shrink.rcc$loadings$X[,i], breaks = length(color_palette), labels = FALSE)]
V(g)$color = node_colors


plot.igraph(g,
            vertex.label = V(g)$name,  # Node labels
            vertex.size = 10,  # Adjust node size
            #vertex.color =,  # Node color
            #palette = "viridis",
            edge.color = "black",  # Edge color
            vertex.label.color = "black",  # Node label color
            edge.label.color = "red",  # Edge label color
            vertex.label.dist = 0,  # Distance of labels from nodes
            edge.label.cex = 0.8,  # Edge label size
            vertex.label.cex = 1.,  # Node label size,
            layout = as.matrix(positions[,c("x", "y", "z")])
)


# Add edges
for(e in E(g)) {
  start <- ends(g, e)[1]
  end <- ends(g, e)[2]
  start_coord <- as.numeric(positions[start,c("x","y", "z")])
  end_coord <- as.numeric(positions[end,c("x","y", "z")])
  p <- add_lines( p,
                 x = c(start_coord[1], end_coord[1]), 
                 y = c(start_coord[2], end_coord[2]), 
                 z = c(start_coord[3], end_coord[3]),
                 line = list(color = 'black', width = 2)
                 )
  
}

p

p = ncol(X_transformed)
q = ncol(Y_transformed)
r = 5
test1<- additional_checks(as.matrix(X_transformed), as.matrix(Y_transformed),  
                          S=NULL, 
                          rank=r, kfolds=3, 
                          method.type = "FIT_SAR_BIC",
                          lambdax= 10^seq(-1,0.5, length.out = 10),
                          lambday = c(0, 0))
test1$u = test1$u %*% sqrtm(t(test1$u) %*% cov(X) %*%test1$u )$Binv
test1$v = test1$v %*% sqrtm(t(test1$v) %*% cov(Y) %*%test1$v )$Binv
Uhat_comp = matrix(0, p, r)


#### Couple of test

ggplot(
  dataset %>% 
  group_by(id) %>%
  summarize(v = mean(Response)),
  aes(x=v))+
  geom_histogram()
##### Maybe there is something to say for the subject effect.
dataset = dataset %>% 
  filter(Schaffer_group >0 ) 

mean_subjects = dataset %>% 
  group_by(id) %>%
  summarize(v = mean(Response)) %>%
  arrange(v)
ggplot(
  dataset %>% 
    group_by(voxel) %>%
    summarize(v = mean(Response)),
  aes(x=v))+
  geom_histogram()

median_subjects = dataset %>% 
  filter(Schaffer_group >0 ) %>%
  group_by(id) %>%
  summarize(v = median(Response)) %>%
  arrange(v)

ggplot(median_subjects,
  aes(x=v))+
  geom_histogram()

n_voxels = length(unique(dataset$voxel))
ggplot(dataset %>% filter(id  %in% sample(unique(dataset$id), 10),
       voxel %in% sample(unique(dataset$voxel), 200)))+
  geom_point(aes(x=voxel, y=Response, colour=id))
ggplot(dataset %>% filter(id  %in% c("CONN137_con_0007.nii", "CONN015_con_0007.nii", "CONN126_con_0007.nii",
                                     "CONN052_con_0007.nii", "CONN079_con_0007.nii", "CONN131_con_0007.nii"),
                          voxel %in% sample(unique(dataset$voxel), 500)))+
  geom_point(aes(x=voxel, y=Response, colour=id))


dataset2 = dataset %>% 
  group_by(id) %>%
  mutate(y = Response- median(Response),
         y.m = Response- mean(Response)) %>%
  ungroup()

##### Not sure that I should actually demean => maybe that's an unreasonable assumption
##### We'd be taking out the strong components....  Not if we demedian though

ggplot(dataset2 %>% filter(id  %in% c("CONN137_con_0007.nii", "CONN015_con_0007.nii", "CONN126_con_0007.nii",
                                     "CONN052_con_0007.nii", "CONN079_con_0007.nii", "CONN131_con_0007.nii"),
                          voxel %in% sample(unique(dataset$voxel), 500)))+
  geom_point(aes(x=voxel, y=y, colour=id))

ggplot(dataset2 %>% filter(id  %in% c("CONN137_con_0007.nii", "CONN015_con_0007.nii", "CONN126_con_0007.nii",
                                      "CONN052_con_0007.nii", "CONN079_con_0007.nii", "CONN131_con_0007.nii")))+
  geom_histogram(aes(x=y, fill=id))

ggplot(dataset2 %>% filter(id  %in% c("CONN137_con_0007.nii", "CONN015_con_0007.nii", "CONN126_con_0007.nii",
                                      "CONN052_con_0007.nii", "CONN079_con_0007.nii", "CONN131_con_0007.nii")))+
  geom_boxplot(aes(x=id, y=Response, fill=id))

ggplot(dataset2 %>% filter(id  %in% c("CONN137_con_0007.nii", "CONN015_con_0007.nii", "CONN126_con_0007.nii",
                                      "CONN052_con_0007.nii", "CONN079_con_0007.nii", "CONN131_con_0007.nii")))+
  geom_boxplot(aes(x=id, y=y, fill=id))

#### 

#### We need to recover a graph

t = dataset2 %>% 
  dplyr::select(Schaffer_group, voxel) %>%
  group_by(Schaffer_group) %>% mutate(n=n_distinct(voxel))


p = exp(sum(log(dim(t1))))
G = make_lattice(dim(t1))  %>%
  set_vertex_attr("name", value = 1:p)

node_list = V(G)$name
#node_list = node_list[which(as.vector(atlas.d)>0)]

G = delete_vertices(G, node_list[which(as.vector(atlas.d)==0)])
#### This gets us 5 different components, interestingly
n_nodes = length(which(as.vector(atlas.d)==0))
ego_size(
  G,
  order = 2,
  nodes = V(G)
)

t = dataset2 %>% 
     dplyr::select(Schaffer_group, voxel) %>%
     group_by(Schaffer_group) %>% 
     tally()

##### it needs to be sparse in some dimension...
#### Maybe the jumps between areas need to be sparse




X = pivot_wider(dataset2 %>% 
                  dplyr::select(id, y, voxel), id_cols=c(id),
                names_from=voxel, values_from=y)

X.m = dataset2 %>%
  group_by(id, Schaffer_group) %>%
  summarise_all(mean)
X.m = pivot_wider(X.m %>% 
                  dplyr::select(id, y, Schaffer_group), id_cols=c(id),
                  names_from=Schaffer_group, values_from=y)
#### Should also be sparse????
#### Let's try our Sparse CCA approach?
#### 


#### Let's try our Group CCA approach?
#### correlation
rho = (t(X)  %*% X)/nrow(X)
#A = as_adjacency_matrix(G, sparse=FALSE)  + diag(rep(1, p))
X_a = X %*% as_adjacency_matrix(G, sparse=FALSE)  + diag(rep(1, p)) ### does not fit into memory
#### This makes things a little bit less sparse and we are going to need to threshold
#### We do the CCA on the smoothed matrix
#### Then essentially we do an 
#### Maybe we need an extra condition
#### This should essentially increase the signal
#### Unfortunately, it renders things less smooth
#### Consequently, we bet on finding the jumps
##### We want to aggregate things by group
sapply(unique(X$id), function(x){str_split(x, "_")[[1]][1]})
write_csv("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/all_brain_data.csv")
