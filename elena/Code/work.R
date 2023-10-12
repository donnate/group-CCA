##### We need a metric to check for the validity of the recovery of these principal components
library(tidyverse)
library(RNifti)
library(oro.nifti)
library(neurobase)
library("scatterplot3d") 

threed_coordinates = read.table("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/atlases/Schaefer2018_600Parcels_7Networks_order.txt", sep="\t")
colnames(threed_coordinates) <- c("id", "name", "x", "y", "z")
atlas = readnii("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/atlases/Schaefer2018_300Parcels_7Networks_order_FSLMNI152_2mm.gz")
atlas.d = atlas@.Data

GM_mask = readnii("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/gm_mask020_bin.nii.gz")
GM.d = GM_mask@.Data
GM.d <- reshape2::melt(GM.d)
colnames(GM.d)<- c("x", "y", "z", "GM")

atlas  = c()
for (i in seq(1000, 100, -100) ){
  atlas_temp = readnii(paste0("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/atlases/Schaefer2018_", as.character(i), 
                         "Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"))
  atlas.d = atlas_temp@.Data
  atlas.d <- reshape2::melt(atlas.d)
  colnames(atlas.d)<- c("x", "y", "z", i)
  if (i ==1000){
    atlas = atlas.d
  }else{
    atlas = atlas %>%
      left_join(atlas.d,
                by = c("x", "y", "z"))
  }
  
}

atlas = atlas %>%
  left_join(GM.d,
            by = c("x", "y", "z"))
atlas["voxel"] = 1:nrow(atlas)
##### Need to construct a distance matrix to group all of these voxels
atlas_agg = atlas %>% 
  dplyr::select(`1000`, `900`, `800`, `700`, `600`,
                `500`, `400`, `300`, `200`, `100`)  %>%
  group_by(`1000`, `900`, `800`, `700`, `600`,
           `500`, `400`, `300`, `200`, `100`) %>%
  summarise(count = n()) 
atlas_agg["id"] = 1:nrow(atlas_agg)


ggplot(atlas_agg %>% dplyr::filter(`1000` > 0), aes(x=count)) +
  geom_histogram() + scale_x_log10()

##### We just want to group regularize 
##### What would this look like?


#atlas.d <- data.frame(atlas.d, row.names= 1:nrow(atlas.d))
atlas.d2 <- (atlas.d)
#colnames(atlas.d2)<- c("x", "y", "z", "value")
atlas.d2 <- atlas.d2 %>%
  mutate(neighbor1 = ifelse(x > min(atlas.d2$x), x-1, NA),
         neighbor2 = ifelse(x < max(atlas.d2$x), x+1, NA),
         neighbor3 = ifelse(y > min(atlas.d2$y), y-1, NA),
         neighbor4 = ifelse(y < max(atlas.d2$y), y+1, NA),
         neighbor5 = ifelse(z> min(atlas.d2$z), z-1, NA),
         neighbor6 = ifelse(z < max(atlas.d2$z), z+1, NA)) 
atlas.d2["id"] = 1:nrow(atlas.d2)
atlas.d2 <-  atlas.d2 %>%
  left_join(atlas.d2 %>% dplyr::select(id, x, y, z) %>% rename(neighbor1=x,
                                                               neighbor1_id = id),
            by=c("neighbor1", "y", "z"))
atlas.d2 <-  atlas.d2 %>%
  left_join(atlas.d2 %>% dplyr::select(id, x, y, z) %>% rename(neighbor2=x,
                                                               neighbor2_id = id),
            by=c("neighbor2", "y", "z"))
atlas.d2 <-  atlas.d2 %>%
  left_join(atlas.d2 %>% dplyr::select(id, x, y, z) %>% rename(neighbor3=y,
                                                               neighbor3_id = id),
            by=c("x", "neighbor3", "z"))
atlas.d2 <-  atlas.d2 %>%
  left_join(atlas.d2 %>% dplyr::select(id, x, y, z) %>% rename(neighbor4=y,
                                                               neighbor4_id = id),
            by=c("x", "neighbor4", "z"))
atlas.d2 <-  atlas.d2 %>%
  left_join(atlas.d2 %>% dplyr::select(id, x, y, z) %>% rename(neighbor5=z,
                                                               neighbor5_id = id),
            by=c("x",  "y", "neighbor5"))
atlas.d2 <-  atlas.d2 %>%
  left_join(atlas.d2 %>% dplyr::select(id, x, y, z) %>% rename(neighbor6=z,
                                                               neighbor6_id = id),
            by=c("x",  "y", "neighbor6"))

##### Make a graph from there


##### Now I can make a adjacency matrix



library(Matrix) 
A = as_adjacency_matrix(G)[which(atlas.d2["value"]>0), which(atlas.d2["value"]>0)] 
A = A + Diagonal(nrow(A))
X_tilde = as.matrix(as.matrix(X[, 2:ncol(X)]) %*% A)/6
X_tilde = X_tilde %*% A/6
X_tilde = X_tilde %*% A/6  #### We smooth everything out...
##### We then want to add more constrasts .... Maybe sparsity is not the best bet here
##### We do need to say something about the graph behaves

# G = make_lattice(dim(t1))  %>%
#   set_vertex_attr("name", value = 1:p)
#### We can also 

scatterplot3d(atlas.d2 %>% dplyr::select(x,y,z), angle = 55)

atlas = readnii("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/atlases/Schaefer2018_800Parcels_7Networks_order_FSLMNI152_2mm.nii.gz")
atlas.d = atlas@.Data








files = list.files("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/wm_ntgtvsbl/")
ids_length = 6
c = c()
c_ids = c()

Y = read_csv("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/cognitive_data.csv")
Y[is.na(Y)] = 0
dataset <- c()
i=1
atlas["voxel"] = 1:nrow(atlas)
for(file in files){
  print(i)
  t1 = readnii(paste0("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/wm_ntgtvsbl/", file))
  test = t1@.Data
  test[is.na(test)] = 0
  test = data.frame(as.matrix(test, nrow=1))
  colnames(test) = "Response"
  test["voxel"] = 1:nrow(test)
  test["id"] = file
  test["id_num"] = i
  test = test %>% 
          left_join(atlas, by="voxel")
  test = test %>% 
    group_by(`1000`, `900`, `800`, `700`,
             `600`, `500`, `400`, `300`, `200`, `100`) %>%
    summarize_all(mean)
  dataset = rbind(dataset, test)
  i <- i+1
}


##### We should try to compute the correlation matrix of all of these patients.



##### Try diagnal thresholding to select a much smaller set of coefficients.
##### Of course, the arguments of Johnstone and Liu does not hold....
##### Maybe we can fix another way of selecting the different coordinates
dataset = dataset %>%
  left_join(atlas_agg %>% rename(voxel_id = id),
            by=c(  "1000","900", "800", "700" ,"600", "500" ,"400", "300", "200" ,"100"))

X = pivot_wider( pivot_longer(dataset %>% filter(GM > 0.1) %>% 
                                select(-c( "id", "x", "y", "z", "GM", "count")),
                               cols = c("id_num", "voxel_id"),
                              ),
                 id_cols = "id",
                 values_from = "Response" )

X = pivot_wider( dataset %>% filter(GM > 0.1) %>% 
                            select(-c( "id", "x", "y", "z", "GM", "count")),
            id_cols = c(  "1000","900", "800", "700" ,"600", "500" ,"400", "300", "200" ,"100", "voxel_id"),
            names_from = "id_num",
            values_from = "Response",
            names_repair =  "universal_quiet")
X  = as.matrix(X[, 12:ncol(X)])
C = (X) %*% t(X)/(nrow(X))
#### Note that the individual voxels, as such, do not have idenitcal variance... 
#### We have to project them unto the basis of interest

##### Use empirical Bayes to threshold the entries of the covariance matrix that are non negative
##### Or simply suse thresholding at level sigma * sqrt(2 log(p)/n)
p = ncol(C)
n = 269
sigma = median(diag(C))
tau = sigma *(1 + sqrt(2* log(p)/n))
s = which(diag(C)>tau)
ggplot(data = data.frame(x= diag(C)), aes(x=x)) +
  geom_histogram()
#


unique_groups = dataset %>%  
  dplyr::select(`1000`, `900`, `800`, `700`,
         `600`, `500`, `400`, `300`, `200`, `100`)
##### Now we make the graph
##### Tree:
atlas_agg 
g = make_empty_graph(n = 1000 )
sim = atlas %>% 
         group_by(200) %>%
         mutate(count_single200= n()) %>%
         ungroup()

sim2 = sim %>% 
  group_by(100) %>%
  mutate(count_single= n()) %>%
  ungroup()%>% 
  group_by(100, 200) %>%
  summarize(count_inter= n(),
            count_single=mean(count_single),
            count_single200=mean(count_single200)) %>%
  ungroup() %>%
  mutate(prop = count_inter/count_single200) %>%
  group_by(200) %>%
  top_n(1, prop)


atlas_agg = atlas %>% 
  dplyr::select(`1000`, `900`, `800`, `700`, `600`,
                `500`, `400`, `300`, `200`, `100`)  %>%
  group_by(`1000`, `900`, `800`, `700`, `600`,
           `500`, `400`, `300`, `200`, `100`) %>%
  summarise(count = n()) 
atlas_agg["id"] = 1:nrow(atlas_agg)

atlas_agg2 = atlas_agg %>%
  filter(`1000`>1)

atlas_agg2= atlas_agg2 %>%
  group_by(`1000`) %>%
  mutate(total_1000 = sum(count),
         distinct_paths = n())
         
atlas_agg2 = atlas_agg2 %>%
  mutate(percentage=count/ total_1000)

atlas_agg3 = pivot_longer(atlas_agg2 %>% 
                            dplyr::select(id, `1000`, `900`, 
                                          `800`, `700`, `600`,
                                                `500`, `400`, 
                                          `300`, `200`, `100`),
                          cols =-c("id"))


library(proxy)
jaccard_distance <- sparseMatrix(i = integer(0), j = integer(0), dims = c(nrow(atlas_agg), (nrow(atlas_agg))))
                                 
dist((atlas_agg4 %>% dplyr::select(-c(id)))[1:1000,], method = "jaccard")



D = rdist(atlas_agg4 , metric = "jaccard")
ggplot(atlas_agg2, aes(x= percentage)) + 
  geom_histogram()
ggplot(atlas_agg2, aes(x= distinct_paths)) + 
  geom_histogram()
#### Tackle the case of 100 first
####
count1 = atlas_agg %>% 
  group_by(`100`) %>%
  mutate(count_single_new_level= n()) %>%
  ungroup()

##### Fix some random order
W = c()
for (i in 1:99){
  ind = count1 %>% 
    filter(`100` == i)
  ind2 = count1 %>% 
    filter(`100` == i + 1)
  vector_a = rep(0, 1000)
  vector_a[ind$id] = 1/sqrt(nrow(ind))
  vector_a[ind2$id] = -1/sqrt(nrow(ind2))
  W = cbind(W, vector_a)
}
prev_t = "100"
ref = count1
for (t in sapply(seq(200, 1000, 100), as.character)) {
  #### find the cluster
  count1 = atlas_agg %>% 
        group_by(!! rlang::sym(t)) %>%
        mutate(count_single_new_level= n_distinct(id)) %>%
        ungroup()
      
  count2 = atlas_agg %>%   ### filter that mess
        group_by(!! rlang::sym(t)) %>%
        mutate(count_single_new_level= n()) %>%
        ungroup()%>% 
        group_by(!! rlang::sym(t), !! rlang::sym(prev_t)) %>%
        summarize(count_inter= n(),
                  count_single_new_level=mean(count_single_new_level)) %>%
        ungroup() %>%
        mutate(prop = count_inter/count_single_new_level) %>%
        group_by(!! rlang::sym(t)) %>%
        top_n(1, prop)
      
      #### group by the largest cluster
    for (i in 1:(as.numeric(prev_t))){
      clusters = count2 %>% filter(!! rlang::sym(prev_t) == i)
      if (nrow(clusters)>1){
        all_clusters = as.numeric(unlist(clusters %>% ungroup()  %>% dplyr::select(!! rlang::sym(t))))
        for (u in 1:(length(all_clusters)-1)){
          ind = atlas_agg %>% 
            filter(!! rlang::sym(t) == all_clusters[u],
            )
          ind = atlas_agg %>% 
            filter(!! rlang::sym(t) == all_clusters[u+1],
            )
          vector_a = rep(0, 1000)
          vector_a[ind$id] = 1/sqrt(nrow(ind))
          vector_a[ind2$id] = -1/sqrt(nrow(ind2))
          W = cbind(W, vector_a)
        }
      }
      
    }
    prev_t = t
}

#### I can then try to express everything in an orthonormal basis.
#### graph
#### 

library(paletteer)
beta = rep(0, ncol(W))
beta[sample(1:ncol(W), 100)] = 3

atlas_agg_coordinates = atlas %>%
  group_by(`1000`) %>%
  summarize(x= median(x),
            y=median(y),
            z=median(z)) %>%
  filter(`1000`>0)
Y = W %*% as.vector(beta)

atlas_agg_coordinates["Y"] = Y
library(paletteer)
nColor <- 5
colors <- as.character(paletteer_c("viridis::inferno", n = nColor))
# Transform the numeric variable in bins
rank <- as.factor( as.numeric( cut(abs(Y), nColor)))
unique(rank)

scatterplot3d(atlas_agg_coordinates$x,atlas_agg_coordinates$y,
              atlas_agg_coordinates$z, pch = 16, color=colors[rank])
library(rgl)
plot3d(heightweight, aes(x = ageYear, y = heightIn, color = weightLb), shade = 0.5, angle = 30)
### jaccard distance
sim3 = pivot_wider(sim2,
                   id_cols = "id",
                   names_from = c(name, value),
                   names_glue = "{name}_{value}",
                   values_from = "v",
                   values_fill=0)

library(vegan)
distance_matrix <- vegdist(as.matrix(sim3 %>% dplyr::select(-c(id))), method = "jaccard")
h_clusters <- hclust(distance_matrix, method="ward.D2")

sim = atlas %>% 
  dplyr::select(1000, 900, 800, 700, 600,
                500, 400, 300, 200, 100)  %>%
  group_by(1000, 900, 800, 700, 600,
           500, 400, 300, 200, 100) %>%
  summarise(count = n()) 
sim["id"] = 1:nrow(sim)
sim2 = pivot_longer(sim, cols=-c("id"))
sim2 = pivot_wider(sim2,
                   id_cols = "id",
                   names_from = "value",
                   values_fill=1)
sim  = sim %>%
  mutate(1000 = 1000/count)

##### Once we have his, we essentially hae something that looks like a tree 
###### Let's make it into a tree ######

  
test = atlas  %>% 
  dplyr::group_by(1000) %>% 
  mutate(count = n())  
test = test  %>% 
  dplyr::group_by(1000, 900) %>% 
  summarise(count = mean(count),
            nb = n())  %>%
  mutate(prop = nb/ count)


##### The graph is not so clean

edges = atlas %>%
  group_by(1000) %>%
  summarize_all(median)


edges = atlas %>%
  group_by(1000) %>%
  summarize_all(max)



##### Analyse the correlation structure between voxels in the ROI.

edges_roi = dataset %>%
  filter(`1000`==1, `900`==1,
         `800`==20, `700`==1, `600`==1,
         `400` == 4) 
  
  
cor_test= cor(pivot_wider(edges_roi, id_cols = c("id_num"),
                          names_from = "voxel" ))
  
  
  
sample_1000 = (atlas_agg %>% filter(count > 50))[, sample(1:nrow(atlas_agg %>% filter(count > 50)), 20)]
corr_all = c()
for(k in 1:5){
  ind = sample(1:200, 1)
  temp = dataset %>% filter(`200`==ind)
  cor_test= cor(pivot_wider(temp%>% dplyr::select(id_num, voxel, Response),
                            id_cols =  "id_num",
                            values_from = "Response",
                            names_from = "voxel" ))
  upper = cor_test[upper.tri(cor_test,diag = F)]
  corr_all = rbind(corr_all, data.frame("corr" = upper,
                                        "i"=1:length(upper),
                                        "region"=ind))
}


#####
corr_all$region =   as.factor(corr_all$region)
ggplot(corr_all, aes(x=corr, y = ..density..,fill=region)) + 
  geom_histogram(alpha = 0.5) + 
  xlab("Correlation") + 
  #ylab("Density") + 
  theme_bw()



##### Step 1: Approximate the Sigma_x
m <- 20
p <- 2
k <- 4
Q1 <- tril(kronecker(Matrix(seq(0.1,p,length=p*p),p,p),diag(m)))
Q2 <- cbind(Q1,Matrix(0,m*p,k))
Q3 <- rbind(Q2,cbind(Matrix(rnorm(k*m*p),k,m*p),Diagonal(k)))
V <- tcrossprod(Q3)
CH <- Cholesky(V)
CH@x


delta = 2
S = 1/n * t(X) %*% X
theta = 1/n * (t(X)%*%X)^2  - S^2
S[which(abs(S) < delta * sqrt(theta * log(p/n)))] = 0
######
#



# Create an empty graph
library(igraph)
atlas_agg4 = atlas_agg %>%
  filter(`400`>0) %>%
  group_by(`400`, `300`, `200`,`100`) %>%
  summarize(n400 = n(),
            n_voxels400 = sum(count))

p = nrow(atlas_agg4)
embs = atlas_agg4
embs["id"] = 1:p
embs = pivot_longer(embs %>% 
                            dplyr::select(-c("n_voxels400")),
                          cols =-c("id"))
embs["fill"] = 1
embs = pivot_wider(embs, 
                         id_cols = c("id"),
                         names_from = c("name", "value"),
                         values_fill=0)
jaccard_dist = proxy::dist(embs, method = "Jaccard")
W = exp(-as.matrix(jaccard_dist)/(2 * (5 * median(jaccard_dist))))
####
library(igraph)
graph <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)

V(G)
library(reticulate)
graph_to_networkx <- function(g) {
  use_python("/Users/cdonnat/anaconda3/bin/python")  # You may need to modify the path to your Python interpreter
  py_networkx <- import("networkx")
  py_g <- py_networkx$Graph()
  
  for (edge in as_edgelist(g)) {
    py_g$add_edge(edge[1], edge[2])
  }
  
  return(py_g)
}

nx_g <- graph_to_networkx(graph)

estimate_treewidth <- function(nx_g) {
  use_python("/usr/bin/python3")  # You may need to modify the path to your Python interpreter
  py_networkx <- import("networkx")
  
  treewidth, _ = py_networkx$approximation$treewidth$min_degree(nx_g)
  return(treewidth)
}

treewidth <- estimate_treewidth(nx_g)
cat("The estimated treewidth of the graph is:", treewidth)


MST <- mst(graph)

# heatmap(as.matrix(jaccard_dist))
G = make_empty_graph(n=p)
G <- add.edges(G, as.matrix(atlas.d2 %>% dplyr::select(id, neighbor1_id) %>% drop_na(neighbor1_id)))
G <- add.edges(G, as.matrix(atlas.d2 %>% dplyr::select(id, neighbor2_id) %>% drop_na(neighbor2_id)))
G <- add.edges(G, as.matrix(atlas.d2 %>% dplyr::select(id, neighbor3_id) %>% drop_na(neighbor3_id)))
G <- add.edges(G, as.matrix(atlas.d2 %>% dplyr::select(id, neighbor4_id) %>% drop_na(neighbor4_id)))
G <- add.edges(G, as.matrix(atlas.d2 %>% dplyr::select(id, neighbor5_id) %>% drop_na(neighbor5_id)))
G <- add.edges(G, as.matrix(atlas.d2 %>% dplyr::select(id, neighbor6_id) %>% drop_na(neighbor6_id)))


