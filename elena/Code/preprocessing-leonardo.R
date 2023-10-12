library(tidyverse)
library(RNifti)
library(oro.nifti)
library(neurobase)
#store all values above diagonal of connectome matrices in matrix c
dir = '../fMRI-data/data/'
files = list.files("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/wm_ntgtvsbl/")
ids_length = 6
c = c()
c_ids = c()

Y = read_csv("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/cognitive_data.csv")
Y[is.na(Y)] = 0
Y = Y %>%
  dplyr::filter(record_id %in% unique(Y$id))
atlas800 = readnii("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/atlases/Schaefer2018_800Parcels_7Networks_order_FSLMNI152_2mm.nii.gz")
atlas.d = atlas@.Data
dataset <- c()
for(file in files){
  print(file)
  t1 = readnii(paste0("/Users/cdonnat/Documents/group-CCA/fMRI-data/data/wm_ntgtvsbl/", file))
  test = t1@.Data
  test[is.na(test)] = 0
  test = data.frame(as.matrix(test, nrow=1))
  colnames(test) = "Response"
  test["voxel"] = 1:nrow(test)
  test["id"] = file
  test["Schaffer_group"] = as.vector(atlas.d)
  dataset = rbind(dataset, test)

}

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
#### Let's try our CCA approach?
#### 



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
