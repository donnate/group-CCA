library(oro.nifti)
library(neurobase)
library(plotly)
library(tidyverse)

dir = 'Data/'
mask_tensor = readnii(paste0(dir, 'gm_mask020_bin.nii'))
loci_3d_coord = which(mask_tensor == 1, arr.ind = TRUE)
colnames(loci_3d_coord) = c("x", "y", "z")
write.csv(loci_3d_coord, file = paste0(dir, 'loci_3d_coord.csv'), row.names = F)

############### Create activation matrix ###############

files = list.files(paste0(dir, 'activations_raw'))
activations = c()
for(file in files){
  activation = readnii(paste0(dir, "activations_raw/", file))
  cat(substr(file, 1, 7), " dim:", dim(activation), " missing:", sum(is.na(activation)), "\n")
  activations = rbind(activations, activation[loci_3d_coord])
}
rownames(activations) = substr(files, 1, 7)
write.csv(t(activations), file = paste0(dir, 'activations.csv'), row.names = F)


############### Smooth imputation of brain activations ###############

loci_index_tensor = mask_tensor
loci_index_tensor[loci_3d_coord] = 1:nrow(loci_3d_coord)
#test: loci_index[loci_3d_coord[1:10,]]

unit_ball_3d_coord = function(r){
  cube = expand.grid((-r):r, (-r):r, (-r):r)
  colnames(cube) = c("x", "y", "z") 
  dist = sqrt(rowSums(cube^2))
  ball = cube[dist <= r & dist != 0,]
  return(ball)
}

# # sanity checks: plot unit ball
# data.frame(unit_ball_3d_coord(3)) %>%
# plot_ly(x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers")

neighbors_index = function(index, r){
  ball = scale(unit_ball_3d_coord(r), center = -loci_3d_coord[index,], scale = F)
  inside = (rowSums(ball < 1) == 0)
  neigh_index = loci_index_tensor[ball[inside, ]]
  return(neigh_index[neigh_index != 0])
}

# #sanity checks: plot the ball + center
# data.frame(rbind(loci_3d_coord[20000,],
#       loci_3d_coord[neighbors_index(20000,2),])) %>%
# plot_ly(x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers")

activations_imp = activations
for(file in substr(files, 1, 7)){
  nas = which(is.na(activations[file,]))
  for(index in nas){
    r = 1
    impute = NA
    while(is.na(impute)){
      impute = mean(activations[file, neighbors_index(index, r)], na.rm = TRUE)
      r = r + 1
    } 
    activations_imp[file, index] = impute
  }
  cat(substr(file, 1, 7), " missing:", sum(is.na(activations_imp)), "\n")
}
rownames(activations_imp) = substr(files, 1, 7)
write.csv(t(activations_imp), file = paste0(dir, 'activations_imputed.csv'), row.names = F)

# #sanity checks: compare signal before and after imputation
# id = 10
# nas = which(is.na(activations[id,]))
# activations_imp[id, nas[1]]
# mean(activations[id, neighbors_index(nas[1], 5)], na.rm = T)
# data.frame(loci_3d_coord[nas,], color = activations_imp[id, nas]) %>%
#   plot_ly(x = ~x, y = ~y, z = ~z, color =~color, size = 2, type = "scatter3d", mode = "markers")


############### Reduce brain data dimensionality ###############

loci_center = round(colMeans(loci_3d_coord))
loci_3d_coord_centered = scale(loci_3d_coord, center = loci_center, scale = F)

# data.frame(loci_3d_coord_centered) %>%
#   plot_ly(x = ~x, y = ~y, z = ~z, size = 1, type = "scatter3d", mode = "markers")

nloci = 10
loci_3d_coord_centered_reduced = round(loci_3d_coord_centered/nloci)

# data.frame(loci_3d_coord_centered_reduced) %>%
#   plot_ly(x = ~x, y = ~y, z = ~z, size = 1, type = "scatter3d", mode = "markers")

activations_reduced = data.frame(t(activations_imp), loci_3d_coord_centered_reduced) %>% 
  group_by(x, y, z) %>% 
  summarise_at(vars(starts_with("CONN")), ~mean(.)) %>%
  arrange(x, y, z) %>% 
  ungroup()
nrow(activations_reduced)

activations_reduced %>%
  plot_ly(x = ~x, y = ~y, z = ~z, color = ~CONN070, size = 1, type = "scatter3d", mode = "markers")

write.csv(activations_reduced %>% dplyr::select(x, y, z), paste0(dir, "loci_3d_coord_reduced", nloci,".csv"), row.names = F)
write.csv(activations_reduced %>% select_at(vars(starts_with("CONN"))), paste0(dir, "activations_reduced", nloci, ".csv"), row.names = F)


############### Behavior scores ######################

bscores = read.csv(paste0(dir, "behavior_scores_raw.csv"), header = T) %>% dplyr::select(-redcap_event_name)  %>% 
  column_to_rownames(var = "record_id")
hist(rowMeans(is.na(bscores)))
bscores = bscores[rowMeans(is.na(bscores)) <= 0.2,]

write.csv(t(bscores), paste0(dir, "behavior_scores.csv"), row.names = T)

