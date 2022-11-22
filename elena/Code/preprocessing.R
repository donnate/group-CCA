
#set up the folder that contains raw data: Rest or Emotion
task = 'Emotion'

#store all values above diagonal of connectome matrices in matrix c
dir = '../Data/HCP Healthy Young Adult/'
files = list.files(paste0(dir, 'Connmats/', task))
ids_length = 6
c = c()
c_ids = c()

for(file in files){
  data = read.csv(paste0(dir, 'Connmats/', task, '/', file), header = FALSE)
  c_ids = append(c_ids, substr(file, 1, ids_length))
  colnames(data) = NULL
  c = rbind(c, data[upper.tri(data, diag = FALSE)])
}
n = ncol(data)
c_names = matrix(paste('c', rep(1:n, n), rep(1:n, rep(n, n)), sep = '_'), n, n)
c_names = c_names[upper.tri(c_names, diag = FALSE)]
rownames(c) = c_ids
colnames(c) = c_names
write.csv(c, file = paste0(dir, 'activation_', tolower(task), '.csv'), row.names = TRUE)


#########################################

dir = '../Data/HCP Disordered Emotional States/'
X =  read.csv(paste0(dir, 'activation.csv'), header = TRUE, row.names = 1)
Y = read.csv(paste0(dir, 'behavior.csv'), header = TRUE, row.names = 1) %>% 
  select(-demo_age, -bio_sex)
groups = unlist(read.csv(paste0(dir, 'activation_groups.csv'), header = FALSE))
X_clean = X[, groups != 0]
groups_clean = groups[groups!=0]

X_mean = data.frame(t(X_clean)) %>%
  mutate(group = groups_clean) %>%
  group_by(group) %>%
  summarise_all(mean) %>%
  select(-group) %>%
  t()
rownames(X_mean) = gsub("[.]", "-", rownames(X_mean))
write.table(X_mean, file = paste0(dir, 'activation_avg.csv'), sep = ",", col.names = FALSE, row.names = TRUE)
