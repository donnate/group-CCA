name_experiment = paste0(namefile, 'environment')

paste0(mysavedir, name_experiment, id_exp,'_all_results.csv')

res = c()
for (u in list.files(path = mysavedir , pattern = paste0("(.*?)", name_experiment, '(.*?)', '_all_results.csv'))){
  res = rbind(res, 
              read_csv(paste0(mysavedir, u)))
}
