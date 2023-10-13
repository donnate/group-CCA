library(tidyverse)

file_list <- list.files(path = "~/Documents/group-CCA/experiments/sparse_CCA/results/", 
                               pattern = "extended_results_exp_sparse_cca_new_exp_*", full.names = TRUE)
# Read and combine CSV files into a single data frame
results <- file_list %>%
  map_dfr(~ read_csv(.x) %>% mutate(filename = .x)) %>%
  mutate("selection" = ifelse(str_detect(filename, "correlation"), "correlation", "prediction"))


##### Analyse this data
res = results %>% 
  mutate(FDR=1 - nb_real_discoveries/max(1,nb_discoveries)) %>%
  group_by(method, selection,criterion, n, nnz, p1, p2, overlapping_amount) %>% 
  summarise_if(is.numeric, mean) %>%
  arrange(n, nnz, p1, p2, distance_tot)
t = as.numeric(sapply(res$method, function(x){ifelse(str_count(x, "folds" )>0, 1, 0)}))
t = t  + as.numeric(sapply(res$method, function(x){ifelse(str_count(x, "CV" )>0, 2, 0)}))
res["shape"] = as.factor(t)


legend_order <- c( "SSVD-theory" ,  "SSVD-method",
                   "adaptive_lasso_Fantope", "adaptive_lasso_Selection",
                  "TG" ,"lasso_Fantope", "lasso_Selection",   "Fantope" ,
                  "thresholded-lasso",
                  "Canonical Ridge-Author",
                  "FIT_SAR_BIC", "FIT_SAR_CV", "SCCA_Parkhomenko",
                  "Witten_Perm", "Witten.CV", "Waaijenborg-Author",
                  "Waaijenborg-CV", "Selection",  "Oracle" )

my_colors <- c( "gray", "black",
  "dodgerblue", "cyan",  "red","orange","orange4","navy", "grey", "yellow",  
              "green", "green3", "navajowhite3", "plum2", "plum3",
              "cadetblue1", "lightskyblue", "brown", "beige","whitesmoke", "white")
              


ggplot(res %>% filter(method %in% c( "SSVD-theory" ,  "SSVD-method",
                                     "adaptive_lasso_Fantope", "adaptive_lasso_Selection",
                                     "TG" ,  "Fantope" ), selection=="prediction"))+
  geom_line(aes(x=p1/n, y=distance_tot, colour=method), linewidth=1.)+
  geom_point(aes(x=p1/n, y=distance_tot, colour=method, size=2.2)+
  #geom_jitter(aes(x=p1/n, y=distance_tot, colour=method, shape=shape), size=2.2)+
  geom_line(data=res %>% group_by(n, nnz, p1, p2) %>% summarise(b=mean(zero_benchmark)),
             aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~sparsity, scales="free") +
  theme_bw()


ggplot(res %>% filter(method %in% c("adaptive_lasso_Fantope", "adaptive_lasso_Selection",
                                    "TG" ,  "Fantope"  )))+
  geom_line(aes(x=p1/n, y=FNR, colour=method, linetype=shape), linewidth=1.)+
  geom_point(aes(x=p1/n, y=FNR, colour=method, shape=shape), size=2.2)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~sparsity, scales="free") +
  theme_bw()



ggplot(res  %>% filter(!method %in% c("lasso_with_folds",
                                      "adaptive_lasso_with_folds")))+
  geom_line(aes(x=p1/n, y=FNR, colour=method), linewidth=1.)+
  geom_point(aes(x=p1/n, y=FNR, colour=method), size=2.2)+
  # geom_line(data=res  %>% filter(nnz<30)%>% 
  #             mutate(fpr_benchmark = (n-nnz)/n) %>%
  #             group_by(n, nnz, p1, p2) %>% 
  #             summarise(b=mean(fpr_benchmark)),
  #           aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~sparsity, scales="free") +
  theme_bw()


ggplot(res %>% filter(!method %in% c( "lasso_with_folds",
                                     "adaptive_lasso_with_folds")))+
  geom_line(aes(x=p1/n, y=FNR, colour=method), linewidth=1.)+
  geom_point(aes(x=p1/n, y=FNR, colour=method), size=2.2)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~sparsity, scales="free") +
  theme_bw()


ggplot(res   %>% filter(!method %in% c( "lasso_with_folds","Fantope",
                                     "adaptive_lasso_with_folds",
                                     "Waaijenborg-Author",
                                     "Waaijenborg-CV")))+
  geom_line(aes(x=p1/n, y=(nb_discoveries), colour=method), linewidth=1.)+
  geom_point(aes(x=p1/n, y=nb_discoveries, colour=method), size=2.2)+
  #geom_jitter(width = 0.2, height = 0) +
  geom_line(data=res %>% 
               mutate(fnr_benchmark = nnz) %>%
               group_by(n, sparsity, p1, p2) %>% 
               summarise(b=mean(ceiling(sparsity * p1))),
              aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  geom_line(data=res %>% 
              mutate(fnr_benchmark = nnz) %>%
              group_by(n, sparsity, p1, p2) %>% 
              summarise(b=2 * mean(ceiling(p1))),
            aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~sparsity, scales="free") +
  #scale_y_log10() + 
  theme_bw()


ggplot(res %>% filter(!method %in% c(  "lasso_with_folds","Fantope",
                                     "adaptive_lasso_with_folds")),
       aes(x=FPR, y=FNR, colour=method, size=p1/n))+
  geom_point()+
  geom_jitter(width = 0.02, height=0.02)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(sparsity~n, scales="free") +
  theme_bw()


ggplot( res%>% filter(!method %in% c("TG_normalized", "lasso_with_folds",
                                  "adaptive_lasso_with_folds")))+
  geom_line(aes(x=p1/n, y=TPR, colour=method))+
  geom_point(aes(x=p1/n, y=TPR, colour=method))+
  geom_line(data=res %>% 
              mutate(tpr_benchmark = (nnz)/n) %>%
              group_by(n, nnz, p1, p2) %>% 
              summarise(b=mean(tpr_benchmark)),
            aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_wrap(n~nnz, scales="free") +
  theme_bw()





############
file_list1 <- list.files(path = "~/Documents/group-CCA/experiments/sparse_CCA/results/results/", pattern = "results_exp_localized_cca_localized_27466*", full.names = TRUE)
# Read and combine CSV files into a single data frame
results_loc <- file_list1 %>%
  map_dfr(~ read_csv(.x) %>% mutate(filename = .x))

##### Try and analyse this data




res_loc = results_loc %>% 
  group_by(method, power, example_type, n, p1, p2) %>% 
  summarise_if(is.numeric, mean) %>%
  arrange(n, power, example_type, nnz, p1, p2, distance)

t = as.numeric(sapply(res_loc$method, function(x){ifelse(str_count(x, "folds" )>0, 1, 0)}))
t = t  + as.numeric(sapply(res_loc$method, function(x){ifelse(str_count(x, "CV" )>0, 2, 0)}))
res_loc["shape"] = as.factor(t)


legend_order <- c(   "adaptive_regularised_lasso" ,    "adaptive_regularised_lasso_with_folds",
                     "adaptive_lasso", "TG", "regularised_lasso",   "regularised_lasso_with_folds",
                     "lasso", "Canonical Ridge-Author",
                  "FIT_SAR_BIC", "FIT_SAR_CV", "SCCA_Parkhomenko",
                  "Witten_Perm", "Witten.CV", "Waaijenborg-Author",
                  "Waaijenborg-CV")
test = res_loc %>% filter(power==0.1)
my_colors <- c(
 "navy", "dodgerblue4", "dodgerblue",  "red", "brown3", "salmon",  "orange",  "yellow", 
              "green", "green3", "navajowhite3", "plum2", "plum3",
              "cadetblue1", "cyan", "whitesmoke")
              
ggplot(res_loc  %>%filter(n==100) %>% filter(!method %in% c("TG_normalized", "lasso_with_folds",
                                                        "adaptive_lasso_with_folds")))+
  geom_line(aes(x=p1/n, y=distance, colour=method, linetype=shape), linewidth=1.)+
  geom_point(aes(x=p1/n, y=distance, colour=method, shape=shape), size=2.2)+
  geom_line(data=res_loc %>%  group_by(n, nnz, p1, p2) %>% summarise(b=mean(zero_benchmark)),
            aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(power~ example_type +n, scales="free") +
  theme_bw()


 ggplot(res_loc  %>% filter(!method %in% c("TG_normalized", "lasso_with_folds",
                                          "adaptive_lasso_with_folds")))+
  geom_line(aes(x=p1/n, y=FNR, colour=method, linetype=shape), linewidth=1.)+
  geom_point(aes(x=p1/n, y=FNR, colour=method, shape=shape), size=2.2)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_wrap(n~power, scales="free") +
  theme_bw()


ggplot(res_loc  %>% filter(!method %in% c("TG_normalized", "lasso_with_folds",
                                          "adaptive_lasso_with_folds")))+
  geom_line(aes(x=p1/n, y=nb_discoveries, colour=method, linetype=shape), linewidth=1.)+
  geom_point(aes(x=p1/n, y=nb_discoveries, colour=method, shape=shape), size=2.2)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_wrap(n~power, scales="free") +
  theme_bw()

ggplot(res_loc%>% filter(!method %in% c("TG_normalized", "lasso_with_folds",
                                                        "adaptive_lasso_with_folds")),
       aes(x=FPR, y=FNR, colour=method, size=p1/n))+
  geom_point()+
  geom_jitter(width = 0.02, height=0.02)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_wrap(power~n, scales="free") +
  theme_bw()









#################################
#################################


library(tidyverse)



file_list <- list.files(path = "~/Documents/group-CCA/experiments/sparse_CCA/results/results/", pattern = "extended_results_exp_sparse_cca_*", full.names = TRUE)
# Read and combine CSV files into a single data frame
results <- file_list %>%
  map_dfr(~ read_csv(.x) %>% mutate(filename = .x))


##### Try and analyse this data




res = results %>% 
  group_by(method, n, nnz, p1, p2) %>% 
  summarise_if(is.numeric, mean) %>%
  arrange(n, nnz, p1, p2, distance)
t = as.numeric(sapply(res$method, function(x){ifelse(str_count(x, "folds" )>0, 1, 0)}))
t = t  + as.numeric(sapply(res$method, function(x){ifelse(str_count(x, "CV" )>0, 2, 0)}))
res["shape"] = as.factor(t)


legend_order <- c("adaptive_lasso", "TG","lasso",  "thresholded-lasso", "Canonical Ridge-Author",
                  "FIT_SAR_BIC", "FIT_SAR_CV", "SCCA_Parkhomenko",
                  "Witten_Perm", "Witten.CV", "Waaijenborg-Author",
                  "Waaijenborg-CV")

my_colors <- c(
  "dodgerblue",  "red","orange",  "black", "yellow", 
              "green", "green3", "navajowhite3", "plum2", "plum3",
              "cadetblue1", "lightskyblue", "whitesmoke")
              


ggplot(res %>% filter(!method %in% c("lasso_with_folds",
                                     "adaptive_lasso_with_folds")))+
  geom_line(aes(x=p1/n, y=distance, colour=method, linetype=shape), linewidth=1.)+
  geom_point(aes(x=p1/n, y=distance, colour=method, shape=shape), size=2.2)+
  geom_line(data=res %>% group_by(n, nnz, p1, p2) %>% summarise(b=mean(zero_benchmark)),
            aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~sparsity, scales="free") +
  theme_bw()


ggplot(res  %>% filter(!method %in% c("lasso_with_folds",
                                      "adaptive_lasso_with_folds")))+
  geom_line(aes(x=p1/n, y=FPR, colour=method), linewidth=1.)+
  geom_point(aes(x=p1/n, y=FPR, colour=method), size=2.2)+
  # geom_line(data=res  %>% filter(nnz<30)%>% 
  #             mutate(fpr_benchmark = (n-nnz)/n) %>%
  #             group_by(n, nnz, p1, p2) %>% 
  #             summarise(b=mean(fpr_benchmark)),
  #           aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~sparsity, scales="free") +
  theme_bw()


ggplot(res %>% filter(!method %in% c( "lasso_with_folds",
                                      "adaptive_lasso_with_folds")))+
  geom_line(aes(x=p1/n, y=FNR, colour=method), linewidth=1.)+
  geom_point(aes(x=p1/n, y=FNR, colour=method), size=2.2)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~sparsity, scales="free") +
  theme_bw()


ggplot(res   %>% filter(!method %in% c( "lasso_with_folds",
                                        "adaptive_lasso_with_folds",
                                        "Waaijenborg-Author",
                                        "Waaijenborg-CV")))+
  geom_line(aes(x=p1/n, y=(nb_discoveries), colour=method), linewidth=1.)+
  geom_point(aes(x=p1/n, y=nb_discoveries, colour=method), size=2.2)+
  #geom_jitter(width = 0.2, height = 0) +
  geom_line(data=res %>% 
              mutate(fnr_benchmark = nnz) %>%
              group_by(n, sparsity, p1, p2) %>% 
              summarise(b=mean(ceiling(sparsity * p1))),
            aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  geom_line(data=res %>% 
              mutate(fnr_benchmark = nnz) %>%
              group_by(n, sparsity, p1, p2) %>% 
              summarise(b=2 * mean(ceiling(p1))),
            aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(n~sparsity, scales="free") +
  scale_y_log10() + 
  theme_bw()


ggplot(res %>% filter(!method %in% c(  "lasso_with_folds",
                                       "adaptive_lasso_with_folds")),
       aes(x=FPR, y=FNR, colour=method, size=p1/n))+
  geom_point()+
  geom_jitter(width = 0.02, height=0.02)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_grid(sparsity~n, scales="free") +
  theme_bw()


ggplot( res%>% filter(!method %in% c("TG_normalized", "lasso_with_folds",
                                     "adaptive_lasso_with_folds")))+
  geom_line(aes(x=p1/n, y=TPR, colour=method))+
  geom_point(aes(x=p1/n, y=TPR, colour=method))+
  geom_line(data=res %>% 
              mutate(tpr_benchmark = (nnz)/n) %>%
              group_by(n, nnz, p1, p2) %>% 
              summarise(b=mean(tpr_benchmark)),
            aes(y=b, x=p1/n, colour="Zero Benchmark"), colour="black", linewidth=0.5)+
  scale_color_manual(values = my_colors, breaks = legend_order) +
  facet_wrap(n~nnz, scales="free") +
  theme_bw()










