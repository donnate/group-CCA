library(tidyverse)
library(ggplot2)
library(pracma)
setwd("~/Documents/group-CCA/elena/")
file_list <- list.files(path = "~/Documents/group-CCA/elena/missing/", 
                        pattern = "simulation-RRR-results99*", full.names = TRUE)

results <- file_list %>%
  map_dfr(~ read_csv(.x) %>% mutate(filename = .x))


results["lambda"] = sapply(results$method, function(x){ifelse(is.null(strfind(x, "Alt-"))==FALSE,
                                                              as.numeric(strsplit(x, "Alt-")[[1]][2]),
                                                              NA)})
results["method_type"] = sapply(results$method, function(x){ifelse(is.null(strfind(x, "Alt-"))==FALSE,
                                                                   "Alt",
                                                                   x)})
summ = results %>% group_by(n, p1, p2, r, r_pca,
                            #nnzeros, 
                            overlapping_amount, noise, 
                            method_type, 
                            lambda,
                            method,
                            #theta_strength,
                            prop_missing) %>% 
  summarize_if(is.numeric, mean) %>% 
  ungroup() %>%
  arrange(n, p1, p2, r, noise, prop_missing, distance_U)

colnames(summ)
unique(summ$nnzeros)
unique(summ$n)
unique(summ$lambda)

ggplot(results %>% filter(n == 1000, method_type=="Alt"),
       aes(x=lambda, y =  1/p1 * distance_tot, 
           colour =method) )+
  geom_point(alpha=0.1)+
  geom_smooth()+
  scale_y_log10()+
  scale_x_log10() 

ggplot(summ %>% filter(n == 1000, 
                       r_pca == 0, r==2,
                       overlapping_amount == 0,
                       lambda %in% c(0.1, NA))) +
  geom_line(aes(x=prop_missing, 
                y = distance_tot, 
                colour =method))+
  facet_grid(theta_strength~ p1, scales = "free")

legend_order <- c("RRR", "Alt", "CCA-mean", "CCA-median") 
my_colors <- c( "red",
                "dodgerblue", "gray", "black")

labels_n = c()

unique(summ$n)
unique(summ$p1)
unique(summ$r_pca)
unique(summ$r)
unique(summ$overlapping_amount)
ggplot(summ %>% filter(n == 500,
                       r==2,
                       r_pca > 0,
                       method %in% legend_order)) +
  geom_line(aes(x=prop_missing, 
                y = distance_tot, 
                colour =method), linewidth=1.2)+
  scale_color_manual(values =my_colors, 
                     breaks =legend_order) +
  scale_y_log10() + 
  facet_grid(p1~n, scales = "free")
