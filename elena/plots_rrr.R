library(tidyverse)
library(ggplot2)
library(pracma)
library(tidyverse)
theme_set(theme_bw(base_size = 14))
setwd("~/Documents/group-CCA/elena/")
file_list <- list.files(path = "~/Documents/group-CCA/elena/missing/results/", 
                        pattern = "new_RRR_efficient_results10*", full.names = TRUE)

results <- file_list %>%
  map_dfr(~ read_csv(.x) %>% mutate(filename = .x))



results["lambda"] = sapply(results$method, function(x){ifelse(is.null(strfind(x, "Alt-"))==FALSE,
                                                              as.numeric(strsplit(x, "Alt-")[[1]][2]),
                                                              NA)})
results["method_type"] = sapply(results$method, function(x){ifelse(is.null(strfind(x, "Alt-"))==FALSE,
                                                                   "Alt",
                                                                   x)})
# summ = results %>% group_by(n, p1, p2, r, r_pca,
#                             nnzeros, 
#                            overlapping_amount, noise, 
#                            #method_type, 
#                            lambda_opt,
#                            method,
#                            theta_strength,
#                            prop_missing) %>% 
#   summarize_if(is.numeric, median) %>% 
#   ungroup() 

summ = results %>% group_by(n, p1, p2, r, r_pca,
                            nnzeros, 
                            overlapping_amount, noise, 
                            #lambda_opt,
                            method,
                            theta_strength,
                            normalize_diagonal,
                            prop_missing) %>% 
  summarise(distance_tot = median(distance_tot),
            distance_tot_mean = mean(distance_tot),
            distance_U = median(distance_U),
            distance_V = median(distance_V),
            distance_tot_q50 = quantile(distance_tot, 0.5),
            distance_tot_q25 = quantile(distance_tot, 0.75),
            distance_tot_q75 = quantile(distance_tot, 0.25),
            prediction_tot_mean= mean(prediction_tot),
            prediction_tot_q50 = quantile(prediction_tot, 0.5),
            prediction_tot_q25 = quantile(prediction_tot, 0.75),
            prediction_tot_q75 = quantile(prediction_tot, 0.25),
            distance_U_q50 = quantile(distance_U, 0.5),
            distance_U_q25 = quantile(distance_U, 0.75),
            distance_U_q75 = quantile(distance_U, 0.25),
            prediction_U_mean= mean(distance_U),
            prediction_U_q50 = quantile(distance_U, 0.5),
            prediction_U_q25 = quantile(distance_U, 0.75),
            prediction_U_q75 = quantile(distance_U, 0.25),
            distance_V_q50 = quantile(distance_V, 0.5),
            distance_V_q25 = quantile(distance_V, 0.75),
            distance_V_q75 = quantile(distance_V, 0.25),
            prediction_V_mean= mean(distance_V),
            prediction_V_q50 = quantile(distance_V, 0.5),
            prediction_V_q25 = quantile(distance_V, 0.75),
            prediction_V_q75 = quantile(distance_V, 0.25),
            time_med = quantile(time, 0.5),
            time_mean = mean(time)
            
  ) %>%
  ungroup() 
colnames(results)
colnames(summ)
unique(summ$nnzeros)

unique(summ$n)


unique(summ$method)
legend_order <- c("Oracle",  "FIT_SAR_CV", 
                  "FIT_SAR_BIC", "Witten_Perm", "Witten.CV",
                  "SCCA_Parkhomenko", "Waaijenborg-CV", "Waaijenborg-Author",
                  "Canonical Ridge-Author", "CCA-mean",
                  "RRR" ,   "Alt-opt", "Gradient-descent",
                  "init-alternating")
my_colors <- c( "black", "chartreuse2", "chartreuse4",
                "orchid1", "orchid3", "indianred",
                "burlywood2", "burlywood4",
                "cyan", "gray", "red",
                "dodgerblue", "orange", "yellow")

labels_n <-   c("Oracle",  "FIT_SAR with CV (Wilms et al)", 
                "FIT_SAR with BIC (Wilms et al)",  
                "Witten et al (with Permutation Test)", "Witten et al.(with CV)",
                "SCCA (Parkhomenko et al)", "SCCA with CV (Waaijenborg et al)", 
                "SCCA with BIC (Waaijenborg et al)",
                "CCA with Ridge Penalty", "CCA-mean",
                "CCA as Reduced Rank Regression" , 
                "Alternating Regression (this paper)", 
                "Gradient Descent (this paper)",
                "Initialization (this paper)")

legend_order <- c("Oracle",  "FIT_SAR_CV", 
                  "FIT_SAR_BIC", "Witten_Perm", "Witten.CV",
                  "SCCA_Parkhomenko", "Waaijenborg-CV", "Waaijenborg-Author",
                  "Canonical Ridge-Author",
                  #"RRR-0.5" ,"RRR-7.5","RRR-10","RRR-12.5",  "RRR-20",   
                  "RRR-opt")
my_colors <- c( "black", "chartreuse2", "chartreuse4",
                "orchid1", "orchid3", "indianred",
                "burlywood2", "burlywood4",
                 "yellow",
                #"lightblue", "lightblue3","cyan", "dodgerblue", "dodgerblue4", 
                "navyblue")

labels_n <-   legend_order

unique(summ$r_pca)
unique(summ$p1)
unique(summ$n)
unique(summ$p2)
unique(summ$nnzeros)
unique(summ$overlapping_amount)
summ %>% filter( r_pca == 5, r==2,  p2==10, p1>50)

ggplot(summ %>% filter( r_pca == 0, r==2,  p2==20,# n<2000,
                        nnzeros==20, overlapping_amount == 0,
                        method %in% legend_order
                        ),
       aes(x=p1, 
           y = distance_U_q50, 
           colour =method)) +
  geom_point()+
  geom_line()+
  scale_color_manual(values = my_colors, breaks = legend_order,
                    labels = labels_n) +
  facet_grid(theta_strength~ n, scales = "free",labeller = as_labeller(c(`20` = "p = 20",
                                                                          `80` = "p = 80",
                                                                          `100` = "n = 100",
                                                                          `200` = "n = 200",
                                                                          `300` = "n = 300",
                                                                         `500` = "n = 500",
                                                                          `high` = "High",
                                                                         `1000` = "n = 1000",
                                                                         `2000` = "n = 2000",
                                                                         `10000` = "n = 10000",
                                                                          `medium` = "Medium",
                                                                          `low` = "Low"
                                                                          
  ))) +
  xlab("p") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  #scale_y_log10()+
  #scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


unique(summ$p1)
ggplot(summ %>% filter( r_pca == 0, r==2,  p2==20,
                        #n<2000,
                        #n>100,
                        nnzeros==20, overlapping_amount == 0,
                        method %in% legend_order
),
aes(x=n, 
    y = distance_U, 
    colour =method)) +
  geom_point()+
  geom_line()+
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  facet_grid(theta_strength~ p1, scales = "free",labeller = as_labeller(c(`20` = "p = 20",
                                                                          `50` = "p = 50",
                                                                         `80` = "p = 80",
                                                                         `100` = "p = 100",
                                                                         `200` = "p = 200",
                                                                         `300` = "n = 300",
                                                                         `500` = "p = 500",
                                                                         `high` = "High",
                                                                         `1000` = "n = 1000",
                                                                         `2000` = "n = 2000",
                                                                         `10000` = "n = 10000",
                                                                         `medium` = "Medium",
                                                                         `low` = "Low"
                                                                         
  ))) +
  xlab("n") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  #scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(summ %>% filter( 
                        r_pca == 0, r==2,
                        nnzeros == 5,
                        overlapping_amount == 1),
       aes(x=p, 
           y = distance_tot_q50, 
           colour =method)) +
  geom_point()+
  geom_line()+
  geom_point(aes(y=distance_tot_q25))+    facet_grid(theta_strength~ p1, scales = "free",
                                                     labeller = as_labeller(c(`20` = "p = 20",
                                                                              `80` = "p = 80",
                                                                              `100` = "p = 100",
                                                                              `200` = "p = 200",
                                                                              `300` = "p = 300",
                                                                              `high` = "High",
                                                                              `medium` = "Medium",
                                                                              `low` = "Low"
                                                                              
                                                     ))) +
  xlab("n (Number of Samples)") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  scale_y_log10()+
  scale_x_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

unique(summ$p1)
unique(summ$r_pca)
unique(summ$nnzeros)
ggplot(summ %>% filter( 
                        r_pca == 5, r==2,
                        overlapping_amount == 1,
                        nnzeros==20
                        ),
       aes(x=n, 
           y = distance_tot, 
           colour =method)) +
  geom_line()+
  geom_point() + 
  scale_x_log10() + 
  #scale_color_manual(values = my_colors, breaks = legend_order,
  #                   labels = labels_n) +
  scale_y_log10() + 
  #geom_errorbar(aes(ymin = distance_tot_q25, 
  #                  ymax = distance_tot_q75), 
  #              width = 0.1, alpha=0.7,
  #              linewidth=1.1, position = position_dodge(0.05))+
  facet_grid(theta_strength~ n, scales = "free",
             labeller = as_labeller(c(`20` = "p = 20",
                                      `80` = "p = 80",
                                      `100` = "p = 100",
                                      `high` = "High",
                                      `medium` = "Medium",
                                      `low` = "Low"
                                      
             ))) +
  xlab("p (Number of Samples)") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



labels_n <-   c("Oracle",  "FIT_SAR with CV (Wilms et al)", 
                "FIT_SAR with BIC (Wilms et al)",  
                "Witten et al (with Permutation Test)", "Witten et al.(with CV)",
                "SCCA (Parkhomenko et al)", "SCCA with CV (Waaijenborg et al)", 
                "SCCA with BIC (Waaijenborg et al)",
                "CCA with Ridge Penalty", "CCA-mean",
                "CCA as Reduced Rank Regression" , 
                "Alternating Regression (this paper)", 
                "Gradient Descent (this paper)",
                "Initialization (this paper)")

my_colors <- c( "black", "chartreuse2", "chartreuse4",
                "orchid1", "orchid3", "indianred",
                "burlywood2", "burlywood4",
                "cyan", "gray", "red",
                "dodgerblue", "orange", "yellow")


legend_order <- c("Oracle",  "FIT_SAR_CV", 
                  "FIT_SAR_BIC", 
                  "CCA-mean",  "Alt-0.1",
                  "Alt-opt", "Gradient-descent",
                  "init-alternating")

labels_n <-   c("Oracle",  "FIT_SAR with CV (Wilms et al)", 
                "FIT_SAR with BIC (Wilms et al)", "CCA-mean",
                "Alt-0.1",
                "Alternating Regression (this paper)", 
                "Gradient Descent (this paper)",
                "Initialization (this paper)")


ggplot(summ %>% filter( nnzeros == 5, method %in% legend_order,
                        r_pca == 5, r==2,
                        overlapping_amount == 0),
       aes(x=p1, 
           y = distance_tot, 
           colour =method)) +
  geom_line()+
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10() + 
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  #geom_errorbar(aes(ymin = distance_tot_q25, 
  #                  ymax = distance_tot_q75), 
  #              width = 0.1, alpha=0.7,
  #              linewidth=1.1, position = position_dodge(0.05))+
  facet_grid(theta_strength~ n, scales = "free") +
  xlab("p (Dimension)") + 
  ylab(expression("Subspace Distance")) +
  labs(colour="Method") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


labels_n = c()

unique(summ$n)
unique(summ$p1)
unique(summ$r_pca)
unique(summ$r)
ggplot(summ %>% filter(theta_strength == "medium",
                       r==2,
                       n==500,
                       overlapping_amount == 0,
                       method %in% legend_order)) +
  geom_line(aes(x=prop_missing, 
                y = distance_tot, 
                colour =method), linewidth=1.2)+
  scale_color_manual(values =my_colors, 
                     breaks =legend_order) +
  scale_y_log10() + 
  facet_grid(p1~r_pca, scales = "free")


ggplot(summ %>% filter(n == 200, 
                       r_pca == 0, r==2,
                       overlapping_amount == 0,
                       lambda<0.5)) +
  geom_line(aes(x=prop_missing, 
                y = sinTheta_tot, 
                colour =method))+
  facet_wrap(theta_strength~ p1, scales = "free")

ggplot(summ) +
  geom_line(aes(x=prop_missing, y = 1/p1 * prediction_U, colour =method))+
  facet_grid(r ~ p1)

ggplot(summ) +
  geom_line(aes(x=prop_missing, y = prediction_U, colour =method))+
  facet_grid(p1~ r)

summ = result %>% group_by(comp, method, prop, p, noise) %>% 
  summarize(mses_sd = sd(mses), mses = mean(mses),
            cors_sd = sd(cors), cors = mean(cors),
            Uangs_sd = sd(Uangs), Uangs = mean(Uangs), 
            Vangs_sd = sd(Vangs), Vangs = mean(Vangs),  
            XUangs_sd = sd(XUangs), XUangs = mean(XUangs),  
            YVangs_sd = sd(YVangs), YVangs = mean(YVangs)) %>% 
  ungroup() %>%
  mutate(component = as.factor(comp)) %>%
  dplyr::select(-comp) %>%
  filter(component %in% c(1,2,3), prop <=0.3)



#plot all metrics
summ_long = cbind(summ %>% dplyr::select(component, method, prop, p, noise, cors, mses, Uangs, Vangs) %>% 
                    rename(`canonical correlation` = cors, `MSE b/w variates` = mses, `angle b/w U spaces` = Uangs, `angle b/w V spaces` = Vangs) %>% 
                    pivot_longer(6:9) %>%
                    rename(metric = name),
                  summ %>% dplyr::select(cors_sd, mses_sd, Uangs_sd, Vangs_sd) %>% 
                    pivot_longer(1:4) %>% 
                    dplyr::select(value) %>%
                    rename(value_sd = value)) %>%
  mutate(metric =  factor(metric, levels = c("canonical correlation", "MSE b/w variates", "angle b/w U spaces", "angle b/w V spaces")))

summ_long %>% filter(p == 10, noise == 0.1) %>%
  ggplot(aes(prop, value, color = method))+
  geom_point()+
  geom_ribbon(aes(ymin = value-value_sd, ymax = value+value_sd, fill = method), alpha = 0.2, color = NA)+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("")+
  facet_grid(metric~component, labeller = labeller(component = label_both), scale = "free")
ggsave("Plots/simulation-RRR-metrics.pdf", width = 6, height = 6)


# #plot mse
# summ %>% filter(p == 10, noise == 0.1) %>%
#   ggplot(aes(prop, mses, color = method))+
#   geom_point()+
#   geom_ribbon(aes(ymin = mses-mses_sd, ymax = mses+mses_sd, fill = method), alpha = 0.2, color = NA)+
#   geom_line()+
#   xlab("mean squared error (MSE)")+
#   ylab("correlation")+
#   facet_wrap(~component, labeller = labeller(component = label_both))
# ggsave("Results/simulation-RRR-mse.pdf", width = 6, height = 2.5)
# 
# 
# #plot correlation
# summ %>% filter(p == 10, noise == 1) %>%
#   ggplot(aes(prop, cors, color = method))+
#   geom_point()+
#   geom_ribbon(aes(ymin = cors-cors_sd, ymax = cors+cors_sd, fill = method), alpha = 0.2, color = NA)+
#   geom_line()+
#   xlab("proportion of missing values")+
#   ylab("correlation")+
#   facet_wrap(~component, labeller = labeller(component = label_both))
# ggsave("Results/simulation-RRR-cor.pdf", width = 6, height = 2.5)
# 
# #plot V and U angles
# summ %>% filter(component %in% c(1,2,3), p == 10, noise == 0.1) %>%
#   ggplot(aes(prop, Uangs, color = method))+
#   geom_point()+
#   geom_ribbon(aes(ymin = Uangs-Uangs_sd, ymax = Uangs+Uangs_sd, fill = method), alpha = 0.2, color = NA)+
#   geom_line()+
#   xlab("proportion of missing values")+
#   ylab("angle between U coefficients")+
#   facet_wrap(~component, labeller = labeller(component = label_both))
# ggsave("Results/simulation-RRR-uangle.pdf", width = 6, height = 2.5)
# 
# summ %>% filter(component %in% c(1,2,3), p == 10, noise == 0.1) %>%
#   ggplot(aes(prop, Vangs, color = method))+
#   geom_point()+
#   geom_ribbon(aes(ymin = Vangs-Vangs_sd, ymax = Vangs+Vangs_sd, fill = method), alpha = 0.2, color = NA)+
#   geom_line()+
#   xlab("proportion of missing values")+
#   ylab("angle between V coefficients")+
#   facet_wrap(~component, labeller = labeller(component = label_both))
# ggsave("Results/simulation-RRR-vangle.pdf", width = 6, height = 2.5)


#plot mse full comparison
summ %>% filter(noise != 10) %>%
  ggplot(aes(prop, mses, color = component, linetype = method))+
  geom_point()+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("mean squared error (MSE)")+
  facet_grid(noise~p, labeller = labeller(p = label_both, noise = label_both), scale = "free")
ggsave("Plots/simulation-RRR-mse-full.pdf",  width = 7, height = 7)

summ %>%filter(noise != 10) %>%
  ggplot(aes(prop, cors, color = component, linetype = method))+
  geom_point()+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("correlation")+
  facet_grid(noise~p, labeller = labeller(p = label_both, noise = label_both), scale = "free")
ggsave("Plots/simulation-RRR-cor-full.pdf", width = 7, height = 7)

summ %>% filter(noise != 10) %>%
  ggplot(aes(prop, Uangs, color = component, linetype = method))+
  geom_point()+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("angle between U coefficients")+
  facet_grid(noise~p, labeller = labeller(p = label_both, noise = label_both), scale = "free")
ggsave("Plots/simulation-RRR-uangle-full.pdf", width = 7, height = 7)

summ %>% filter(noise != 10) %>%
  ggplot(aes(prop, Vangs, color = component, linetype = method))+
  geom_point()+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("angle between V coefficients")+
  facet_grid(noise~p, labeller = labeller(p = label_both, noise = label_both), scale = "free")
ggsave("Plots/simulation-RRR-vangle-full.pdf", width = 7, height = 7)


