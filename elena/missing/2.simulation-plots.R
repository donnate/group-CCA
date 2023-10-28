library(tidyverse)
library(ggplot2)
#setwd("~/Documents/group-CCA/elena/")
result = read.csv("missing CCA code/simulation-RRR-results-sparse.csv", header = TRUE)

result["lambda"] = sapply(result$method, function(x){ifelse(is.null(strfind(x, "Alt-"))==FALSE,
                                                           as.numeric(strsplit(x, "Alt-")[[1]][2]),
                                                           NA)})
result["method_type"] = sapply(result$method, function(x){ifelse(is.null(strfind(x, "Alt-"))==FALSE,
                                                            "Alt",
                                                            x)})
summ = result %>% group_by(n, p1, p2, r, r_pca, 
                           overlapping_amount, noise, method_type, lambda,
                           prop_missing) %>% 
  summarize_all(median) %>% 
  ungroup() %>%
  arrange(n, p1, p2, r, noise, prop_missing, distance_U)


ggplot(result %>% filter(n == 100, method_type=="Alt"),
       aes(x=lambda, y =  1/p1 * distance_tot, 
           colour =method)) +
  geom_point(alpha=0.1)+
  geom_smooth()+
  scale_x_log10()
  facet_wrap(r_pca ~ p1, scales = "free")


ggplot(summ %>% filter(n == 200)) +
  geom_line(aes(x=prop_missing, y = sinTheta_tot, colour =method))+
  facet_wrap(r_pca + r ~ p1, scales = "free")

ggplot(summ) +
  geom_line(aes(x=prop_missing, y = 1/p1 * prediction_U, colour =method))+
  facet_grid(r_pca + r ~ p1)

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


