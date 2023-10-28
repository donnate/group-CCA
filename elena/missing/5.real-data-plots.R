library(tidyr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(scales)

#nutrimouse
result = read.csv("Fits/real-data-nutrimouse-ranks-RRR.csv", header = TRUE)
seeds = max(result$seed)
result %>% filter(comp == 1, lambda2 == 0.001, prop == 0.3) %>%
  ggplot(aes(log(lambda1,10), cors, color = method))+
  geom_point() # check variation in scores - large for nutrimouse

result %>% group_by(comp, method, lambda1, lambda2, prop) %>% filter(comp <= 5) %>%
  summarize(cors = mean(cors)) %>%
  ggplot(aes(log(lambda1,10), cors, color = factor(lambda2), linetype = method))+
  facet_grid(prop~comp, labeller = labeller(prop = label_both, component = label_both))+
  geom_point()+
  geom_line()

result %>% group_by(comp, method, lambda1, lambda2, prop) %>% filter(comp <= 3) %>%
  summarize(cors = mean(cors)) %>%
  group_by(comp, method, prop) %>%
  summarize(cors = max(cors)) %>%
  ggplot(aes(prop, cors, color = method))+
  facet_grid(~comp, labeller = labeller(prop = label_both, component = label_both))+
  geom_point()+
  geom_line()

summ = result %>% group_by(comp, method, lambda1, lambda2, prop) %>% 
  summarize(cors_sd = sd(cors)/sqrt(seeds), cors = mean(cors),
            mses_sd = sd(mses)/sqrt(seeds), mses = mean(mses), 
            angles_sd = sd(angles)/sqrt(seeds), angles = mean(angles)) %>% 
  ungroup() %>%
  mutate(component = comp) %>%
  dplyr::select(-comp)

#metabrics
result = read.csv("Fits/real-data-metabrics-RRR-cv-randomized.csv", header = TRUE)
result %>% filter(comp == 1, prop == 0.25) %>%
  ggplot(aes(log(lambda1,10), cors, color = method))+
  geom_point()
nfold = max(result$fold)

resultsumm = result %>% filter(comp <= 3) %>% 
  group_by(method, lambda1, comp, prop) %>% 
  summarise(cors_sd = sd(cors)/sqrt(nfold), cors = mean(cors), 
            angles_sd = sd(angles)/sqrt(nfold), angles = mean(angles)) 

resultsumm %>%
  ggplot(aes(log(lambda1,10), cors, color = method, fill = method))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin = cors-cors_sd, ymax = cors+cors_sd), color = NA, alpha = 0.2)+
  facet_grid(prop~comp, scales = "free")
resultsumm %>%
  ggplot(aes(log(lambda1,10), angles, color = method, fill = method))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin = angles-angles_sd, ymax = angles+angles_sd), color = NA, alpha = 0.2)+
  facet_wrap(~comp, scales = "free")


