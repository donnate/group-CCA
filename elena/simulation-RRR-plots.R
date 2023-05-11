
result = read.csv("Results/simulation-RRR.csv", header = TRUE)

summ = result %>% group_by(comp, method, prop, p, noise) %>% 
  summarize(mses_sd = sd(mses), mses = mean(mses),
            cors_sd = sd(cors), cors = mean(cors),
            Uangs_sd = sd(Uangs), Uangs = mean(Uangs), 
            Vangs_sd = sd(Vangs), Vangs = mean(Vangs),  
            XUangs_sd = sd(XUangs), XUangs = mean(XUangs),  
            YVangs_sd = sd(YVangs), YVangs = mean(YVangs)) %>% 
  ungroup() %>%
  mutate(component = as.factor(comp)) %>%
  filter(component %in% c(1,2,3), prop <=0.3)


#plot mse
summ %>% filter(p == 10, noise == 1) %>%
  ggplot(aes(prop, mses, color = method))+
  geom_point()+
  geom_ribbon(aes(ymin = mses-mses_sd, ymax = mses+mses_sd, fill = method), alpha = 0.2, color = NA)+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("mean squared error (MSE)")+
  facet_wrap(~component, labeller = labeller(component = label_both))
ggsave("Results/simulation-RRR-mse.pdf", width = 6, height = 3)

#plot correlation
summ %>% filter(p == 10, noise == 1) %>%
  ggplot(aes(prop, cors, color = method))+
  geom_point()+
  geom_ribbon(aes(ymin = cors-cors_sd, ymax = cors+cors_sd, fill = method), alpha = 0.2, color = NA)+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("correlation")+
  facet_wrap(~component, labeller = labeller(component = label_both))
ggsave("Results/simulation-RRR-mse.pdf", width = 6, height = 3)

#plot V and U angles
summ %>% filter(component %in% c(1,2,3), p == 10, noise == 1, prop <=0.3) %>%
  ggplot(aes(prop, Uangs, color = method))+
  geom_point()+
  geom_ribbon(aes(ymin = Uangs-Uangs_sd, ymax = Uangs+Uangs_sd, fill = method), alpha = 0.2, color = NA)+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("angle between U coefficients")+
  facet_wrap(~component, labeller = labeller(component = label_both))
ggsave("Results/simulation-RRR-uangle.pdf", width = 6, height = 3)

summ %>% filter(component %in% c(1,2,3), p == 10, noise == 1, prop <=0.3) %>%
  ggplot(aes(prop, Vangs, color = method))+
  geom_point()+
  geom_ribbon(aes(ymin = Vangs-Vangs_sd, ymax = Vangs+Vangs_sd, fill = method), alpha = 0.2, color = NA)+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("angle between V coefficients")+
  facet_wrap(~component, labeller = labeller(component = label_both))
ggsave("Results/simulation-RRR-vangle.pdf", width = 6, height = 3)


#plot mse full comparison
summ %>% filter(p != 30, noise %in% c(0.1, 1, 2)) %>%
  ggplot(aes(prop, mses, color = component, linetype = method))+
  geom_point()+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("mean squared error (MSE)")+
  facet_grid(noise~p, labeller = labeller(p = label_both, noise = label_both), scale = "free")
ggsave("Results/simulation-RRR-mse-full.pdf", width = 8, height = 6)

summ %>% filter(p != 30, noise %in% c(0.1, 1, 2)) %>%
  ggplot(aes(prop, cors, color = component, linetype = method))+
  geom_point()+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("correlation")+
  facet_grid(noise~p, labeller = labeller(p = label_both, noise = label_both), scale = "free")
ggsave("Results/simulation-RRR-cor-full.pdf", width = 8, height = 6)

summ %>% filter(p != 30, noise %in% c(0.1, 1, 2)) %>%
  ggplot(aes(prop, Uangs, color = component, linetype = method))+
  geom_point()+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("angle between U coefficients")+
  facet_grid(noise~p, labeller = labeller(p = label_both), scale = "free")
ggsave("Results/simulation-RRR-uangle-full.pdf", width = 8, height = 6)

summ %>% filter(p != 30, noise %in% c(0.1, 1, 2)) %>%
  ggplot(aes(prop, Vangs, color = component, linetype = method))+
  geom_point()+
  geom_line()+
  xlab("proportion of missing values")+
  ylab("angle between V coefficients")+
  facet_grid(noise~p, labeller = labeller(p = label_both), scale = "free")
ggsave("Results/simulation-RRR-vangle-full.pdf", width = 8, height = 6)


