library(ggplot2)
library(dplyr)
library(fields)
library(CCA)
library(tidyr)
library(zoo)
library(pracma)


generate =  function(n, p, q, s, prop){
  X = matrix(rnorm(n*p), n, p)
  R = matrix(rnorm(p*q), p, q)
  E = matrix(rnorm(n*q, 0, s), n, q)
  Y = X %*% R + E
  X = scale(X, center = TRUE, scale = FALSE)
  Y = scale(Y, center = TRUE, scale = FALSE)
  maskX = sample(1:(n*p), n*p*prop)
  Xna = X
  Xna[maskX] = NA
  maskY = sample(1:(n*q), n*q*prop)
  Yna = Y
  Yna[maskY] = NA
  return(list(X = X, Y = Y, Xna = Xna, Yna = Yna))
}

# #Sanity checks
# #generate data
# gen = generate(n = 100, p = 10, q = 5, s = 1, prop = 0.1)
#
# #For non-missing values: CCA and RRR are equivalent
# X = gen$X
# Y = gen$Y
# rrr = CCA_RRR(X, Y, ncol(Y))
# rrr$cors
# cca = cc(X, Y)
# cca$cor
# U0 = cca$xcoef
# V0 = cca$ycoef
# evaluate(X, Y, rrr$U, rrr$V, U0, V0)
# 
# #For missing values
# Xna = gen$Xna
# Yna = gen$Yna
# rrr = CCA_RRR(Xna, Yna)
# U = rrr$U
# V = rrr$V
# rrr$cors
# evaluate(X, Y, rrr$U, rrr$V, U0, V0)


#Simulation for missing values in both X and Y

n = 100
seeds = 1:100
ps = c(5, 10, 20, 30, 50, 70)
ss = c(0.01, 0.1, 0.5, 1, 2)
props = seq(0, 0.7, 0.05)

result = c()

for(seed in 1:100){
for(s in ss){
for(p in ps){
  q = p
  for(prop in props){
      set.seed(seed)
      gen = generate(n, p, q, s, prop)
      X = gen$X
      Y = gen$Y
      Xna = gen$Xna
      Yna = gen$Yna
      Ximp = data.frame(Xna) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE)))
      Yimp = data.frame(Yna) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE)))
      
      cca0 = cc(X, Y)
      rrr = CCA_RRR(Xna, Yna, k = 5, eps = 1e-10)
      simple = CCA_simple(Xna, Yna, k = 5, eps = 1e-10)
      result = rbind(result, data.frame(evaluate(X, Y, rrr$U, rrr$V, cca0$xcoef, cca0$ycoef), 
                              p,  noise = s, method = "RRR",  prop, seed))
      result = rbind(result, data.frame(evaluate(X, Y, simple$U,simple$V, cca0$xcoef, cca0$ycoef), 
                                        p,  noise = s, method = "simple",  prop, seed))
      ### Imputed
      cca = cc(Ximp, Yimp)
      result = rbind(result, data.frame(evaluate(X, Y, cca$xcoef, cca$ycoef, cca0$xcoef, cca0$ycoef),
                              p, noise = s, method = "CCA", prop, seed))
    }
  }
}
}
write.csv(result, "Results/simulation-RRR.csv", row.names = F)

test = result %>% 
  group_by(noise, method, prop ) %>%
  summarise_all(mean)


