# Test multivariate distribution:

library(soilfoodwebs)
library(tidyverse)

set.seed(2323)
Prey1 = runif(10, min = 0.5, max = 1.5)
Prey2 = runif(10, min = 0.5, max = 1.5)
Predator = Prey1*0.2 + Prey2*0.01 + runif(10, min = -0.01, max = 0.01)

biomaesses = data.frame(Prey1 = Prey1, Prey2 = Prey2, Predator = Predator)

cov(biomaesses)

plot(Predator~Prey1)

# Create a data frame of feeding relationships:
feedinglist = data.frame(
Predator = c("Pred", "Pred"),
Prey = c("Prey1", "Prey2"),
Preference = c(1,1.2))

# Create a data frame of properties for each species:
properties = data.frame(ID = c("Pred", "Prey1", "Prey2"), # Name
d = c(1,3,0.5), # Death rate
a = c(0.61,0.65,0.45), # Assimilation efficiency
p = c(0.5,0.4,0.3), # Production efficiency for carbon
B = c(mean(Predator),mean(Prey1),mean(Prey2)), # Biomass
CN = c(4.5,4.8,5), # Carbon to nitrogen ratio
DetritusRecycling = c(0,0,0), # proportion of detritus recycling
isDetritus = c(0,0,0), # Boolean: Is this pool detritus?
isPlant = c(0,0,0), # Boolean: Is this pool a plant?
canIMM = c(0,0,0)) # Boolean: Can the pool immobilize inorganic nitrogen?

# Build the food web:
comm = build_foodweb(feedinglist, properties)

comana(comm)

output = exp(MASS::mvrnorm(n = 1000, mu = c(mean(log(Predator)),mean(log(Prey1)),mean(log(Prey2))), Sigma = cov(log(biomaesses))))

results = vector('list', 1000)

for(i in 1:dim(output)[1]){
  comm2 = comm
  comm2$prop$B = output[i,]
  results[[i]] = comana(comm2)$consumption
}

output2 = cbind(
  cbind(exp(rnorm(1000, mean = mean(log(Predator)), sd = sd(log(Predator)))),
        exp(rnorm(1000, mean = mean(log(Prey1)), sd = sd(log(Prey1))))),
  exp(rnorm(1000, mean = mean(log(Prey2)), sd = sd(log(Prey2))))
)

results2 = vector('list', 1000)

for(i in 1:dim(output2)[1]){
  comm2 = comm
  comm2$prop$B = output2[i,]
  results2[[i]] = comana(comm2)$consumption
}

data.frame(
  rbind(
    cbind(do.call('rbind', results), Type = 1),
    cbind(do.call('rbind', results2), Type = 2)
  )
) %>%
  tibble() %>%
  pivot_longer(!Type) %>%
  ggplot(aes(x = value, color = paste(Type))) + geom_density() + facet_wrap(.~name, scales = "free")

data.frame(
  rbind(
    cbind(do.call('rbind', results), Type = 1),
    cbind(do.call('rbind', results2), Type = 2)
  )
) %>%
  tibble() %>%
  pivot_longer(!Type) %>%
  group_by(Type, name) %>%
  summarize(mu = mean(value), std = sd(value))
