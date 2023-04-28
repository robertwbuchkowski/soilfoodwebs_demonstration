# Script for the paper presenting soilfoodwebs
# Author: Robert W. Buchkowski
# Date created: April, 28/ 2021

# If need you to install new version from github:
if(F){
  devtools::install_github("robertwbuchkowski/soilfoodwebs", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual"))
}

# Load in the libraries:
if (!require("pacman")) install.packages("pacman")
p_load(soilfoodwebs,tidyverse)

# Calculate the contributions to mineralization of trophic groups across the food webs----

# Create a list to store all the results:
results = vector("list", length = 16)

reps_per_web = 100

set.seed(102132)

# ... Andres et al. 2016----

# Because only biomass changes between the communities, we can use this baseline community and cycle through the mean biomass random draws.

# Load in the file of biomasses and standard errors:
biomasses = read_csv("Data/Andres/biomass_sd.csv") %>%
  # Convert to gamma distribution:
  mutate(n = ifelse(grepl("SOM", ID),3,10)) %>%
  mutate(Bvar = n*SE^2) %>%
  filter(ID != "Totalinvertebrate biomass")

biomasses = biomasses %>%
  filter(ID == "ActiveSOM") %>%
  mutate(ID = "Roots") %>%
  mutate(Bvar = (3000*0.2)^2) %>%
  bind_rows(
    biomasses
  ) %>%
  select(ID, name, Bvar) %>% rename(B = Bvar)

# SITE A
biomassescur = biomasses %>% filter(name == "A_G") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$GA,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

pullconsump <- function(x){c(
  Cmin = sum(x$Cmin),
  Cflux = sum(x$fmat),
  BEC = x$consumption["Bacteria"],
  FEC = x$consumption["Fungi"]
  
)
  
}

puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[1]] = cbind(as.data.frame(puGA), Web = "GA")


biomassescur = biomasses %>% filter(name == "A_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGA,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[2]] = cbind(as.data.frame(puGA), Web = "UGA")

# SITE B
biomassescur = biomasses %>% filter(name == "B_G") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$GB,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[3]] = cbind(as.data.frame(puGA), Web = "GB")


biomassescur = biomasses %>% filter(name == "B_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGB,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[4]] = cbind(as.data.frame(puGA), Web = "UGB")

# SITE C
biomassescur = biomasses %>% filter(name == "C_G") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$GC,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[5]] = cbind(as.data.frame(puGA), Web = "GC")


biomassescur = biomasses %>% filter(name == "C_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGC,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[6]] = cbind(as.data.frame(puGA), Web = "UGC")

# META-SITE
tempsite = Andres2016$GA
tempsite$prop$B = apply(
  matrix(c(Andres2016$GA$prop$B,Andres2016$GB$prop$B,Andres2016$GC$prop$B), nrow = 21, ncol = 3),
  1,
  mean
)

temperror = matrix(data = apply(
  matrix(c(Andres2016$GA$prop$B,Andres2016$GB$prop$B,Andres2016$GC$prop$B), nrow = 21, ncol = 3),
  1,var), nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(Andres2016$GA$prop$ID, "B"))
temperror["Roots","B"] = (3000*0.2)^2

puGA = parameter_uncertainty(usin = tempsite,
                             errormeasure = temperror,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[7]] = cbind(as.data.frame(puGA), Web = "GS")

tempsite = Andres2016$UGA
tempsite$prop$B = apply(
  matrix(c(Andres2016$UGA$prop$B,Andres2016$UGB$prop$B,Andres2016$UGC$prop$B), nrow = 21, ncol = 3),
  1,
  mean
)

temperror = matrix(data = apply(
  matrix(c(Andres2016$UGA$prop$B,Andres2016$UGB$prop$B,Andres2016$UGC$prop$B), nrow = 21, ncol = 3),
  1,var), nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(Andres2016$GA$prop$ID, "B"))
temperror["Roots","B"] = (3000*0.2)^2

puGA = parameter_uncertainty(usin = tempsite,
                             errormeasure = temperror,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[8]] = cbind(as.data.frame(puGA), Web = "UGS")

rm(tempsite, temperror)

results2 = do.call("rbind", results)

png("Plots/demonstration_parameteruncertainty.png", width = 8, height = 5, units = "in", res = 600)
results2 %>%
  tibble() %>%
  pivot_longer(-Web) %>%
  separate(Web, into = c("Grazed2", "Site2"), sep = -1) %>%
  left_join(
    tibble(Grazed2 = c("UG", "G"),
           Treatment = c("Grazed", "Ungrazed")) 
  ) %>%
  left_join(
    tibble(Site2 = c("A", "B", "C", "S"),
           Site = factor(c("A", "B", "C", "All"), levels = c("A", "B", "C", "All")))
  ) %>%
  ggplot(aes(x = Site, y = value, color = Treatment)) + geom_boxplot() + facet_wrap(.~name, scales = "free") + theme_classic() + ylab(parse(text = "Carbon~consumption~(kg[C]~ha^-1~yr^-1)"))
dev.off()
