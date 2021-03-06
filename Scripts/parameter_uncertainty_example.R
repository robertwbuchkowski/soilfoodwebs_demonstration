# Script for the paper presenting soilfoodwebs
# Author: Robert W. Buchkowski
# Date created: Sept. 29/ 2021

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

reps_per_web = 1000

setseed(102132)

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

pullconsump <- function(x) x$consumption

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

results2 = do.call("rbind", results)

png("Plots/demonstration_parameteruncertainty.png", width = 8, height = 5, units = "in", res = 600)
results2 %>%
  tibble() %>%
  pivot_longer(-Web) %>%
  filter(grepl("nematode", name)) %>%
  separate(Web, into = c("Grazed", "Site"), sep = -1) %>%
  separate(name, into = c("A", "B"), sep = -9) %>%
  mutate(name = paste(A, B)) %>%
  ggplot(aes(x = Site, y = value, color = Grazed)) + stat_summary(fun.data = "mean_cl_boot") + facet_wrap(.~name, scales = "free") + theme_classic() + ylab(parse(text = "Carbon~consumption~(kg[C]~ha^-1~yr^-1)"))
dev.off()
