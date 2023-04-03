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

# Calculate total carbon and nitrogen mineralization across the webs ----

# Create a list to store all the results:
results = vector("list", length = 16)

reps_per_web = 10

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

pullconsump <- function(x) c(Cmin = sum(x$Cmin), Nmin = sum(x$Nmin))

pullwhomins <- function(x){
  t1 = whomineralizes(x$usin)
  colSums(t1[t1$ID %in% c("Predaceousmites"),c(2,4)])
  }

comtrosp(puGA[[1]]$usin, selected = c("Bacteriophagousnematodes","Fungivorousnematodes","Omnivorousnematodes", "Phytophagousnematodes"))

puGA = cbind(do.call("rbind",lapply(puGA, pullconsump)),
             do.call("rbind",lapply(puGA, pullwhomins)))

results[[1]] = cbind(as.data.frame(puGA), Web = "GA")


biomassescur = biomasses %>% filter(name == "A_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGA,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = cbind(do.call("rbind",lapply(puGA, pullconsump)),
             do.call("rbind",lapply(puGA, pullwhomins)))
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


puGA = cbind(do.call("rbind",lapply(puGA, pullconsump)),
             do.call("rbind",lapply(puGA, pullwhomins)))
results[[3]] = cbind(as.data.frame(puGA), Web = "GB")


biomassescur = biomasses %>% filter(name == "B_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGB,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = cbind(do.call("rbind",lapply(puGA, pullconsump)),
             do.call("rbind",lapply(puGA, pullwhomins)))
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


puGA = cbind(do.call("rbind",lapply(puGA, pullconsump)),
             do.call("rbind",lapply(puGA, pullwhomins)))
results[[5]] = cbind(as.data.frame(puGA), Web = "GC")


biomassescur = biomasses %>% filter(name == "C_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGC,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = cbind(do.call("rbind",lapply(puGA, pullconsump)),
             do.call("rbind",lapply(puGA, pullwhomins)))
results[[6]] = cbind(as.data.frame(puGA), Web = "UGC")

# Holtkamp ------


# Load in errors, convert SD to VAR
errmes <- as.matrix(read.csv("Data/Holtkamp/biomass_sd.csv", header = T, row.names = 1))^2


# Simulate the direct and indirect effects:
errmes2 <- as.matrix(errmes[,"Young"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Young,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

puGA = cbind(do.call("rbind",lapply(puGA, pullconsump)),
             do.call("rbind",lapply(puGA, pullwhomins)))
results[[7]] = cbind(as.data.frame(puGA), Web = "Young")

errmes2 <- as.matrix(errmes[,"Mid"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Mid,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

puGA = cbind(do.call("rbind",lapply(puGA, pullconsump)),
             do.call("rbind",lapply(puGA, pullwhomins)))
results[[8]] = cbind(as.data.frame(puGA), Web = "Mid")

errmes2 <- as.matrix(errmes[-3,"Old"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Old,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

puGA = cbind(do.call("rbind",lapply(puGA, pullconsump)),
             do.call("rbind",lapply(puGA, pullwhomins)))
results[[9]] = cbind(as.data.frame(puGA), Web = "Old")

errmes2 <- as.matrix(errmes[,"Health"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Heathland,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

puGA = cbind(do.call("rbind",lapply(puGA, pullconsump)),
             do.call("rbind",lapply(puGA, pullwhomins)))
results[[10]] = cbind(as.data.frame(puGA), Web = "Heathland")

# Summarize

results2 = do.call("rbind", results)

results2 = results2 %>% 
  mutate(Manuscript = ifelse(Web %in% c("Young", "Mid", "Old", "Heathland"), "Holtkamp et al.", "Andres et al.")) %>%
  mutate(Web = factor(Web, levels = c("Young", "Mid", "Old", "Heathland", "GA", "UGA", "GB", "UGB", "GC", "UGC")))

cowplot::plot_grid(
  results2 %>%
    tibble() %>%
    ggplot(aes(x = Web, y = Cmin, color = Manuscript)) + geom_boxplot() + theme_classic(),
  
  results2 %>%
    tibble() %>%
    ggplot(aes(x = Web, y = Cmin, color = Manuscript)) + geom_boxplot() + theme_classic(),
  
  results2 %>%
    tibble() %>%
    ggplot(aes(x = Web, y = Cmin, color = Manuscript)) + geom_boxplot() + theme_classic(),
  
  results2 %>%
    tibble() %>%
    ggplot(aes(x = Web, y = Cmin, color = Manuscript)) + geom_boxplot() + theme_classic()
)

png("Plots/demonstration_parameteruncertainty.png", width = 8, height = 5, units = "in", res = 600)
results2 %>%
  tibble() %>%
  pivot_longer(-Web) %>%
  separate(Web, into = c("Grazed2", "Site2"), sep = -1) %>%
  left_join(
    tibble(Grazed2 = c("G", "UG"),
           Treatment = c("Grazed", "Ungrazed")) 
  ) %>%
  left_join(
    tibble(Site2 = c("A", "B", "C", "S"),
           Site = factor(c("A", "B", "C", "All"), levels = c("A", "B", "C", "All")))
  ) %>%
  mutate(Site = paste0(Site, "-", Grazed2)) %>%
  ggplot(aes(x = Site, y = value, color = Treatment)) + geom_boxplot() + facet_wrap(.~name, scales = "free") + theme_classic() + ylab(parse(text = "Nutrient~mineralization~(kg[Nutrient]~ha^-1~yr^-1)"))
dev.off()