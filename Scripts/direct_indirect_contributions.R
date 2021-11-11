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
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)


puGA = do.call("rbind",puGA)
results[[1]] = cbind(puGA, Web = "GA")


biomassescur = biomasses %>% filter(name == "A_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGA,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)


puGA = do.call("rbind",puGA)
results[[2]] = cbind(puGA, Web = "UGA")

# SITE B
biomassescur = biomasses %>% filter(name == "B_G") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$GB,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)


puGA = do.call("rbind",puGA)
results[[3]] = cbind(puGA, Web = "GB")


biomassescur = biomasses %>% filter(name == "B_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGB,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)


puGA = do.call("rbind",puGA)
results[[4]] = cbind(puGA, Web = "UGB")

# SITE C
biomassescur = biomasses %>% filter(name == "C_G") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$GC,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)


puGA = do.call("rbind",puGA)
results[[5]] = cbind(puGA, Web = "GC")


biomassescur = biomasses %>% filter(name == "C_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGC,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)


puGA = do.call("rbind",puGA)
results[[6]] = cbind(puGA, Web = "UGC")

# ... Koltz et al. web ----
errmes <- as.matrix(read.csv("Data/Koltz/biomass_sd.csv", header = T, row.names = 1))^2

# Identify and remove zero biomass nodes
puGA = parameter_uncertainty(usin = checkcomm(Koltz2018),
                             errormeasure = errmes,
                             errortype = "Variance",
                             fcntorun = "whomineralizes",
                             returnprops = F, 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[7]] = cbind(puGA, Web = "Koltz")

# ... Hunt1987 web ----

puGA = parameter_uncertainty(usin = Hunt1987,
                             errormeasure = 0.2,
                             errortype = "CV",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[8]] = cbind(puGA, Web = "CPER")


# ... Holtkamp et al. 2011 webs ----

# Load in errors, convert SD to VAR
errmes <- as.matrix(read.csv("Data/Holtkamp/biomass_sd.csv", header = T, row.names = 1))^2

errmes2 <- as.matrix(errmes[,"Young"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Young,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[9]] = cbind(puGA, Web = "Young")

errmes2 <- as.matrix(errmes[,"Mid"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Mid,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[10]] = cbind(puGA, Web = "Mid")

errmes2 <- as.matrix(errmes[-3,"Old"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Old,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[11]] = cbind(puGA, Web = "Old")

errmes2 <- as.matrix(errmes[,"Health"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Heathland,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[12]] = cbind(puGA, Web = "Heathland")


# ... deRuiter et al. 1994 webs ----

puGA = parameter_uncertainty(usin = deRuiter1994$CON,
                             errormeasure = 0.2,
                             errortype = "CV",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[13]] = cbind(puGA, Web = "CON")

puGA = parameter_uncertainty(usin = deRuiter1994$INT,
                             errormeasure = 0.2,
                             errortype = "CV",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[14]] = cbind(puGA, Web = "INT")

puGA = parameter_uncertainty(usin = deRuiter1994$CON10,
                             errormeasure = 0.2,
                             errortype = "CV",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[15]] = cbind(puGA, Web = "CON10")

puGA = parameter_uncertainty(usin = deRuiter1994$INT10,
                             errormeasure = 0.2,
                             errortype = "CV",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[16]] = cbind(puGA, Web = "INT10")

# Analyze the results of the full study ----

results2 <- do.call("rbind", results)

# Organize trophic groups:

rbind(
  cbind(TL = TLcheddar(Andres2016$GA$imat), Web = "Andres2016"),
  cbind(TL = TLcheddar(Hunt1987$imat), Web = "CPER"),
  cbind(TL = TLcheddar(Koltz2018$imat), Web = "Koltz2018"),
  cbind(TL = TLcheddar(deRuiter1994$CON$imat), Web = "deRuiter1994"),
  cbind(TL = TLcheddar(Holtkamp2011$Young$imat), Web = "Holtkamp2011")
)

results3 = results2 %>% tibble() %>%
  left_join(
    tibble(Web = unique(results2$Web),
           Web2 = c(rep("Andres2016", 6), "Koltz2018", "CPER", rep("Holtkamp2011",4), rep("deRuiter1994", 4))),
    by = "Web"
  ) %>%
  left_join(
    read.csv("Data/matching_trophicspecies.csv"),
    by = c("Web2", "ID")
  ) %>%
  pivot_longer(DirectC:IndirectN) %>%
  filter(TL > 1.0)

# Check which types are only positive
results3 %>%
  mutate(value = ifelse(abs(value) < 1e-8, 0, value)) %>%
  group_by(Group2, name) %>%
  summarize(a = min(value), b = max(value)) %>%
  filter(a*b >=0) %>% View()

png("Plots/demonstration_effects.png", width = 8, height = 5, units = "in", res = 600)
results3 %>%
  group_by(Group2, name) %>%
  summarize(N = n()) %>%
  filter(N > 1000) %>%
  filter(Group2 != 'Herbivores') %>%
  select(-N) %>%
  left_join(
    results3, by = c("Group2", "name")
  ) %>%
  left_join(
    results3 %>%
      mutate(value = ifelse(abs(value) < 1e-8, 0, value)) %>%
      group_by(Group2, name) %>%
      summarize(a = min(value), b = max(value)) %>%
      filter(a*b >=0) %>%
      select(Group2, name) %>%
      mutate(Positive2 = "Yes"), by = c("Group2", "name")
  ) %>%
  mutate(Positive = ifelse(!is.na(Positive2), "Yes", "No")) %>%
  separate(name, into = c("Effect", "Element"), sep = -1) %>%
  ggplot(aes(x = paste(Effect, Element), y = value*100, color = Positive)) + geom_hline(yintercept = 0, linetype =2) + geom_boxplot(outlier.size = 0.5) + facet_wrap(.~Group2, scales = "free") + ylab("Effect on C or N mineralization  (%)") + xlab("Effect") + scale_color_manual(values = c("blue", "orange"), name = "Strictly Positive or Negative") + theme_classic() + theme(legend.position = "top")
dev.off()
