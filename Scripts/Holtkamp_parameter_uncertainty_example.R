# Script for the paper presenting soilfoodwebs
# Author: Robert W. Buchkowski
# Date created: Aug 8/ 2022

# If need you to install new version from github:
if(F){
  devtools::install_github("robertwbuchkowski/soilfoodwebs", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual"))
}

# Load in the libraries:
if (!require("pacman")) install.packages("pacman")
p_load(soilfoodwebs,tidyverse)

reps_per_web = 500
results = vector("list",4)

# Load in errors, convert SD to VAR
errmes <- as.matrix(read.csv("Data/Holtkamp/biomass_sd.csv", header = T, row.names = 1))^2


# Simulate the direct and indirect effects:
errmes2 <- as.matrix(errmes[,"Young"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Young,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[1]] = cbind(puGA, Web = "Young")

errmes2 <- as.matrix(errmes[,"Mid"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Mid,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[2]] = cbind(puGA, Web = "Mid")

errmes2 <- as.matrix(errmes[-3,"Old"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Old,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[3]] = cbind(puGA, Web = "Old")

errmes2 <- as.matrix(errmes[,"Health"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Heathland,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "whomineralizes", 
                             replicates = reps_per_web)

puGA = do.call("rbind",puGA)
results[[4]] = cbind(puGA, Web = "Heathland")

results = do.call("rbind", results)
# View(results)

png("Plots/holtkamp_effects_parameteruncertainty.png", width = 8, height = 4, units = "in", res = 600)
results %>%
  tibble() %>%
  filter(DirectC != 0) %>%
  filter(grepl("BactNem", ID) | grepl("FungNem", ID)) %>%
  pivot_longer(!ID & !Web) %>%
  filter(grepl("tN", name)) %>%
  separate(name, into = c('name', NA), sep = -1) %>%
  mutate(ID = ifelse(ID == "BactNem", "Bactivorous", "Fungivorous")) %>%
  mutate(Web = factor(Web, levels = c("Young", "Mid", "Old", "Heathland"))) %>%
  filter(name == "Direct") %>%
  ggplot(aes(x = Web, y = value, color = Web)) + geom_boxplot() + facet_wrap(.~ID) + theme_classic() + ylab("Effect magnitude") + ggtitle("Nematode effects on nitrogen mineralization") + xlab("Field Age") + scale_color_discrete(guide = "none")
dev.off()  
  