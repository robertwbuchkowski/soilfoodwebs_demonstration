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

reps_per_web = 1000

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

pullconsump <- function(x){x$consumption}

puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[1]] = cbind(as.data.frame(puGA), Web = "GA", MS = "Andres2016") %>% tibble() %>% pivot_longer(!Web & !MS)


biomassescur = biomasses %>% filter(name == "A_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGA,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[2]] = cbind(as.data.frame(puGA), Web = "UGA", MS = "Andres2016") %>% tibble() %>% pivot_longer(!Web & !MS)

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
results[[3]] = cbind(as.data.frame(puGA), Web = "GB", MS = "Andres2016") %>% tibble() %>% pivot_longer(!Web & !MS)


biomassescur = biomasses %>% filter(name == "B_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGB,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[4]] = cbind(as.data.frame(puGA), Web = "UGB", MS = "Andres2016") %>% tibble() %>% pivot_longer(!Web & !MS)

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
results[[5]] = cbind(as.data.frame(puGA), Web = "GC", MS = "Andres2016") %>% tibble() %>% pivot_longer(!Web & !MS)


biomassescur = biomasses %>% filter(name == "C_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGC,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[6]] = cbind(as.data.frame(puGA), Web = "UGC", MS = "Andres2016") %>% tibble() %>% pivot_longer(!Web & !MS)


# Load in errors, convert SD to VAR
errmes <- as.matrix(read.csv("Data/Holtkamp/biomass_sd.csv", header = T, row.names = 1))^2

errmes2 <- as.matrix(errmes[,"Young"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Young,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[7]] = cbind(as.data.frame(puGA), Web = "Young", MS = "Holtkamp2011") %>% tibble() %>% pivot_longer(!Web & !MS)

errmes2 <- as.matrix(errmes[,"Mid"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Mid,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[8]] = cbind(as.data.frame(puGA), Web = "Mid", MS = "Holtkamp2011") %>% tibble() %>% pivot_longer(!Web & !MS)

errmes <- as.matrix(read.csv("Data/Holtkamp/biomass_sd.csv", header = T, row.names = 1))^2

errmes2 <- as.matrix(errmes[-3,"Old"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Old,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[9]] = cbind(as.data.frame(puGA), Web = "Old", MS = "Holtkamp2011") %>% tibble() %>% pivot_longer(!Web & !MS)

errmes2 <- as.matrix(errmes[,"Health"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Heathland,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

puGA = do.call("rbind",lapply(puGA, pullconsump))
results[[10]] = cbind(as.data.frame(puGA), Web = "Heathland", MS = "Holtkamp2011")  %>% tibble() %>% pivot_longer(!Web & !MS)

results2 = do.call("rbind", results)


results2 = results2 %>%
  left_join(
    read_csv("Data/matching_trophicspecies.csv"), by = c("name" = "ID", "MS" = "Web2"),relationship = "many-to-many"
  ) %>%
  filter(Group2 != "SOM") %>%
  filter(Group2 != "Roots") %>%
  filter(!(Group3 %in% c("Proturans", "Symphyla", "Dilurans", "Enchytraeids", "Ciliates", "CryJuvMite", "PredCol")))

results2 %>%
  filter(Group3 %in% c("Bacteria", "Fungi")) %>%
  pull(value) %>% min()



png("Plots/demonstration_parameteruncertainty.png", width = 12, height = 6, units = "in", res = 600)
results2 %>%
  mutate(Group3 = factor(Group3, levels = c("Bacteria", "Fungi", "Amoebae", "Flagellates", "BactNem", "FungNem", "OmniNem", "PhytoNem", "PredNem", "FungCol", "CryMite", "NonCryMite", "NemMite","PredMite"))) %>%
  mutate(Web = factor(Web, levels = c("GA","UGA","GB","UGB","GC","UGC","Young","Mid","Old","Heathland"))) %>%
  left_join(
    tibble(MS = c("Andres2016", "Holtkamp2011"),
           Manuscript = c("Andrés et al. 2016", "Holtkamp et al. 2011"))
  ) %>%
  ggplot(aes(x = Group3, y = value, color = Manuscript, group = paste0(Group3, Web))) + geom_hline(yintercept = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000), linetype = 2, color = 'grey') + geom_boxplot(outlier.size = 0.5) + theme_classic() + ylab(parse(text = "Carbon~consumption~(kg[C]~ha^-1~yr^-1)")) + scale_y_continuous(trans = "log", breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)) + geom_vline(xintercept = c(2.5, 4.5, 9.5), linetype = 2) + xlab("Trophic Group") + annotate(geom = "text", x = c(1.5, 3.5, 7, 12), y = 20000, label = c("Microbes", "Single-cell Pred.", "Nematodes", "Microarthropods"))
dev.off()
