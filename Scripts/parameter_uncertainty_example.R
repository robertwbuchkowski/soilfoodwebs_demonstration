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

pulldata <- function(x){
  tempcomb = comtrosp(x$usin, selected = c("FungivorousCryptostigmata", "FungivorousProstigmata", "Nematophagousmites", "Predaceousmites"), newname = "Mite", allFEEDING1 = T,deleteCOMBOcannibal = T)
  
  t1 = whomineralizes(tempcomb, selected = "Mite")
  names(t1) = paste0("Base_", names(t1))
  
  t1a = decompexpt(x$usin, selected = "Proturans")$basedecomp
  names(t1a) = c("Base_Fast", "Base_Slow")
  t1a = data.frame(t(t1a))
  
  t1b = decompexpt(tempcomb, selected = "Mite")$decompeffects %>% pivot_longer(!ID & !DetritusID) %>% mutate(DetritusID = ifelse(DetritusID == "ActiveSOM", "Fast", "Slow")) %>% mutate(name = paste0("Base_",DetritusID, "/", name)) %>% select(-DetritusID) %>% pivot_wider() %>% as.data.frame()
  
  tempcomb = comtrosp(tempcomb, selected = c("Bacteriophagousnematodes","Fungivorousnematodes","Omnivorousnematodes", "Phytophagousnematodes", "Predaceousnematodes"), newname = "Nematode", allFEEDING1 = T,deleteCOMBOcannibal = T)
  
  tempcomb = comtrosp(tempcomb, selected = c("Ciliates", "Flagellates", "Amoeba"), newname = "Protist", allFEEDING1 = T,deleteCOMBOcannibal = T)
  
  t2 = whomineralizes(tempcomb, selected = "Mite")
  names(t2) = paste0("Mod_", names(t2))
  
  t2a = decompexpt(tempcomb, selected = "Proturans")$basedecomp
  names(t2a) = c("Mod_Fast", "Mod_Slow")
  t2a = data.frame(t(t2a))
  
  t2b = decompexpt(tempcomb, selected = "Mite")$decompeffects %>% pivot_longer(!ID & !DetritusID) %>% mutate(DetritusID = ifelse(DetritusID == "ActiveSOM", "Fast", "Slow")) %>% mutate(name = paste0("Mod_",DetritusID, "/", name)) %>% select(-DetritusID) %>% pivot_wider() %>% as.data.frame()
  
  t3 = comana(tempcomb)

  cbind(Base_Cmin = sum(x$Cmin), Base_Nmin = sum(x$Nmin),Mod_Cmin = sum(t3$Cmin), Mod_Nmin = sum(t3$Nmin),t1[,-1],t2[,-1], t1a, t1b[,-1], t2a, t2b[,-1])
}

CNsim(Andres2016$GB, start_mod = c(rep(1, 19), 1.1, 1), TIMES = 1:10)

CPER_CNsim <- function(COMMin){
  baseline <- CNsim(COMMin, start_mod = c(rep(1, 19), 1.1, 1))
  
  ddrun <- CNsim(COMMin, start_mod = c(rep(1, 20), 1.1, 1),
                 densitydependence = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0))
  
  return(baseline %>% tibble() %>%
           pivot_longer(-Day) %>%
           filter(grepl("_Carbon", name)) %>%
           mutate(run = "base") %>%
           bind_rows(
             ddrun %>% tibble() %>%
               pivot_longer(-Day) %>%
               filter(grepl("_Carbon", name)) %>%
               mutate(run = "Density-dependent")
           )) 
}

puGA = do.call("rbind",lapply(puGA, pulldata))

results[[1]] = cbind(as.data.frame(puGA), Web = "GA")


biomassescur = biomasses %>% filter(name == "A_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGA,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pulldata))
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


puGA = do.call("rbind",lapply(puGA, pulldata))
results[[3]] = cbind(as.data.frame(puGA), Web = "GB")


biomassescur = biomasses %>% filter(name == "B_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGB,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pulldata))
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


puGA = do.call("rbind",lapply(puGA, pulldata))
results[[5]] = cbind(as.data.frame(puGA), Web = "GC")


biomassescur = biomasses %>% filter(name == "C_UG") %>%
  select(-name) %>% data.frame()

puGA = parameter_uncertainty(usin = Andres2016$UGC,
                             errormeasure = 
                               matrix(data = biomassescur$B, nrow = dim(biomassescur)[1], ncol = 1, dimnames = list(biomassescur$ID, "B")),
                             
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)


puGA = do.call("rbind",lapply(puGA, pulldata))
results[[6]] = cbind(as.data.frame(puGA), Web = "UGC")

# Holtkamp ------

pulldata <- function(x){
  tempcomb = comtrosp(x$usin, selected = c("FungCryJuvMite","FungCryAdultMite", "FungNoncryMite", "NemMite", "PredMite"), newname = "Mite", allFEEDING1 = T,deleteCOMBOcannibal = T)
  
  t1 = whomineralizes(tempcomb, selected = "Mite")
  names(t1) = paste0("Base_", names(t1))
  
  t1a = decompexpt(x$usin, selected = "Amoebae")$basedecomp[-1] # Get rid of the sugars
  names(t1a) = c("Base_Fast", "Base_Slow")
  t1a = data.frame(t(t1a))
  
  t1b = decompexpt(tempcomb, selected = "Mite")$decompeffects[-1,] %>% pivot_longer(!ID & !DetritusID) %>% mutate(DetritusID = ifelse(DetritusID == "LabileOM", "Fast", "Slow")) %>% mutate(name = paste0("Base_",DetritusID, "/", name)) %>% select(-DetritusID) %>% pivot_wider() %>% as.data.frame()
  
  tempcomb = comtrosp(tempcomb, selected = c("BactNem","FungNem","OmniNem", "PhytoNem", "PredNem"), newname = "Nematode", allFEEDING1 = T,deleteCOMBOcannibal = T)
  
  tempcomb = comtrosp(tempcomb, selected = c("PredColl", "FungColl"), newname = "Collembola", allFEEDING1 = T,deleteCOMBOcannibal = T)
  
  tempcomb = comtrosp(tempcomb, selected = c("Flagellates", "Amoebae"), newname = "Protist", allFEEDING1 = T,deleteCOMBOcannibal = T)
  
  t2 = whomineralizes(tempcomb, selected = "Mite")
  names(t2) = paste0("Mod_", names(t2))
  
  t2a = decompexpt(tempcomb, selected = "Protist")$basedecomp[-1]# Get rid of the sugars
  names(t2a) = c("Mod_Fast", "Mod_Slow")
  t2a = data.frame(t(t2a))
  
  t2b = decompexpt(tempcomb, selected = "Mite")$decompeffects[-1,] %>% pivot_longer(!ID & !DetritusID) %>% mutate(DetritusID = ifelse(DetritusID == "LabileOM", "Fast", "Slow")) %>% mutate(name = paste0("Mod_",DetritusID, "/", name)) %>% select(-DetritusID) %>% pivot_wider() %>% as.data.frame()
  
  t3 = comana(tempcomb)
  
  cbind(Base_Cmin = sum(x$Cmin), Base_Nmin = sum(x$Nmin),Mod_Cmin = sum(t3$Cmin), Mod_Nmin = sum(t3$Nmin),t1[,-1],t2[,-1], t1a, t1b[,-1], t2a, t2b[,-1])
}


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

puGA = do.call("rbind",lapply(puGA, pulldata))
results[[7]] = cbind(as.data.frame(puGA), Web = "Young")

errmes2 <- as.matrix(errmes[,"Mid"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Mid,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

puGA = do.call("rbind",lapply(puGA, pulldata))
results[[8]] = cbind(as.data.frame(puGA), Web = "Mid")

errmes2 <- as.matrix(errmes[,"Health"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Heathland,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

puGA = do.call("rbind",lapply(puGA, pulldata))
results[[10]] = cbind(as.data.frame(puGA), Web = "Heathland")

pulldata <- function(x){
  tempcomb = comtrosp(x$usin, selected = c("FungCryJuvMite","FungCryAdultMite", "FungNoncryMite", "NemMite", "PredMite"), newname = "Mite", allFEEDING1 = T,deleteCOMBOcannibal = T)
  
  t1 = whomineralizes(tempcomb, selected = "Mite")
  names(t1) = paste0("Base_", names(t1))
  
  t1a = decompexpt(x$usin, selected = "Amoebae")$basedecomp[-1] # Get rid of the sugars
  names(t1a) = c("Base_Fast", "Base_Slow")
  t1a = data.frame(t(t1a))
  
  t1b = decompexpt(tempcomb, selected = "Mite")$decompeffects[-1,] %>% pivot_longer(!ID & !DetritusID) %>% mutate(DetritusID = ifelse(DetritusID == "LabileOM", "Fast", "Slow")) %>% mutate(name = paste0("Base_",DetritusID, "/", name)) %>% select(-DetritusID) %>% pivot_wider() %>% as.data.frame()
  
  tempcomb = comtrosp(tempcomb, selected = c("BactNem","FungNem","OmniNem", "PhytoNem", "PredNem"), newname = "Nematode", allFEEDING1 = T,deleteCOMBOcannibal = T)
  
  # No need to to Collembola because there are no predators at this site.
  
  tempcomb = comtrosp(tempcomb, selected = c("Flagellates", "Amoebae"), newname = "Protist", allFEEDING1 = T,deleteCOMBOcannibal = T)
  
  t2 = whomineralizes(tempcomb, selected = "Mite")
  names(t2) = paste0("Mod_", names(t2))
  
  t2a = decompexpt(tempcomb, selected = "Protist")$basedecomp[-1]# Get rid of the sugars
  names(t2a) = c("Mod_Fast", "Mod_Slow")
  t2a = data.frame(t(t2a))
  
  t2b = decompexpt(tempcomb, selected = "Mite")$decompeffects[-1,] %>% pivot_longer(!ID & !DetritusID) %>% mutate(DetritusID = ifelse(DetritusID == "LabileOM", "Fast", "Slow")) %>% mutate(name = paste0("Mod_",DetritusID, "/", name)) %>% select(-DetritusID) %>% pivot_wider() %>% as.data.frame()
  
  t3 = comana(tempcomb)
  
  cbind(Base_Cmin = sum(x$Cmin), Base_Nmin = sum(x$Nmin),Mod_Cmin = sum(t3$Cmin), Mod_Nmin = sum(t3$Nmin),t1[,-1],t2[,-1], t1a, t1b[,-1], t2a, t2b[,-1])
}

errmes2 <- as.matrix(errmes[-3,"Old"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Old,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             fcntorun = "comana", 
                             replicates = reps_per_web)

puGA = do.call("rbind",lapply(puGA, pulldata))
results[[9]] = cbind(as.data.frame(puGA), Web = "Old")

# Summarize

results2 = do.call("rbind", results)

results2 = results2 %>%
  mutate(UID = 1:dim(.)[1]) %>%
  pivot_longer(!Web & !UID) %>%
  separate(name, into = c("Structure", "name"), sep = "_") %>%
  mutate(Manuscript = ifelse(Web %in% c("Young", "Mid", "Old", "Heathland"), "Holtkamp et al.", "Andres et al.")) %>%
  mutate(Web = factor(Web, levels = c("Young", "Mid", "Old", "Heathland", "GA", "UGA", "GB", "UGB", "GC", "UGC")))

png("Plots/demonstration_parameteruncertainty.png", width = 10, height = 11, units = "in", res = 600)
ggpubr::ggarrange(
  results2 %>%
    filter(name %in% c("Cmin", "Nmin")) %>%
    mutate(name = ifelse(name == "Cmin", "Carbon", "Nitrogen")) %>%
    mutate(Structure = ifelse(Structure == "Base", "Original", "Grouped")) %>%
    mutate(Structure = factor(Structure, levels = c("Original", "Grouped"))) %>%
    ggplot(aes(x = Web, y = value, color = Structure)) + geom_boxplot() + theme_classic() + facet_wrap(.~name, scales = "free") +  ylab(parse(text = "Mineralization~(kg[Nutrient]~ha^-1~yr^-1)")) + geom_vline(xintercept = 4.5, linetype = 2, color = "grey"),
  
  results2 %>%
    filter(name %in% c("Fast", "Slow")) %>%
    mutate(name = ifelse(name == "Fast", "Fast litter decomposition", "Slow litter decomposition")) %>%
    mutate(Structure = ifelse(Structure == "Base", "Original", "Grouped")) %>%
    mutate(Structure = factor(Structure, levels = c("Original", "Grouped"))) %>%
    ggplot(aes(x = Web, y = value, color = Structure)) + geom_boxplot() + theme_classic() + facet_wrap(.~name, scales = "free") +  ylab(parse(text = "Decomposition~(yr^-1)")) + geom_vline(xintercept = 4.5, linetype = 2, color = "grey"),
  
  results2 %>%
    filter(name %in% c("DirectC", "IndirectC")) %>%
    mutate(name = ifelse(name == "DirectC", "Mite direct effect on carbon", "Mite indirect effect on carbon")) %>%
    mutate(Structure = ifelse(Structure == "Base", "Original", "Grouped")) %>%
    mutate(Structure = factor(Structure, levels = c("Original", "Grouped"))) %>%
    ggplot(aes(x = Web, y = value, color = Structure)) + geom_boxplot() + theme_classic() + facet_wrap(.~name, scales = "free") +  ylab(parse(text = "Effect~(proportion)")) + geom_vline(xintercept = 4.5, linetype = 2, color = "grey"),
  
  results2 %>%
    filter(name %in% c("Fast/Indirect", "Slow/Indirect")) %>%
    mutate(name = ifelse(name == "Fast/Indirect", "Mite indirect effect on fast litter decomposition", "Mite indirect effect on slow litter decomposition")) %>%
    mutate(Structure = ifelse(Structure == "Base", "Original", "Grouped")) %>%
    mutate(Structure = factor(Structure, levels = c("Original", "Grouped"))) %>%
    ggplot(aes(x = Web, y = value, color = Structure)) + geom_boxplot() + theme_classic() + facet_wrap(.~name, scales = "free") +  ylab(parse(text = "Effect~(proportion)")) + geom_vline(xintercept = 4.5, linetype = 2, color = "grey"),
  ncol = 1,common.legend = T
)
dev.off()
