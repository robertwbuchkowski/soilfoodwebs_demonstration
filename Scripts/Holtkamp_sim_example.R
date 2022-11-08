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

reps_per_web = 5
results = vector("list",2)

# Load in errors, convert SD to VAR
errmes <- as.matrix(read.csv("Data/Holtkamp/biomass_sd.csv", header = T, row.names = 1))^2

# Prepare function to simulate:
Holtkamp_CNsim <- function(COMMin, ddvec, startvec){
  baseline <- CNsim(COMMin, start_mod = startvec, TIMES = 1:50)
  
  ddrun <- CNsim(COMMin, start_mod = startvec,
                 densitydependence = ddvec, TIMES = 1:50)
  
  return(baseline %>% tibble() %>%
           pivot_longer(-Day) %>%
           filter(grepl("_Carbon", name)) %>%
           mutate(run = "DI") %>%
           bind_rows(
             ddrun %>% tibble() %>%
               pivot_longer(-Day) %>%
               filter(grepl("_Carbon", name)) %>%
               mutate(run = "DD")
           ) %>%
           tibble()%>%
           separate(name, into = c("name", NA),sep = "_") %>%
           left_join(
             TLgroups %>%
               filter(TL %in% c("3+", "3")) %>%
               tibble(), by = c("name"="ID")
           )%>%
           mutate(TL = ifelse(is.na(TL), name, TL)) %>%
           select(-name) %>% rename(name = TL) %>%
           group_by(Day, name, run) %>%
           summarize(value = sum(value), .groups = "drop") %>%
           mutate(RunID = paste0(runif(1)*1e8,paste0(sample(letters,10), collapse = "")))
         ) 
}


# Create a grouped data frame:
TLcheddar(Holtkamp2011$Young$imat)

TLgroups = data.frame(ID = Holtkamp2011$Young$prop$ID,
                      TL = c("3+", "3+", "3+", "3+", "3+", "3+",
                             3,3,3,3,3,3,3,
                             2,2,2,2,
                             1,1,1,1))

Holtmod = Holtkamp2011$Young

Holtmod = comtrosp(Holtmod, c("PredMite", "NemMite", "PredColl","PredNem", "OmniNem", "Amoebae"), newname = "TL3+", deleteCOMBOcannibal = T)

Holtmod = comtrosp(Holtmod, c("FungNem", "FungColl", "FungCryJuvMite","FungCryAdultMite", "FungNoncryMite", "BactNem","Flagellates"), newname = "TL3", deleteCOMBOcannibal = T)



# Simulate:
errmes2 <- as.matrix(errmes[,"Young"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Young,
                             errormeasure = errmes2,
                             errortype = "Variance", 
                             replicates = reps_per_web,
                             returnprops = T,
                             returnresults = F)


puGA = lapply(puGA, FUN = Holtkamp_CNsim, ddvec = c(rep(1, 17), rep(0,4)), startvec = c(rep(1, 17), 1, 1.1, 1.1, 1))


puGA = do.call("rbind",puGA)
results[[1]] = cbind(puGA, Web = "Young", Type = "Complete")

# RESULTS 2
errmes2 <- as.matrix(errmes[,"Young"])
colnames(errmes2) = "B"
errmes2 = 
  tibble(ID = unname(Holtmod$prop$ID)) %>%
  left_join(
    errmes2 %>%
  data.frame() %>%
  rownames_to_column('ID') %>%
  tibble() %>%
  left_join(
    TLgroups
  ) %>%
  mutate(TL = ifelse(TL %in% c("3", "3+"), TL, ID)) %>%
  group_by(TL) %>%
  summarize(B = sum(B)) %>%
  mutate(TL = ifelse(TL == "3", "TL3", TL))%>%
  mutate(TL = ifelse(TL == "3+", "TL3+", TL)) %>%
  rename(ID = TL)
  )

errmes2 = matrix(errmes2 %>%pull(B), dimnames = list(errmes2 %>%pull(ID), "B"))


puGA = parameter_uncertainty(usin = Holtmod,
                             errormeasure = errmes2,
                             errortype = "Variance", 
                             replicates = reps_per_web,
                             returnprops = T,
                             returnresults = F)


puGA = lapply(puGA, FUN = Holtkamp_CNsim, ddvec = c(rep(1, 6), rep(0,4)), startvec = c(rep(1, 6), 1, 1.1, 1.1, 1))


puGA = do.call("rbind",puGA)
results[[2]] = cbind(puGA, Web = "Young", Type = "Simplified")

results2 = do.call("rbind", results)

results2 = results2 %>%
  tibble() %>%
  group_by(name, run, Web, Type, RunID) %>%
  summarise(value = sum(value)) %>%
  group_by(name, run, Web, Type)


results2 %>%
  slice_min(order_by = value) %>%
  select(-value) %>%
  mutate(Lvl = "min") %>%
  bind_rows(
    results2 %>%
      slice_max(order_by = value) %>%
      select(-value) %>%
      mutate(Lvl = "max")
  ) %>%
  left_join(
    do.call("rbind", results), by = c("name", "run", "Web", "Type", "RunID")
  ) %>%
  ungroup() %>%
  mutate(name  = ifelse(name == "3", "TL3", name))%>%
  mutate(name  = ifelse(name == "3+", "TL3+", name)) %>%
  filter(grepl("TL3", name)) %>%
  ggplot(aes(x = Day, y = value, group = Lvl, color = run)) + geom_line() + facet_grid(name~Type, scales = "free") + theme_classic()
  
  
  summarize(avg = sum(value), sd = sd(value), .groups = "drop") %>%
  ungroup() %>%
  ggplot(aes(x = Day, y = avg, color = run)) + geom_line() + facet_wrap(.~name, scales ="free")




















errmes2 <- as.matrix(errmes[,"Young"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Young,
                             errormeasure = errmes2,
                             errortype = "Variance", 
                             replicates = reps_per_web,
                             returnprops = T,
                             returnresults = F)


puGA = lapply(puGA, FUN = Holtkamp_CNsim, ddvec = c(rep(1, 17), rep(0,4)), startvec = c(rep(1, 17), 1, 1.1, 1.1, 1))


puGA = do.call("rbind",puGA)
results[[1]] = cbind(puGA, Web = "Young", Type = "Complete")






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
  ggplot(aes(x = name, y = value, color = Web)) + geom_boxplot() + facet_wrap(.~ID, scales = "free") + theme_classic() + ylab("Effect magnitude") + ggtitle("Nematode effects on nitrogen mineralization") + xlab("Effect type")
dev.off()  
