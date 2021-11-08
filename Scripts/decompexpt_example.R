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

# Decomposition experiments ---- 

results = vector("list", 4)

# Load in errors, convert SD to VAR
errmes <- as.matrix(read.csv("Data/Holtkamp/biomass_sd.csv", header = T, row.names = 1))^2

errmes2 <- as.matrix(errmes[,"Young"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Young,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             returnprops = T,
                             returnresults = F)

for(i in 1: length(puGA)){
  temp =  decompexpt(corrstoich(puGA[[i]]), overtime = 50)$overtime$LabileOM %>%
    select(Day, Original) %>%
    full_join(
      decompexpt(corrstoich(removenodes(puGA[[i]], c("FungCryJuvMite", "FungCryAdultMite"))), overtime = 50)$overtime$LabileOM %>%
        select(Day, Original) %>%
        rename(Nomite = Original), by = "Day"
    )  %>%
    mutate(RunID = i)
  puGA[[i]] = temp
}

puGA = do.call("rbind",puGA)
results[[1]] = cbind(puGA, Web = "Young")

errmes2 <- as.matrix(errmes[,"Mid"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Mid,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             returnprops = T,
                             returnresults = F
)

for(i in 1: length(puGA)){
  temp =  decompexpt(corrstoich(puGA[[i]]), overtime = 50)$overtime$LabileOM %>%
    select(Day, Original) %>%
    full_join(
      decompexpt(corrstoich(removenodes(puGA[[i]], c("FungCryJuvMite", "FungCryAdultMite"))), overtime = 50)$overtime$LabileOM %>%
        select(Day, Original) %>%
        rename(Nomite = Original), by = "Day"
    )  %>%
    mutate(RunID = i)
  puGA[[i]] = temp
}

puGA = do.call("rbind",puGA)
results[[2]] = cbind(puGA, Web = "Mid")

errmes2 <- as.matrix(errmes[-3,"Old"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Old,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             returnprops = T,
                             returnresults = F
)

for(i in 1: length(puGA)){
  temp =  decompexpt(corrstoich(puGA[[i]]), overtime = 50)$overtime$LabileOM %>%
    select(Day, Original) %>%
    full_join(
      decompexpt(corrstoich(removenodes(puGA[[i]], c("FungCryJuvMite", "FungCryAdultMite"))), overtime = 50)$overtime$LabileOM %>%
        select(Day, Original) %>%
        rename(Nomite = Original), by = "Day"
    )  %>%
    mutate(RunID = i)
  puGA[[i]] = temp
}

puGA = do.call("rbind",puGA)
results[[3]] = cbind(puGA, Web = "Old")

errmes2 <- as.matrix(errmes[,"Health"])
colnames(errmes2) = "B"

puGA = parameter_uncertainty(usin = Holtkamp2011$Heathland,
                             errormeasure = errmes2,
                             errortype = "Variance",
                             returnprops = T,
                             returnresults = F
)

for(i in 1: length(puGA)){
  temp =  decompexpt(corrstoich(puGA[[i]]), overtime = 50)$overtime$LabileOM %>%
    select(Day, Original) %>%
    full_join(
      decompexpt(corrstoich(removenodes(puGA[[i]], c("FungCryJuvMite", "FungCryAdultMite"))), overtime = 50)$overtime$LabileOM %>%
        select(Day, Original) %>%
        rename(Nomite = Original), by = "Day"
    )  %>%
    mutate(RunID = i)
  puGA[[i]] = temp
}
puGA = do.call("rbind",puGA)
results[[4]] = cbind(puGA, Web = "Heathland")

results2 = do.call("rbind",results) %>%
  tibble() %>%
  group_by(Day, Web) %>%
  summarize(m = mean(100*Original), sd = sd(100*Original)) 

Figure3A = results2 %>%
  ggplot(aes(x = Day, y = m, color = Web, fill = Web, group = Web)) + geom_line(lwd = 1.5) + theme_classic() + xlab("Year") + ylab("Detritus Remaining (%)") +
  geom_errorbar(aes(ymin = m -sd, ymax = m+sd),data = results2 %>%
                  filter(Day == 10 & Web == "Young"|
                           Day == 11 & Web == "Mid"|
                           Day == 12 & Web == "Old"|
                           Day == 13 & Web == "Heathland"|
                           Day == 25 & Web == "Young"|
                           Day == 26 & Web == "Mid"|
                           Day == 27 & Web == "Old"|
                           Day == 28 & Web == "Heathland"|
                           Day == 40 & Web == "Young"|
                           Day == 41 & Web == "Mid"|
                           Day == 42 & Web == "Old"|
                           Day == 43 & Web == "Heathland") )

# Look at figure
Figure3A

# Example for effect of oribatids

results2 = do.call("rbind",results) %>%
  tibble() %>%
  mutate(EOO = 100*(Original - Nomite)/Original) %>%
  group_by(Day, Web) %>%
  summarize(m = mean(EOO), sd = sd(EOO)) 

Figure3B = results2 %>%
  ggplot(aes(x = Day, y = m, color = Web, fill = Web, group = Web)) + geom_line(lwd = 1.5) + theme_classic() + xlab("Year") + ylab("Effect of oribatids on remaining litter (%)") +
  geom_errorbar(aes(ymin = m -sd, ymax = m+sd),data = results2 %>%
                  filter(Day == 10 & Web == "Young"|
                           Day == 11 & Web == "Mid"|
                           Day == 12 & Web == "Old"|
                           Day == 13 & Web == "Heathland"|
                           Day == 25 & Web == "Young"|
                           Day == 26 & Web == "Mid"|
                           Day == 27 & Web == "Old"|
                           Day == 28 & Web == "Heathland"|
                           Day == 40 & Web == "Young"|
                           Day == 41 & Web == "Mid"|
                           Day == 42 & Web == "Old"|
                           Day == 43 & Web == "Heathland") )

png("Plots/demonstration_decomp.png", width = 8, height = 4, units = "in", res = 600)
ggpubr::ggarrange(Figure3A,Figure3B, labels = "AUTO", common.legend = T)
dev.off()
# Get the biomass and relative biomass of oribatid mites in the different Holtkamp webs:

mite_prop <- function(COMM){
  COMM$prop %>%
    select(ID, B) %>%
    mutate(Mite = ID %in% c("FungCryJuvMite", "FungCryAdultMite", "FungNoncryMite")) %>%
    group_by(Mite) %>%
    summarise(B = sum(B))
}

do.call("rbind",lapply(Holtkamp2011, mite_prop))
