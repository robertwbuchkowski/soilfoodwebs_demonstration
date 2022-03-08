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

# Create a series of communities of different sizes to run ----

# Base these communities of the Koltz web--the larget one in the package

curcomm = Holtkamp2011$Heathland
smod = c(rep(1, 17),1, 1, 0.9,1)

dim(curcomm$imat)[1]

results = vector(mode = "list", length = dim(curcomm$imat)[1])

tt = system.time(CNsim(curcomm, start_mod = smod))

results[[1]] = c(Size = dim(curcomm$imat)[1], Time = unname(tt["elapsed"]))

IDS = curcomm$prop$ID

curcomm = renamenode(curcomm,"PredMite", "GrNode")

for(i in 2:17){
  curcomm = comtrosp(curcomm, selected = c("GrNode", IDS[i]), newname = "GrNode")
  smod = smod[-1]
  # p1 = CNsim(curcomm, start_mod = smod) %>%
  #   select(contains("Day")| contains("_Carbon")) %>%
  #   pivot_longer(-Day) %>%
  #   ggplot(aes(x = Day, y = value)) + geom_line() + facet_wrap(.~name, scales = "free")
  # plot(p1)
  tt = system.time(CNsim(curcomm, start_mod = smod))
  results[[i]] = c(Size = dim(curcomm$imat)[1], Time = unname(tt["elapsed"]))
  print(results[[i]])
  print(curcomm$prop$ID)
}

results2 = do.call("rbind",results)

tibble(as.data.frame(results2)) %>%
  ggplot(aes(x = Size, y = Time)) + geom_line() + theme_classic()

# CPER food web ---

curcomm = Hunt1987
smod = c(rep(1, 15),0.9,1)

dim(curcomm$imat)[1]

Aresults = vector(mode = "list", length = dim(curcomm$imat)[1])

tt = system.time(CNsim(curcomm, start_mod = smod))

Aresults[[1]] = c(Size = dim(curcomm$imat)[1], Time = unname(tt["elapsed"]))

IDS = curcomm$prop$ID

curcomm = renamenode(curcomm,"Predatorymites", "GrNode")

for(i in 2:15){
  curcomm = comtrosp(curcomm, selected = c("GrNode", IDS[i]), newname = "GrNode")
  smod = smod[-1]
  # p1 = CNsim(curcomm, start_mod = smod) %>%
  #   select(contains("Day")| contains("_Carbon")) %>%
  #   pivot_longer(-Day) %>%
  #   ggplot(aes(x = Day, y = value)) + geom_line() + facet_wrap(.~name, scales = "free")
  # plot(p1)
  tt = system.time(CNsim(curcomm, start_mod = smod))
  Aresults[[i]] = c(Size = dim(curcomm$imat)[1], Time = unname(tt["elapsed"]))
  print(Aresults[[i]])
  # print(curcomm$prop$ID)
}

Aresults2 = do.call("rbind",Aresults)

tibble(as.data.frame(Aresults2)) %>%
  ggplot(aes(x = Size, y = Time)) + geom_line() + theme_classic()

png("Plots/sim_time.png", width = 5, height = 3, units = "in", res = 600)
tibble(as.data.frame(results2))  %>%
  mutate(Web = "Holtkamp et al. 2011") %>%
  bind_rows(
    tibble(as.data.frame(Aresults2)) %>%
      mutate(Web = "Hunt et al. 1987")
  ) %>%
  ggplot(aes(x = Size, y = Time, color = Web)) + geom_line() + theme_classic() + scale_y_log10(name = "Seconds")
dev.off()
