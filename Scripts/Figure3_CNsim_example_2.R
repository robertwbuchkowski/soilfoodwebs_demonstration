# CNsim example script to make Figure 3:
# Author: Robert W. Buchkowski
# Date created: April 30/2023

# If need you to install new version from github:
if(F){
  devtools::install_github("robertwbuchkowski/soilfoodwebs", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual"))
}

# Load in the libraries:
if (!require("pacman")) install.packages("pacman")
p_load(soilfoodwebs,tidyverse)

comm1 = Andres2016$UGA

comm2 = comtrosp(comm1, selected = c("ActiveSOM", "SlowSOM"), newname = "SOM")
comm2 = comtrosp(comm2, selected = c("FungivorousCryptostigmata", "FungivorousProstigmata", "Nematophagousmites", "Predaceousmites"), newname = "Mites")

comm2 = comtrosp(comm2, selected = c("Bacteriophagousnematodes", "Omnivorousnematodes", "Phytophagousnematodes", "Fungivorousnematodes", "Predaceousnematodes"), newname = "Nematodes")

comm2 = comtrosp(comm2, selected = c("Ciliates", "Flagellates", "Amoeba"), newname = "Protists")


comm2

out1 = CNsim(comm1, TIMES = 1:10,start_mod = c(0.95,rep(1, 18), 0.95, 1))

out2 = CNsim(comm2, TIMES = 1:100,start_mod = c(0.95,rep(1, 10)))

out1 %>%
  tibble() %>%
  pivot_longer(-Day) %>%
  separate(name, into = c("State", "Element"), sep ="_") %>%
  filter(Element %in% c("Carbon","Nitrogen")) %>%
  pivot_wider(names_from = Element) %>%
  mutate(CN = Carbon/Nitrogen) %>%
  ggplot(aes(x = Day, y = CN)) + geom_line() +facet_wrap(.~State, scales = "free")

comm1$prop$ID
