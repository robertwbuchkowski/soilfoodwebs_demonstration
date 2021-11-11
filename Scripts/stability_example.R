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

# Calculate Smin for each community
calc_smin(checkcomm(Koltz2018))
calc_smin(checkcomm(Holtkamp2011$Young))
calc_smin(checkcomm(Holtkamp2011$Mid))
calc_smin(checkcomm(Holtkamp2011$Old))
calc_smin(checkcomm(Holtkamp2011$Heathland))
calc_smin(checkcomm(Hunt1987))
calc_smin(checkcomm(Andres2016$GA))
calc_smin(checkcomm(Andres2016$UGA))
calc_smin(checkcomm(Andres2016$GB))
calc_smin(checkcomm(Andres2016$UGB))
calc_smin(checkcomm(Andres2016$GC))
calc_smin(checkcomm(Andres2016$UGC))
calc_smin(checkcomm(deRuiter1994$CON))
calc_smin(checkcomm(deRuiter1994$CON10))
calc_smin(checkcomm(deRuiter1994$INT))
calc_smin(checkcomm(deRuiter1994$INT10))

# See whether detritus recycling is the issue:
temp1 <- checkcomm(Holtkamp2011$Young)
temp1$prop$isDetritus = 0
temp1$prop$DetritusRecycling = 0
temp1$prop$DetritusRecycling[18] = 10
temp1$prop$DetritusRecycling[19] = 0
temp1$prop$DetritusRecycling[20] = 10

calc_smin(checkcomm(temp1))
