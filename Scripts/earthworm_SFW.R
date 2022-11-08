# Earthworm soilfoodwebs example:

# Author: Robert W. Buchkowski
# Date created: Aug 8/ 2022

# If need you to install new version from github:
if(F){
  devtools::install_github("robertwbuchkowski/soilfoodwebs", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual"))
}

# Load in the libraries:
if (!require("pacman")) install.packages("pacman")
p_load(soilfoodwebs,tidyverse)

# Create the solution:

comm1 = list(imat = matrix(c(0,0,1000,1,
                     0,0,0,1,
                     0,0,0,0,
                     0,0,0,0),
                   nrow = 4, ncol = 4, byrow = T, dimnames = list(c("Epi","Endo", "L", "S"),c("Epi","Endo","L", "S"))),
     
     props = data.frame(ID = c("Epi","Endo", "L", "S"),
                        a = c(0.01, 0.01, 1,1),
                        p = c(0.17, 0.35, 1,1),
                        d = c(0.014, 0.014, 0,0),
                        CN = c(5,5, 30, 10),
                        B = c(1.4,2, 9900, 77800),
                        canIMM = c(0,0,0,0),
                        isDetritus = c(0,0,1,1),
                        DetritusRecycling = c(0,0,1,0),
                        isPlant = c(0,0,0,0))
     
     )

comana(comm1)

comana(corrstoich(comm1, dietlimits = matrix(c(1,1,1,0.1,
                                               1,1,1,1,
                                               1,1,1,1,
                                               1,1,1,1),
                                             nrow = 4, ncol = 4, byrow = T, dimnames = list(c("Epi","Endo", "L", "S"),c("Epi","Endo","L", "S")))))
