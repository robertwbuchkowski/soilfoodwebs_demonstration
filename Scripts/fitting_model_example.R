# Script for the paper presenting soilfoodwebs
# Author: Robert W. Buchkowski
# Date created: Nov. 8/ 2021

# If need you to install new version from github:
if(F){
  devtools::install_github("robertwbuchkowski/soilfoodwebs", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual"))
}

# Load in the libraries:
if (!require("pacman")) install.packages("pacman")
p_load(soilfoodwebs,tidyverse)

# Adjusting microbial production efficiency ----

# Create function to modify fungal production efficiency
ccxf <- function(p, ccx){
  
  ccx$prop$p[14] = p[1]
  ccx$prop$p[15] = p[2]
  
  return(c(Cmin = sum(comana(corrstoich(ccx))$Cmin), Nmin = sum(comana(corrstoich(ccx))$Nmin)))
}

# Return carbon and nitrogen mineralization across these gradients
res1 = expand.grid(1:100/100,1:100/100)
res2 = res1

for(i in 1:dim(res1)[1]){
  res2[i,] = ccxf(as.numeric(res1[i,]), ccx = Hunt1987)
}

res = cbind(res1, res2)
colnames(res) = c("p1", "p2", "Cmin", "Nmin")

# Build range:

buildrange1 = function(X, Y){
  (ccxf(c(X, Y), Hunt1987)[1] - ccxf(c(0.3, 0.3), Hunt1987)[1])^2
}

res3 = data.frame(p2 = 1:100/100,
                  p1 = NA,
                  diff = NA)
for(i in 1:dim(res3)[1]){
  remp1 = optimize(buildrange1, c(0,1), Y = res3[i,"p2"])
  res3[i,"p1"] = remp1$minimum
  res3[i,"diff"] = sqrt(remp1$objective)
}

# Get rid of the cases where there is no solution
res3 = subset(res3, res3$diff < 1)


gg = tibble(res) %>%
  ggplot(aes(x = p1, y = p2)) + geom_contour_filled(aes(z = log(Cmin))) +
  geom_line(data = res3, size = 3) + theme_classic() + xlab("Production efficiency Bacteria") + ylab("Production efficiency Fungi") + scale_fill_viridis_d(name = "Cmin (log)")

ff = tibble(res) %>%
  filter(p1 == round(p1, digits = 1) & 
           p2 == round(p2, digits = 1)) %>%
  ggplot(aes(x = Cmin, y = Nmin, color = p2, fill = p1))+ geom_hline(yintercept = sum(comana(corrstoich(Hunt1987))$Nmin), linetype = 2)+ geom_vline(xintercept = sum(comana(corrstoich(Hunt1987))$Cmin), linetype = 2) + geom_point(shape = 21, size = 2) + theme_classic() + scale_color_viridis_c(name = "Prod. eff.\nFungi (boarder)", option = "C") + scale_fill_gradient(name = "Prod. eff.\nBacteria (fill)", low = "grey80", high = "grey10")

png("Plots/demonstration_fitting.png", width = 10, height = 4, units = "in", res = 600)
ggpubr::ggarrange(
  gg,
  ff,
  labels = "AUTO", ncol = 2
)
dev.off()