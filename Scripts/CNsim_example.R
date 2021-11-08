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

# CNsim and the effects of density-dependence ----

# Use the CPER food web as an example:

# Create a list of 100 communities for CPER:
several_comms <- parameter_uncertainty(CPER, returnprops = T) # Return the communities for further analyses.

# Select only the communities:
several_comms = several_comms$communitylist

# Run the communities through the CNsim function using lapply:

CPER_CNsim <- function(COMMin){
  baseline <- CNsim(COMMin, start_mod = c(rep(1, 15), 1.1, 1))
  
  ddrun <- CNsim(COMMin, start_mod = c(rep(1, 15), 1.1, 1),
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

# Combine by trophic level:
groupdf = data.frame(TLcheddar(CPER$imat))

groupdf = cbind(groupdf,ID = rownames(groupdf))
rownames(groupdf)= NULL

colnames(groupdf) = c("TL", "name")

groupdf$TL = ifelse(round(groupdf$TL, 1) > 2.5, "3+", round(groupdf$TL, 1))

CPERmod = comtrosp(CPER, selected = groupdf$name[groupdf$TL == 2.5])

CPERmod = comtrosp(CPERmod, selected = groupdf$name[groupdf$TL == "3+"])

CPERmod = comtrosp(CPERmod, selected = c(groupdf$name[groupdf$TL == 2][2:3],groupdf$name[groupdf$TL == 1][1]))

comana(CPERmod, mkplot = T, whattoplot = "web")

groupsCPER = tibble(ID = groupdf$name[groupdf$TL == 2.5], TL = "TL = 2.5") %>% bind_rows(
  tibble(ID = groupdf$name[groupdf$TL == "3+"], TL = "TL = 3+")
) %>%
  bind_rows(
    tibble(ID = c(groupdf$name[groupdf$TL == 2][2:3],groupdf$name[groupdf$TL == 1][1]), TL = "Microbes")
  )

# Create a list of 100 communities for CPER:
several_comms2 <- parameter_uncertainty(CPERmod, returnprops = T) # Return the communities for further analyses.

# Select only the communities:
several_comms2 = several_comms2$communitylist


CPERmod_CNsim <- function(COMMin){
  baseline <- CNsim(COMMin, start_mod = c(1,1,1,1,1.1,1))
  
  ddrun <- CNsim(COMMin, start_mod = c(1,1,1,1,1.1,1),
                 densitydependence = c(1,1,1,1,0,0))
  
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

# Run the versions through the simulator

CPERres <- lapply(several_comms[1:100],CPER_CNsim)
# dim(CPERres[[1]]) # 3400
CPERres <- do.call("rbind",CPERres)
CPERres[,"RunID"] = rep(1:100, each = 3400)
write_rds(CPERres,"Data/soilfoodwebs_demonstration/CPERres.rds")

CPERmodres <- lapply(several_comms2[1:100],CPERmod_CNsim)
# dim(CPERmodres[[1]]) # 1200
CPERmodres <- do.call("rbind",CPERmodres)
CPERmodres[,"RunID"] = rep(1:100, each = 1200)
write_rds(CPERmodres,"Data/soilfoodwebs_demonstration/CPERmodres.rds")

# Read back in the data:
CPERres <- read_rds("Data/soilfoodwebs_demonstration/CPERres.rds")
CPERmodres <- read_rds("Data/soilfoodwebs_demonstration/CPERmodres.rds")

# Create variation plots:

CPERfull <- CPERres %>%
  separate(name, into = c('ID', NA), sep = "_") %>%
  left_join(
    groupsCPER
  ) %>%
  mutate(TL = ifelse(is.na(TL), ID, TL)) %>%
  group_by(Day, run, TL, RunID) %>%
  summarize(value = sum(value)) %>%
  mutate(Model = "Complete") %>%
  bind_rows(
    CPERmodres %>%
      left_join(
        tibble(name = unique(CPERmodres$name),
               TL = c("TL = 3+","TL = 2.5", "Phytophagousnematodes", "Microbes", "Detritus", "Roots")), by = "name"
      ) %>%
      select(-name) %>%
      mutate(Model = "Simple")
  )

CPERexample <- CPER_CNsim(CPER) %>%
  separate(name, into = c('ID', NA), sep = "_") %>%
  left_join(
    groupsCPER
  ) %>%
  mutate(TL = ifelse(is.na(TL), ID, TL)) %>%
  group_by(Day, run, TL) %>%
  summarize(value = sum(value)) %>%
  mutate(Model = "Complete") %>%
  bind_rows(
    CPERmod_CNsim(CPERmod) %>%
      left_join(
        tibble(name = unique(CPERmodres$name),
               TL = c("TL = 3+","TL = 2.5", "Phytophagousnematodes", "Microbes", "Detritus", "Roots")), by = "name"
      ) %>%
      select(-name) %>%
      mutate(Model = "Simple")
  )

png("Plots/demonstration_simulation.png", width = 8, height = 5, units = "in", res = 600)
CPERfull %>%
  group_by(run, TL, RunID, Model) %>%
  summarise(value = sum(value)) %>%
  group_by(run, TL, Model) %>%
  slice_min(order_by = value) %>%
  mutate(Type = "Min") %>%
  bind_rows(
    CPERfull %>%
      group_by(run, TL, RunID, Model) %>%
      summarise(value = sum(value)) %>%
      group_by(run, TL, Model) %>%
      slice_max(order_by = value)%>%
      mutate(Type = "Max")
  ) %>%
  select(-value) %>%
  ungroup() %>%
  left_join(
    CPERfull, by = c("run", "TL", "RunID", "Model")
  )%>%
  filter(!(TL %in% c("Phytophagousnematodes", "Roots"))) %>%
  select(-RunID) %>%
  pivot_wider(names_from = Type) %>%
  ggplot(aes(x = Day, fill = run, group = paste0(run, Model))) + geom_ribbon(aes(ymin = Min, ymax = Max), alpha = 0.4) + facet_grid(TL~Model, scales = "free") + theme_classic() + scale_fill_manual(values = c("blue", "orange")) + scale_color_manual(values = c("blue", "orange")) +
  geom_line(aes(y = value, color = run),
            data = CPERexample %>%
              filter(!(TL %in% c("Phytophagousnematodes", "Roots"))) ) + ylab("Biomass (Kg[C] per ha)") + xlab("Year")
dev.off()  

# Explore the oscillations by individual species:
CPER_explore_ind <- CPER_CNsim(CPER) %>%
  separate(name, into = c('ID', NA), sep = "_")

CPER_explore_ind %>%
  filter(run == "base") %>%
  left_join(
    groupsCPER
  ) %>%
  left_join(
    CPER_explore_ind %>% filter(Day == 1) %>% select(-Day) %>% rename(base = value)
  ) %>%
  mutate(TL = ifelse(is.na(TL), ID, TL)) %>%
  ggplot(aes(x = Day, y = value/base, color = ID)) + geom_line(size = 1.4) + facet_wrap(.~TL, scales = "free_y") + theme(legend.position = "none")
