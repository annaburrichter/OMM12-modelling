install.packages("remotes")
remotes::install_github("SysBioChalmers/sybil")
install.packages(c("ReacTran", "deSolve", "glpkAPI", "R.matlab"))
install.packages(c('RcppArmadillo', 'RcppEigen'))
install.packages("BacArena-1.8.1.tar.gz", repos = NULL, type = "source")

library(tidyverse)
library(sybil)
library(BacArena)
library(igraph)
library(data.table)

#read medium
medium_aam <- fread("AAM_nomuc_rev5.csv")
medium_m10 <- fread("m10_vj_tre.csv")


#extract models as list, keep only relevant ones
csm_aam <- `OMM_munich-csm-aam_20241205`[1:120,1]

#extract individual models, keep only _all, not _strict
for (i in 1:length(csm_aam)) {
a <- csm_aam[[i]]
name <- paste(a@mod_desc)
b <- Bac(a, typename = name)
if (grepl("_all", name) == TRUE)
  assign(name, b)
else
  rm(b)
}
rm(a,i, name)

#create arenas

arena_aam <- Arena(n=50,m=50)
arena_aam <- addOrg(arena_aam, `KB1-AAM_all`, amount=26)
#arena_aam <- addOrg(arena_aam, `YL2-AAM_all`, amount=0)
#arena_aam <- addOrg(arena_aam, `KB18-AAM_all`, amount=0)
arena_aam <- addOrg(arena_aam, `YL27-AAM_all`, amount=14)
arena_aam <- addOrg(arena_aam, `YL31-AAM_all`, amount=24)
arena_aam <- addOrg(arena_aam, `YL32-AAM_all`, amount=2)
arena_aam <- addOrg(arena_aam, `YL44-AAM_all`, amount=9)
arena_aam <- addOrg(arena_aam, `YL45-AAM_all`, amount=14)
arena_aam <- addOrg(arena_aam, `I46-AAM_all`, amount=6)
arena_aam <- addOrg(arena_aam, `I48-AAM_all`, amount=26)
#arena_aam <- addOrg(arena_aam, `I49-AAM_all`, amount=0)
arena_aam <- addOrg(arena_aam, `YL58-AAM_all`, amount=134)
arena_aam <- addSubs(arena_aam, mediac=paste0("EX_",medium_aam$compounds,"_e0"), smax=medium_aam$maxFlux, unit="mM")

arena_m10 <- Arena(n=50,m=50)
arena_m10 <- addOrg(arena_m10, `KB1-M10_all`, amount=11)
#arena_m10 <- addOrg(arena_m10, `YL2-M10_all`, amount=0)
#arena_m10 <- addOrg(arena_m10, `KB18-M10_all`, amount=0)
#arena_m10 <- addOrg(arena_m10, `YL27-M10_all`, amount=0)
arena_m10 <- addOrg(arena_m10, `YL31-M10_all`, amount=15)
arena_m10 <- addOrg(arena_m10, `YL32-M10_all`, amount=3)
arena_m10 <- addOrg(arena_m10, `YL44-M10_all`, amount=6)
arena_m10 <- addOrg(arena_m10, `YL45-M10_all`, amount=10)
arena_m10 <- addOrg(arena_m10, `I46-M10_all`, amount=3)
arena_m10 <- addOrg(arena_m10, `I48-M10_all`, amount=132)
#arena_m10 <- addOrg(arena_m10, `I49-M10_all`, amount=0)
arena_m10 <- addOrg(arena_m10, `YL58-M10_all`, amount=75)
arena_m10 <- addSubs(arena_m10, mediac=paste0("EX_",medium_m10$compounds,"_e0"), smax=medium_m10$maxFlux, unit="mM")

arena_mouse <- Arena(n=50,m=50)
arena_mouse <- addOrg(arena_mouse, `KB1-mouse_all`, amount=1)
#arena_mouse <- addOrg(arena_mouse, `YL2-mouse_all`, amount=0)
#arena_mouse <- addOrg(arena_mouse, `KB18-mouse_all`, amount=0)
arena_mouse <- addOrg(arena_mouse, `YL27-mouse_all`, amount=10)
arena_mouse <- addOrg(arena_mouse, `YL31-mouse_all`, amount=3)
arena_mouse <- addOrg(arena_mouse, `YL32-mouse_all`, amount=12)
arena_mouse <- addOrg(arena_mouse, `YL44-mouse_all`, amount=4)
arena_mouse <- addOrg(arena_mouse, `YL45-mouse_all`, amount=6)
arena_mouse <- addOrg(arena_mouse, `I46-mouse_all`, amount=3)
arena_mouse <- addOrg(arena_mouse, `I48-mouse_all`, amount=89)
arena_mouse <- addOrg(arena_mouse, `I49-mouse_all`, amount=1)
arena_mouse <- addOrg(arena_mouse, `YL58-mouse_all`, amount=126)
arena_mouse <- addSubs(arena_mouse, mediac=paste0("EX_",medium_m10$compounds,"_e0"), smax=medium_m10$maxFlux, unit="mM")

#simulate 3 steps
sim_OMM12_aam <- simEnv(arena_aam, time=3, sec_obj='mtf')
sim_OMM12_m10 <- simEnv(arena_m10, time=3, sec_obj='mtf')
sim_OMM12_mouse <- simEnv(arena_mouse, time=3, sec_obj='mtf')

`OMM_munich-csm-aam_20241114` <- readRDS("OMM_munich-csm-aam_20241114.RDS")
csm_alt <- `OMM_munich-csm-aam_20241114`[1:120,1]
in_alt<- setdiff(csm_alt[[77]]@react_id, csm_aam[["result.77"]]@react_id)%>% str_replace_all("_c0","")

`OMM_munich-csm-m10_20241119` <- readRDS("OMM_munich-csm-m10_20241119.RDS")
csm_m10 <- `OMM_munich-csm-aam_20241119`[1:120,1]
in_m10<- setdiff(csm_m10[[77]]@react_id, csm_aam[["result.77"]]@react_id)%>% str_replace_all("_c0","")

#extract trehalose uptake fluxes

trehalose_mouse <- sim_OMM12_mouse@exchangeslist[[4]]%>%
  select(species, EX_cpd00794_e0)%>%
  separate_wider_delim(species, delim = "-", names = c("species", "model"), cols_remove = TRUE)

mouse_tre_summary <- trehalose_mouse %>%
  group_by(species) %>%
  summarise(
    n_cells_mouse = n(),
    total_flux_mouse = sum(abs(EX_cpd00794_e0), na.rm = TRUE),
    flux_per_cell_mouse = total_flux_mouse / n_cells_mouse,
    .groups = "drop") %>%
  mutate(
    flux_fraction_mouse = total_flux_mouse / sum(total_flux_mouse))
    


trehalose_apf <- sim_OMM12_m10@exchangeslist[[4]]%>%
  select(species, EX_cpd00794_e0)%>%
  separate_wider_delim(species, delim = "-", names = c("species", "model"), cols_remove = TRUE)

apf_tre_summary <- trehalose_apf %>%
  group_by(species) %>%
  summarise(
    n_cells_apf = n(),
    total_flux_apf = sum(abs(EX_cpd00794_e0), na.rm = TRUE),
    flux_per_cell_apf = total_flux_apf / n_cells_apf,
    .groups = "drop") %>%
  mutate(
    flux_fraction_apf = total_flux_apf / sum(total_flux_apf))



trehalose_af <- sim_OMM12_aam@exchangeslist[[4]]%>%
  select(species, EX_cpd00794_e0)%>%
  separate_wider_delim(species, delim = "-", names = c("species", "model"), cols_remove = TRUE)

af_tre_summary <- trehalose_af %>%
  group_by(species) %>%
  summarise(
    n_cells_af = n(),
    total_flux_af = sum(abs(EX_cpd00794_e0), na.rm = TRUE),
    flux_per_cell_af = total_flux_af / n_cells_af,
    .groups = "drop") %>%
  mutate(
    flux_fraction_af = total_flux_af / sum(total_flux_af))

all_summary<- full_join(mouse_tre_summary, apf_tre_summary)%>%full_join(af_tre_summary)%>%
  select(species, flux_fraction_af, flux_fraction_apf, flux_fraction_mouse)%>%
  pivot_longer(!species, names_to = "simulation", values_to = "flux_fraction")


all_summary$species <- factor(all_summary$species, 
                              levels=c("KB1", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58"))


col_annaw <- load("L:/Members/!Anna Burrichter/R/Colorscheme_AnnaW_ordered.rda")
col_reduced <- c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc", "#871900", "#365934",
                 "#ff9b05", "#004000", "#0b6760")

ggplot(all_summary, aes(x = simulation, y = flux_fraction, fill=species)) + 
  geom_col() + scale_fill_manual(values = col_reduced) +
  ylab("predicted flux fraction") + 
  theme_bw()+
  theme(legend.title=element_blank(), axis.title.x =element_blank(), legend.position = "left", axis.text.x = element_text(angle=90))


save(trehalose_mouse, file = "trehalose_uptake_mouse_model.rda")
save(trehalose_apf, file = "trehalose_uptake_apf_model.rda")
save(trehalose_af, file = "trehalose_uptake_af_model.rda")
