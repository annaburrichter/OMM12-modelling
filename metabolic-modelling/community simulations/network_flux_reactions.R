install.packages("devtools")
devtools::install_github("Waschina/gapseq.biochem")
library(gapseq.biochem)
library(tidyverse)
DB<-initDB()

setdiff(met_af, met_apf)
setdiff(met_apf, met_af)
setdiff(met_af, met_mouse)
setdiff(met_mouse, met_af)
setdiff(met_apf, met_mouse)
met_vitro <- rbind(met_af, met_apf)
exchanges_mouse <- setdiff(met_mouse, met_apf)
vivo_only <- setdiff(met_mouse, met_vitro) %>% left_join(gapseq_compounds) %>% select(id, name)
vitro_only <- setdiff(met_vitro, met_mouse) %>% left_join(gapseq_compounds) %>% select(id, name)
write.csv(vivo_only, "exchanges_vivo_only.csv")
write.csv(vitro_only, "exchanges_vitro_only.csv")
exchanges_mouse <- left_join(exchanges_mouse, gapseq_compounds) %>% select(id, name)

for (i in (1:18)) {

reactions <- cpdRxns(DB, exchanges_mouse[i,1])

for (j in (1:10)) {
  
  species <- names(sim_OMM12_mouse@mfluxlist[[3]])[j]
  
  fluxes_names <- names(sim_OMM12_mouse@mfluxlist[[3]][[j]])
  fluxes_names <- str_replace_all(fluxes_names, "_c0", "")
  fluxes_mouse <- data.frame(sim_OMM12_mouse@mfluxlist[[3]][[j]]) %>% mutate(rxn = fluxes_names)
  flux_mouse_con <- filter(fluxes_mouse, rxn %in% reactions$consuming)
  flux_mouse_prod <- filter(fluxes_mouse, rxn %in% reactions$producing)
  flux_mouse_conprod <- filter(fluxes_mouse, rxn %in% reactions$`consuming+producing`)
  flux_mouse_trans <- filter(fluxes_mouse, rxn %in% reactions$transporters)
 
  if (nrow(flux_mouse_con) >0){
    write.csv(flux_mouse_con, paste(species, exchanges_mouse[i,2], "flux_mouse_con", "csv", sep = "."))}
  if (nrow(flux_mouse_prod) >0){
  write.csv(flux_mouse_prod, paste(species, exchanges_mouse[i,2], "flux_mouse_prod", "csv", sep = "."))}
  if (nrow(flux_mouse_conprod) >0){
  write.csv(flux_mouse_conprod, paste(species, exchanges_mouse[i,2], "flux_mouse_conprod", "csv", sep = "."))}
  if (nrow(flux_mouse_trans) >0){
  write.csv(flux_mouse_trans, paste(species, exchanges_mouse[i,2], "flux_mouse_trans", "csv", sep = "."))}
  }

for (h in (1:8)) {
  
  species <- names(sim_OMM12_m10@mfluxlist[[3]])[h]
  
  fluxes_names <- names(sim_OMM12_m10@mfluxlist[[3]][[h]])
  fluxes_names <- str_replace_all(fluxes_names, "_c0", "")
  fluxes_m10 <- data.frame(sim_OMM12_m10@mfluxlist[[3]][[h]]) %>% mutate(rxn = fluxes_names)
  flux_m10_con <- filter(fluxes_m10, rxn %in% reactions$consuming)
  flux_m10_prod <- filter(fluxes_m10, rxn %in% reactions$producing)
  flux_m10_conprod <- filter(fluxes_m10, rxn %in% reactions$`consuming+producing`)
  flux_m10_trans <- filter(fluxes_m10, rxn %in% reactions$transporters)
  
  if (nrow(flux_m10_con) >0){
    write.csv(flux_m10_con, paste(species, exchanges_mouse[i,2], "flux_m10_con", "csv", sep = "."))}
    if (nrow(flux_m10_prod) >0){
      write.csv(flux_m10_prod, paste(species, exchanges_mouse[i,2], "flux_m10_prod", "csv", sep = "."))}
      if (nrow(flux_m10_conprod) >0){
        write.csv(flux_m10_conprod, paste(species, exchanges_mouse[i,2], "flux_m10_conprod", "csv", sep = "."))}
        if (nrow(flux_m10_trans) >0){
          write.csv(flux_m10_trans, paste(species, exchanges_mouse[i,2], "flux_m10_trans", "csv", sep = "."))}
}

for (k in (1:9)) {
  
  species <- names(sim_OMM12_aam@mfluxlist[[3]])[k]
  
  fluxes_names <- names(sim_OMM12_aam@mfluxlist[[3]][[k]])
  fluxes_names <- str_replace_all(fluxes_names, "_c0", "")
  fluxes_aam <- data.frame(sim_OMM12_aam@mfluxlist[[3]][[k]]) %>% mutate(rxn = fluxes_names)
  flux_aam_con <- filter(fluxes_aam, rxn %in% reactions$consuming)
  flux_aam_prod <- filter(fluxes_aam, rxn %in% reactions$producing)
  flux_aam_conprod <- filter(fluxes_aam, rxn %in% reactions$`consuming+producing`)
  flux_aam_trans <- filter(fluxes_aam, rxn %in% reactions$transporters)
  
  if (nrow(flux_aam_con) >0){
    write.csv(flux_aam_con, paste(species, exchanges_mouse[i,2], "flux_aam_con", "csv", sep = "."))}
  if (nrow(flux_aam_prod) >0){
    write.csv(flux_aam_prod, paste(species, exchanges_mouse[i,2], "flux_aam_prod", "csv", sep = "."))}
  if (nrow(flux_aam_conprod) >0){
    write.csv(flux_aam_conprod, paste(species, exchanges_mouse[i,2], "flux_aam_conprod", "csv", sep = "."))}
  if (nrow(flux_aam_trans) >0){
    write.csv(flux_aam_trans, paste(species, exchanges_mouse[i,2], "flux_aam_trans", "csv", sep = "."))}
}

assign(paste("reactions", exchanges_mouse[i,2], sep = "_"), reactions)

}


