install.packages("devtools")
devtools::install_github("Waschina/gapseq.biochem")
library(gapseq.biochem)
library(tidyverse)
DB<-initDB()


NH3_reactions <- cpdRxns(DB, "cpd00013")
NH3_consuming <- NH3_reactions$consuming
NH3_con_prod <- NH3_reactions$`consuming+producing`

I46_mouse_3_mfluxlist_names <- names(sim_OMM12_mouse@mfluxlist[[3]][["I46-mouse_all"]])
I46_mouse_3_mfluxlist <- data.frame(sim_OMM12_mouse@mfluxlist[[3]][["I46-mouse_all"]]) %>% mutate(rxn = I46_mouse_3_mfluxlist_names)

I46_mouse_4_NH3 <- filter(I46_mouse_4_mfluxlist, rxn %in% NH3_reactions$consuming)
I46_mouse_4_NH3 <- filter(I46_mouse_4_mfluxlist, rxn %in% NH3_reactions$`consuming+producing`)

ornithine_reactions <- cpdRxns(DB, "cpd00064")
YL45_mouse_4_mfluxlist_names <- names(sim_OMM12_mouse@mfluxlist[[4]][["YL45-mouse_all"]])
YL45_mouse_4_mfluxlist <- data.frame(sim_OMM12_mouse@mfluxlist[[4]][["YL45-mouse_all"]]) %>% mutate(rxn = YL45_mouse_4_mfluxlist_names)
YL45_mouse_4_NH3 <- filter(YL45_mouse_4_mfluxlist, rxn %in% ornithine_reactions$`consuming+producing`)
