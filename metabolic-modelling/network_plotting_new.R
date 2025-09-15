library(BacArena)
library(foreach)
library(doParallel)
library(stringr)
library(data.table)
library(tidyverse)
library(igraph)
library(readxl)
library(inlcolor)

gapseq_compounds <- read_excel("gapseq compound IDs.xlsx")


# interactions
network_mouse <- findFeeding3(sim_OMM12_mouse, mets = sim_OMM12_mouse@mediac, time = 3)
network_af <- findFeeding3(sim_OMM12_aam,  mets = sim_OMM12_aam@mediac, time = 3)
network_apf <- findFeeding3(sim_OMM12_m10,  mets = sim_OMM12_aam@mediac, time = 3)




# colour codes species
V(network_mouse[[2]])
colors_org_mouse <- c("#6b00bc","#871900", "#ff9b05", "#0b6760", "#004000", "#11a196", "#0e8b18", "#cc6a00", "#365934", "#7bc57b")
V(network_af[[2]])
colors_org_af <- c("#871900", "#cc6a00", "#7bc57b","#0e8b18","#6b00bc","#ff9b05", "#0b6760", "#11a196", "#365934")
V(network_apf[[2]])
colors_org_apf <- c("#871900", "#0b6760", "#7bc57b", "#ff9b05", "#0e8b18","#6b00bc",  "#365934", "#11a196")

# colour codes and legends metabolites
met_af <- as.data.frame(matrix(levels(network_af[[1]][["met"]]), ncol=1))
met_apf <- as.data.frame(matrix(levels(network_apf[[1]][["met"]]), ncol=1))
met_mouse <- as.data.frame(matrix(levels(network_mouse[[1]][["met"]]), ncol=1))
met_all <- unique(rbind(met_af, met_apf, met_mouse))
met_all$V1 <- str_replace_all(met_all$V1, "EX_", "")
met_all$V1 <- str_replace_all(met_all$V1, "_e0", "")
names(met_all) <- c("id")
met_all <- left_join(met_all, gapseq_compounds)
colors <- as.character(get_colors(nrow(met_all), scheme = "plasma"))
met_all$color <- colors

#get consistent coordinates
testplot_mouse <- tkplot(network_mouse[[2]])
layout_mouse <- tk_coords(testplot_mouse)
met_mouse$V1 <- str_replace_all(met_mouse$V1, "EX_", "")
met_mouse$V1 <- str_replace_all(met_mouse$V1, "_e0", "")
names(met_mouse) <- c("id")
legend_mouse <- left_join(met_mouse, met_all)
legend_mouse_text <- legend_mouse$name
colors_met_mouse <- as.character(legend_mouse$color)
plot(network_mouse[[2]], layout=layout_mouse, edge.color = colors_met_mouse[network_mouse[[1]][["met"]]], edge.arrow.size =0.75, edge.width =2, vertex.color =colors_org_mouse, vertex.size = 35, vertex.label.color ="white", vertex.label = c("YL44", "YL45", "I48", "YL58", "I49", "YL32", "KB1", "YL27", "I46", "YL31")) +
 legend(x=-2, y=1, legend=legend_mouse_text, fill=colors_met_mouse, cex =0.8)



testplot_apf <- tkplot(network_apf[[2]])
layout_apf <- tk_coords(testplot_apf)
met_apf$V1 <- str_replace_all(met_apf$V1, "EX_", "")
met_apf$V1 <- str_replace_all(met_apf$V1, "_e0", "")
names(met_apf) <- c("id")
legend_apf <- left_join(met_apf, met_all)
legend_apf_text <- legend_apf$name
colors_met_apf <- as.character(legend_apf$color)
plot(network_apf[[2]], layout=layout_apf, edge.color = colors_met_apf[network_apf[[1]][["met"]]], edge.arrow.size =0.75, edge.width =2, vertex.color =colors_org_apf, vertex.size = 35, vertex.label.color ="white", vertex.label = c("YL45", "YL58", "YL31", "I48", "KB1", "YL44", "I46", "YL32")) +
  legend(x=-2, y=1, legend=legend_apf_text, fill=colors_met_apf, cex =0.8)



testplot_af <- tkplot(network_af[[2]])
layout_af <- tk_coords(testplot_af)
met_af$V1 <- str_replace_all(met_af$V1, "EX_", "")
met_af$V1 <- str_replace_all(met_af$V1, "_e0", "")
names(met_af) <- c("id")
legend_af <- left_join(met_af, met_all)
legend_af_text <- legend_af$name
colors_met_af <- as.character(legend_af$color)
plot(network_af[[2]], layout=layout_mouse, edge.color = colors_met_af[network_af[[1]][["met"]]], edge.arrow.size =0.75, edge.width =2, vertex.color =colors_org_af, vertex.size = 35, vertex.label.color ="white", vertex.label = c("YL45", "YL27", "YL31", "KB1", "YL44", "I48", "YL58", "YL32", "I46")) +
  legend(x=-2, y=1, legend=legend_af_text, fill=colors_met_af, cex =0.8)



## graph theory outputs

degree_af <- degree(network_af[[2]])
degree_in_af <- degree(network_af[[2]], mode = "in")
degree_out_af <- degree(network_af[[2]], mode = "out")

degree_apf <- degree(network_apf[[2]])
degree_in_apf <- degree(network_apf[[2]], mode = "in")
degree_out_apf <- degree(network_apf[[2]], mode = "out")

degree_mouse <- degree(network_mouse[[2]])
degree_in_mouse <- degree(network_mouse[[2]], mode = "in")
degree_out_mouse <- degree(network_mouse[[2]], mode = "out")

#print everything (no better solution yet)
degree_af
degree_in_af
degree_out_af
degree_apf
degree_in_apf
degree_out_apf
degree_mouse
degree_in_mouse
degree_out_mouse


write.csv(network_af[[1]], file="network_fluxes_af.csv")
write.csv(network_apf[[1]], file="network_fluxes_apf.csv")
write.csv(network_mouse[[1]], file="network_fluxes_mouse.csv")

#network fluxes with compound names
network_fluxes_af <- network_af[[1]] %>% 
  separate_wider_delim(met, delim = "_", names = c("EX", "compound", "c0"), cols_remove = FALSE)%>%
  mutate(EX = NULL) %>% mutate(c0=NULL)%>%
  left_join(gapseq_compounds, by = c('compound' = 'id'))%>%
  mutate(MNX_ID = NULL)%>% select(name, formula, charge, compound, sim_step, prod, cons, prod.flux, cons.flux)
write.csv(network_fluxes_af, file="network_fluxes_af_named.csv")

network_fluxes_apf <- network_apf[[1]] %>% 
  separate_wider_delim(met, delim = "_", names = c("EX", "compound", "c0"), cols_remove = FALSE)%>%
  mutate(EX = NULL) %>% mutate(c0=NULL)%>%
  left_join(gapseq_compounds, by = c('compound' = 'id'))%>%
  mutate(MNX_ID = NULL)%>% select(name, formula, charge, compound, sim_step, prod, cons, prod.flux, cons.flux)
write.csv(network_fluxes_apf, file="network_fluxes_apf_named.csv")

network_fluxes_mouse <- network_mouse[[1]] %>% 
  separate_wider_delim(met, delim = "_", names = c("EX", "compound", "c0"), cols_remove = FALSE)%>%
  mutate(EX = NULL) %>% mutate(c0=NULL)%>%
  left_join(gapseq_compounds, by = c('compound' = 'id'))%>%
  mutate(MNX_ID = NULL)%>% select(name, formula, charge, compound, sim_step, prod, cons, prod.flux, cons.flux)
write.csv(network_fluxes_mouse, file="network_fluxes_mouse_named.csv")

