install.packages("remotes")
remotes::install_github("SysBioChalmers/sybil")
install.packages(c("ReacTran", "deSolve", "glpkAPI", "R.matlab"))
install.packages("//nas.ads.mwn.de/lmmepki$/AG_Stecher/!Anna Burrichter/metabolic models AG Kaleta/BacArena-1.8.1.tar.gz", repos = NULL, type = "source")

library(tidyverse)
library(sybil)
library(BacArena)
library(igraph)
library(data.table)

#read medium
medium_aam <- fread("//nas.ads.mwn.de/lmmepki$/AG_Stecher/!Anna Burrichter/metabolic models AG Kaleta/AAM_nomuc_rev5.csv")
medium_m10 <- fread("//nas.ads.mwn.de/lmmepki$/AG_Stecher/!Anna Burrichter/metabolic models AG Kaleta/m10_vj_tre.csv")


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

`OMM_munich-csm-aam_20241114` <- readRDS("//nas.ads.mwn.de/lmmepki$/AG_Stecher/!Anna Burrichter/metabolic models AG Kaleta/community modelling 11-2024/OMM_munich-csm-aam_20241114.RDS")
csm_alt <- `OMM_munich-csm-aam_20241114`[1:120,1]
in_alt<- setdiff(csm_alt[[77]]@react_id, csm_aam[["result.77"]]@react_id)%>% str_replace_all("_c0","")

`OMM_munich-csm-m10_20241119` <- readRDS("//nas.ads.mwn.de/lmmepki$/AG_Stecher/!Anna Burrichter/metabolic models AG Kaleta/community modelling 11-2024/OMM_munich-csm-m10_20241119.RDS")
csm_m10 <- `OMM_munich-csm-aam_20241119`[1:120,1]
in_m10<- setdiff(csm_m10[[77]]@react_id, csm_aam[["result.77"]]@react_id)%>% str_replace_all("_c0","")

Yl31_maus_m10 <- Bac(csm_m10[[77]])
arena_mouse_test <- Arena(n=50,m=50)
arena_mouse_test <- addOrg(arena_mouse_test, Yl31_maus_m10, amount=3)
arena_mouse_test <- addSubs(arena_mouse_test, mediac=paste0("EX_",medium_m10$compounds,"_e0"), smax=medium_m10$maxFlux, unit="mM")
sim_test <- simEnv(arena_mouse_test, time=3, sec_obj='mtf')

library(gapseq.biochem)
library(tidyverse)
DB<-initDB()
# show reactions -> TCA cycle reactions?
for (i in 11:20) {
  summary(DB, in_alt[i])
}
