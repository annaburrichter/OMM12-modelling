library(clusterProfiler)
library(readxl)
library(sybil)
library(ggplot2)
library(gapseq.biochem)

db <- initDB("gapseq/") # needs gapseq directory (git clone https://github.com/jotech/gapseq)
models <- readRDS("dat/OMM_munich-csm-aam_20241205.RDS")
rm_org <- c("KB18", "I49", "YL2")
models.rel <- models[grep(paste0(rm_org, collapse="|"), names(models), invert = TRUE)]

biosynthesis <- c("Aminoacyl-tRNAs-Charging","AROMATIC-COMPOUNDS-BIOSYN","Carbohydrates-Biosynthesis","Cell-Structure-Biosynthesis","Cofactor-Biosynthesis","Lipid-Biosynthesis","Metabolic-Regulators","PYR-NUC-SYN","PUR-NUC-SYN","Other-biosynthesis","Polyprenyl-Biosynthesis","SECONDARY-METABOLITE-BIOSYNTHESIS","Storage-Compounds-Biosynthesis","Tetrapyrrole-Biosynthesis")
degradation <- c("Alcohol-Degradation","Aldehyde-Degradation","AMINE-DEG","Amino-Acid-Degradation","AROMATIC-COMPOUNDS-DEGRADATION","C1-COMPOUNDS","Carbohydrates-Degradation","CARBOXYLATES-DEG","CHLORINATED-COMPOUNDS-DEG","COFACTOR-DEGRADATION","Other-Degradation","Fatty-Acid-and-Lipid-Degradation","HORMONE-DEG","Noncarbon-Nutrients","NUCLEO-DEG","Polymer-Degradation", "Protein-Degradation","SECONDARY-METABOLITE-DEGRADATION")
energy <- c("Acetyl-CoA-Biosynthesis","CHEMOAUTOTROPHIC-ENERGY-METABOLISM","Electron-Transfer","Entner-Duodoroff-Pathways","Fermentation","GLYCOLYSIS-VARIANTS","Hydrogen-Production","OTHER-ENERGY", "Pentose-Phosphate-Cycle","Photosynthesis","Respiration","TCA-VARIANTS")
detox <- c("8-Oxo-GTP-Detoxification","Acid-Resistance","Antibiotic-Resistance","Arsenic-Detoxification","Cyanide-Detoxification","Mercury-Detoxification","Methylglyoxal-Detoxification","REACTIVE-OXYGEN-SPECIES-DEGRADATION","Seleno-Amino-Acid-Detoxification")
glycan <- c("Glycan-Biosynthesis","Glycan-Degradation")
macro <- c("Nucleic-Acid-Processing","Protein-Modification")
activation <- c("Activation","Inactivation","Interconversion")
#
sub_level1 <- c("Activation-Inactivation-Interconversion","Bioluminescence","Biosynthesis","Degradation","Detoxification","Energy-Metabolism","Glycan-Pathways","Macromolecule-Modification","Metabolic-Clusters","Signaling-Pathways","Transport-Pathways")
sub_level2 <- c(biosynthesis, degradation,energy,detox,glycan,macro,activation, "Transport")

# reaction list from models										  
mod.rxn.dt <- rbindlist(lapply(seq_along(models), function(i) { data.table(id=names(models)[i], rxn=models[[i]]@react_id) } ))
mod.rxn.rel.dt <- rbindlist(lapply(seq_along(models.rel), function(i) { data.table(id=names(models)[i], rxn=models[[i]]@react_id) } ))
mod.rxn.dt[!grepl("batch|mouse|AAM|M10|strict", id) ,source:="mono"]
mod.rxn.dt[grepl("AAM_all", id) ,source:="AF"]
mod.rxn.dt[grepl("M10_all", id) ,source:="APF"]
mod.rxn.dt[grepl("mouse_all", id) ,source:="vivo"]
mod.rxn.dt <- mod.rxn.dt[!is.na(source)]
mod.rxn.count.dt <- mod.rxn.dt[grepl("^rxn[0-9]+", rxn), list(N.mono=sum(source=="mono"), N.af=sum(source=="AF"), N.apf=sum(source=="APF"), N.vivo=sum(source=="vivo")),by=rxn]
mod.rxn.rel.dt[!grepl("batch|mouse|AAM|M10|strict", id) ,source:="mono"]
mod.rxn.rel.dt[grepl("AAM_all", id) ,source:="AF"]
mod.rxn.rel.dt[grepl("M10_all", id) ,source:="APF"]
mod.rxn.rel.dt[grepl("mouse_all", id) ,source:="vivo"]
mod.rxn.rel.dt <- mod.rxn.dt[!is.na(source)]
mod.rxn.rel.count.dt <- mod.rxn.rel.dt[grepl("^rxn[0-9]+", rxn), list(N.mono=sum(source=="mono"), N.af=sum(source=="AF"), N.apf=sum(source=="APF"), N.vivo=sum(source=="vivo")),by=rxn]

# subsystems
meta_pwy_dt <- db@pwy$main[, list(seed = unlist(str_split(id, " "))), by = .(meta.rxn.id, meta.rxn.name, meta.pwy, meta.sub)]
meta_pwy_dt[, meta.sub.lst := strsplit(meta.sub, ",")]
meta_pwy_dt <- meta_pwy_dt[, .(sub = unlist(meta.sub.lst)), by = .(seed, meta.rxn.id, meta.rxn.name, meta.pwy, meta.sub)]
meta_pwy_dt <- meta_pwy_dt[seed %in% str_extract(mod.rxn.dt$rxn, "rxn[0-9]+")] # reference subsystems are the ones available in all genomes
sub1_dt <- meta_pwy_dt[sub %in% sub_level1]
sub2_dt <- meta_pwy_dt[sub %in% sub_level2]

# vitro vs. mono
ora.mono.lst <- mod.rxn.rel.count.dt[N.af / N.mono > 1 , str_extract(rxn,"rxn[0-9]+")] # vitro:=AF (monocultures were also done with AF)
ora.mono <- enricher(gene=ora.mono.lst, TERM2GENE=unique(sub2_dt[, .(sub,seed)]))
clusterProfiler::dotplot(ora.mono, showCategory=100) + xlab("Reaction ratio") + ggtitle("vitro vs. mono (up)")
ggsave("img/ora.vitro_mono-up.pdf", width=7, height=2.4)
ggsave("img/ora.vitro_mono-up.png", width=7, height=2.4)
#
ora.mono.lst <- mod.rxn.rel.count.dt[N.af / N.mono < 1 , str_extract(rxn,"rxn[0-9]+")]
ora.mono <- enricher(gene=ora.mono.lst, TERM2GENE=unique(sub2_dt[, .(sub,seed)])) 
clusterProfiler::dotplot(ora.mono, showCategory=100) + xlab("Reaction ratio") + ggtitle("vitro vs. mono (down)")
ggsave("img/ora.vitro_mono-down.pdf", width=7, height=6)
ggsave("img/ora.vitro_mono-down.png", width=7, height=6)


# vivo vs. vitro
ora.vivo.lst <- mod.rxn.rel.count.dt[N.vivo / N.apf > 1, str_extract(rxn,"rxn[0-9]+")] # vitro:=APF 
ora.vivo <- enricher(gene=ora.vivo.lst, TERM2GENE=unique(sub2_dt[, .(sub,seed)]))
clusterProfiler::dotplot(ora.vivo, showCategory=100) + xlab("Reaction ratio") + ggtitle("vivo vs. vitro (up)")
ggsave("img/ora.vivo_vitro-up.pdf", width=7, height=3.2)
ggsave("img/ora.vivo_vitro-up.png", width=7, height=3.2)
#
ora.vivo.lst <- mod.rxn.rel.count.dt[N.vivo / N.apf < 1 , str_extract(rxn,"rxn[0-9]+")]
ora.vivo <- enricher(gene=ora.vivo.lst, TERM2GENE=unique(sub2_dt[, .(sub,seed)]))
clusterProfiler::dotplot(ora.vivo, showCategory=100) + xlab("Reaction ratio") + ggtitle("vivo vs. vitro (down)")
ggsave("img/ora.vivo_vitro-down.pdf", width=7, height=6)
ggsave("img/ora.vivo_vitro-down.png", width=7, height=6)

