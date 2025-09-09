library(tidyverse)
library(readxl)
library(prism2R)
# read data from working directory
annotations <- read.delim("protein_annotation.tsv")%>% mutate(prodigal=NULL, id=NULL, ffn=NULL) %>% rename(protein_id_new = protein_id)
data <- read_excel("reanalysed_mouse_vs_invitro.xlsx", 
                                 sheet = 1)%>% rename(locus_tag = `locus tag_old`)

species <- c("E. faecalis KB1", "B. animalis YL2", "A. muris KB18", "M. intestinale YL27", "F. plautii YL31",
             "E. clostridioformis YL32", "A. muciniphila YL44", "T. muris YL45", "C. innocuum I46", 
             "B. caecimuris I48", "L. reuteri I49", "B. coccoides YL58")
species_short <- c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58")

##join annotations and proteome data
data_ann <- left_join(data, annotations, by = "locus_tag") %>%
            select(locus_tag, species, `refseq annotation`, koeggNOG, KEGG_ko, AAM_1:mouse_3)

#filter by species and normalise columns by column total
for (i in 1:12) {
  org <- species[i]
  suppressWarnings(
  pd <- data_ann %>% 
    filter(species == org) %>%
    mutate(across(AAM_1:mouse_3, as.numeric)) %>% 
    mutate(across(AAM_1:mouse_3, ~10^.x)) %>% 
    mutate(across(AAM_1:mouse_3, ~./sum(., na.rm = TRUE))))
  assign(paste(species_short[i], "normalised", sep = "_"), pd)
  rm(org, pd)
}
all_normalised <- bind_rows(KB1_normalised, YL2_normalised, KB18_normalised, YL27_normalised,
                           YL31_normalised, YL32_normalised, YL44_normalised, YL45_normalised,
                           I46_normalised, I48_normalised, I49_normalised, YL58_normalised) %>% 
                  select(locus_tag, AAM_1:mouse_3)

write.csv(all_normalised, "all_normalised.csv")

all_normalised_annotated <- bind_rows(KB1_normalised, YL2_normalised, KB18_normalised, YL27_normalised,
                            YL31_normalised, YL32_normalised, YL44_normalised, YL45_normalised,
                            I46_normalised, I48_normalised, I49_normalised, YL58_normalised)

write.csv(all_normalised_annotated, "all_normalised_annotated.csv")

## -> t-tests/ANOVA in PRISM