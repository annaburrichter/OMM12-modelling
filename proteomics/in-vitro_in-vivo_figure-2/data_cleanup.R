library(tidyverse)
library(readxl)
library(prism2R)
# read data from working directory
annotations <- read.delim("protein_annotation.tsv")%>% mutate(prodigal=NULL, id=NULL, ffn=NULL) %>% rename(protein_id_new = protein_id)


data <- read_excel("reanalysed_mouse_vs_invitro.xlsx", sheet = 1,
                   col_types = c("skip", "text", "text", "text", "text", "text", "text", "text", 
                                 "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric"))%>% 
        rename(locus_tag = `locus tag_old`)


species_long <- c("E. faecalis KB1", "B. animalis YL2", "A. muris KB18", "M. intestinale YL27", "F. plautii YL31",
             "E. clostridioformis YL32", "A. muciniphila YL44", "T. muris YL45", "C. innocuum I46", 
             "B. caecimuris I48", "L. reuteri I49", "B. coccoides YL58")
species_short <- c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58")

##join annotations and proteome data
data_ann <- left_join(data, annotations, by = "locus_tag") %>%
            select(locus_tag, species, `refseq annotation`, koeggNOG, KEGG_ko, AAM_1:mouse_3)

species_normalised_list <- list()
all_normalised <- NULL

#filter by species and normalise columns by column total
for (i in 1:12) {
  pd <-   filter(data, species == species_long[i])
  
  #if there is only one value for a protein, set this to NA also  
  cols <- c(6,7,8)  # AAM
  pd <- pd %>% mutate(across(all_of(cols),
                                   ~ replace(.,rowSums(!is.na(pick(all_of(cols)))) < 2,NA)))
  cols <- c(9, 10, 11)  # M10
  pd <- pd %>% mutate(across(all_of(cols),
                                   ~ replace(.,rowSums(!is.na(pick(all_of(cols)))) < 2,NA)))
  cols <- c(12, 13, 14)  # mouse
  pd <- pd %>% mutate(across(all_of(cols),
                                   ~ replace(.,rowSums(!is.na(pick(all_of(cols)))) < 2,NA)))
  rm(cols)
  pd <- pd %>%
    mutate(across(AAM_1:mouse_3, ~10^.x)) %>% 
    mutate(across(AAM_1:mouse_3, ~./sum(., na.rm = TRUE)))
  
all_normalised <- bind_rows(all_normalised, pd)
species_normalised_list[[species_short[i]]] <- pd
pd <- select(pd, locus_tag, AAM_1:mouse_3)
filename <- paste(species_short[i], "_normalised.csv")
write.csv(pd, filename, row.names = FALSE)
rm(filename, pd)  
}

write.csv(all_normalised, "all_normalised.csv", row.names = FALSE)

## -> t-tests/ANOVA in PRISM