library(tidyverse)
library(readxl)
library(openxlsx)


# read and clean up data --------------------------------------------------

annotations <- read.delim2("protein_annotation.tsv")%>% mutate(prodigal=NULL, id=NULL, ffn=NULL) %>% rename(protein_id_new = protein_id)

#read single vs. community data, replace I48 dataset with separate one (contamination of original)

single_vs_community <- read_excel("reanalysed_single_vs_community.xlsx", 
                                  sheet = "single vs community",
                                  col_types = c("skip", "text", "text", "skip", "skip", "text", "text", "skip", 
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))%>%
  separate_wider_delim(2, delim = " ", names = c("genus", "species", "strain"))%>%
  filter(strain!= "I48")

i48_vs_community <- read_excel("reanalysed_single_vs_community.xlsx", 
                               sheet = "single vs community_I48", 
                               col_types = c("skip", "text", "text", "skip", "skip", "text", "text", "skip", 
                                             "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))%>%
  separate_wider_delim(2, delim = " ", names = c("genus", "species", "strain"))

single_vs_community<- bind_rows(single_vs_community, i48_vs_community)
rm(i48_vs_community)

single_vs_community <- single_vs_community %>%
  pivot_longer(
    cols = matches("^(KB1|KB18|YL2|YL27|YL31|YL32|YL44|YL45|YL58|I46|I48|I49)_[123]$"),
    names_to = c("monoculture_species", ".value"),
    names_pattern = "(.*)_([123])"
  ) %>%
  rename(
    monoculture_1 = `1`,
    monoculture_2 = `2`,
    monoculture_3 = `3`
  ) %>%
  filter(monoculture_species == strain)%>%
  rename(locus_tag = `locus tag_old`)%>%
  rename(community_1 = `community(batch)_1`)%>%
  rename(community_2 = `community(batch)_2`)%>%
  rename(community_3 = `community(batch)_3`)%>%
  left_join(annotations)

#read in vitro vs. in vivo data

vitro_vs_vivo <- read_excel("reanalysed_mouse_vs_invitro.xlsx", sheet = 1,
                   col_types = c("skip", "text", "text", "text", "text", "text", "text", "text", 
                                 "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric"))%>% 
  separate_wider_delim(2, delim = " ", names = c("genus", "species", "strain"))%>%
  rename(locus_tag = `locus tag_old`)%>%
  left_join(annotations)

#cleanup: if there is only one value for a protein, set this to NA also, de-log data

single_vs_community <- single_vs_community %>% 
  mutate(across(all_of(c(7,8,9)),~ replace(.,rowSums(!is.na(pick(all_of(c(7,8,9))))) < 2,NA)))%>% #community data columns
  mutate(across(all_of(c(11,12,13)),~ replace(.,rowSums(!is.na(pick(all_of(c(11,12,13))))) < 2,NA)))%>% #monoculture data columns
  mutate(across(c(7:9, 11:13), ~ 10^.))
  
vitro_vs_vivo <- vitro_vs_vivo %>% 
  mutate(across(all_of(c(10,11,12)),~ replace(.,rowSums(!is.na(pick(all_of(c(10,11,12))))) < 2,NA)))%>% #AF/AAM data columns
  mutate(across(all_of(c(13,14,15)),~ replace(.,rowSums(!is.na(pick(all_of(c(13,14,15))))) < 2,NA)))%>% #APF/M10 data columns
  mutate(across(all_of(c(16,17,18)),~ replace(.,rowSums(!is.na(pick(all_of(c(16,17,18))))) < 2,NA)))%>% #mouse data columns
  mutate(across(c(10:18), ~ 10^.))


# calculate mass fraction per COG category --------------------------------


# single vs. community ----------------------------------------------------


#separate rows with multiple COG categories (those will be duplicated so they count twice)
single_vs_community_COG <- single_vs_community %>%
  mutate(COG_category = if_else(is.na(COG_category) | COG_category == "","-",COG_category))%>%
  mutate(COG_category = str_split(COG_category, "")) %>%unnest(COG_category)

value_cols <- c("community_1", "community_2", "community_3","monoculture_1", "monoculture_2", "monoculture_3")

#calculate fractions

COG_fraction_species_single_comm <- single_vs_community_COG %>%
  group_by(strain, COG_category) %>%
  summarise(across(all_of(value_cols), ~ sum(.x, na.rm = TRUE)),.groups = "drop") %>%
  group_by(strain) %>%
  mutate(across(all_of(value_cols),~ .x / sum(.x, na.rm = TRUE),.names = "{.col}_fraction")) %>%
  ungroup()%>%
  select(strain,COG_category,
    monoculture_1_fraction,monoculture_2_fraction,monoculture_3_fraction,
    community_1_fraction,community_2_fraction,community_3_fraction)

COG_fraction_all_single_comm <- single_vs_community_COG %>%
  group_by(COG_category) %>%
  summarise(across(all_of(value_cols), ~ sum(.x, na.rm = TRUE)),.groups = "drop") %>%
  mutate(across(all_of(value_cols),~ .x / sum(.x, na.rm = TRUE),.names = "{.col}_fraction")) %>%
  select(COG_category,
         monoculture_1_fraction,monoculture_2_fraction, monoculture_3_fraction, 
         community_1_fraction,community_2_fraction,community_3_fraction)

species_list <- split(COG_fraction_species_single_comm, COG_fraction_species_single_comm$strain)
species_list <- lapply(species_list,function(x) select(x, -strain))
species_list <- c(species_list,
  list(all_species = COG_fraction_all_single_comm))
species_list <- map(species_list, function(df) {
    missing_categories <- tibble(COG_category = c("A", "Z", "W")) %>%
    anti_join(df, by = "COG_category")
  bind_rows(df,missing_categories %>%
      mutate(across(where(is.numeric),~ 0))) %>%
    arrange(COG_category)})

write.xlsx(species_list,
  file = "COG_mass_fractions_single_comm.xlsx", overwrite = TRUE)

#group categories

species_list_grouped <- map(species_list, function(df) {
  
  df %>%mutate(
      group = case_when(
        COG_category == "C" ~ "catabolism",
        COG_category %in% c("D", "J", "K", "L") ~ "growth",
        COG_category %in% c("G", "N", "T") ~ "resource acquisition",
        TRUE ~ NA_character_)) %>%
    filter(!is.na(group)) %>%
    group_by(group) %>%
    summarise(
      across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop") %>%
    rename(COG_category = group)})

write.xlsx(species_list_grouped,file = "COG_mass_fractions_single_comm_grouped.xlsx", overwrite = TRUE)


# vitro vs. vivo ----------------------------------------------------


#separate rows with multiple COG categories (those will be duplicated so they count twice)
vitro_vs_vivo_COG <- vitro_vs_vivo %>%
  mutate(COG_category = if_else(is.na(COG_category) | COG_category == "","-",COG_category))%>%
  mutate(COG_category = str_split(COG_category, "")) %>%unnest(COG_category)

value_cols <- c("AAM_1", "AAM_2", "AAM_3","M10_1", "M10_2", "M10_3", "mouse_1", "mouse_2", "mouse_3")

#calculate fractions

COG_fraction_species_vitro_vivo <- vitro_vs_vivo_COG %>%
  group_by(strain, COG_category) %>%
  summarise(across(all_of(value_cols), ~ sum(.x, na.rm = TRUE)),.groups = "drop") %>%
  group_by(strain) %>%
  mutate(across(all_of(value_cols),~ .x / sum(.x, na.rm = TRUE),.names = "{.col}_fraction")) %>%
  ungroup()%>%
  select(strain,COG_category,
         AAM_1_fraction,AAM_2_fraction,AAM_3_fraction,M10_1_fraction,M10_2_fraction,M10_3_fraction,
         mouse_1_fraction,mouse_2_fraction,mouse_3_fraction)

COG_fraction_all_vitro_vivo <- vitro_vs_vivo_COG %>%
  group_by(COG_category) %>%
  summarise(across(all_of(value_cols), ~ sum(.x, na.rm = TRUE)),.groups = "drop") %>%
  mutate(across(all_of(value_cols),~ .x / sum(.x, na.rm = TRUE),.names = "{.col}_fraction")) %>%
  select(COG_category,
         AAM_1_fraction,AAM_2_fraction,AAM_3_fraction,M10_1_fraction,M10_2_fraction,M10_3_fraction,
         mouse_1_fraction,mouse_2_fraction,mouse_3_fraction)

species_list <- split(COG_fraction_species_vitro_vivo, COG_fraction_species_vitro_vivo$strain)
species_list <- lapply(species_list,function(x) select(x, -strain))
species_list <- c(species_list,
                  list(all_species = COG_fraction_all_vitro_vivo))
species_list <- map(species_list, function(df) {
  missing_categories <- tibble(COG_category = c("A", "Z", "W")) %>%
    anti_join(df, by = "COG_category")
  bind_rows(df,missing_categories %>%
              mutate(across(where(is.numeric),~ 0))) %>%
    arrange(COG_category)})
write.xlsx(species_list,
           file = "COG_mass_fractions_vitro_vivo.xlsx", overwrite = TRUE)

#group categories

species_list_grouped <- map(species_list, function(df) {
  
  df %>%mutate(
    group = case_when(
      COG_category == "C" ~ "catabolism",
      COG_category %in% c("D", "J", "K", "L") ~ "growth",
      COG_category %in% c("G", "N", "T") ~ "resource acquisition",
      TRUE ~ NA_character_)) %>%
    filter(!is.na(group)) %>%
    group_by(group) %>%
    summarise(
      across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop") %>%
    rename(COG_category = group)})

write.xlsx(species_list_grouped,file = "COG_mass_fractions_vitro_vivo_grouped.xlsx", overwrite = TRUE)
