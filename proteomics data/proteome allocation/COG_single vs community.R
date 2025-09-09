library(tidyverse)
library(readxl)
setwd("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/data by COG categories/single vs community")

#load data
single_vs_community <- read_xlsx("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/reanalysed_single_vs_community.xlsx", sheet = "single vs community")
annotations <-read.csv (file = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/protein_annotation.tsv", header = TRUE, sep = "\t")

#add further annotations
single_vs_community <- rename(single_vs_community, "locus_tag" = "locus tag_old" )
single_vs_community <- single_vs_community %>% left_join(annotations, by = "locus_tag") 
single_vs_community_unique <- distinct(single_vs_community, locus_tag, .keep_all = TRUE) 
single_vs_community_short <- select(single_vs_community_unique, "locus_tag", "species", "refseq annotation", "koeggNOG", "prodigal", "cds", "rxn", "ec", "pathway", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "COG_category", "community(batch)_1":"YL58_3") %>% mutate_at (c('KB1_1', 'KB1_2', 'KB1_3','KB18_1', 'KB18_2', 'KB18_3','YL2_1', 'YL2_2', 'YL2_3','YL27_1', 'YL27_2', 'YL27_3','YL31_1', 'YL31_2', 'YL31_3','YL32_1', 'YL32_2', 'YL32_3','YL44_1', 'YL44_2', 'YL44_3','YL45_1', 'YL45_2', 'YL45_3','I46_1', 'I46_2', 'I46_3','I49_1', 'I49_2', 'I49_3','YL58_1', 'YL58_2', 'YL58_3','community(batch)_1', 'community(batch)_2', 'community(batch)_3'), as.numeric)

 
#extract all community proteins and each single strain for graph via COG categories (keep only proteins in 2 out of 3 replicates)

community_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `community(batch)_1`, `community(batch)_2`, `community(batch)_3`) %>% 
  filter((!is.na(`community(batch)_1`) & !is.na(`community(batch)_2`)) | (!is.na(`community(batch)_1`) & !is.na(`community(batch)_3`)) | (!is.na(`community(batch)_2`) & !is.na(`community(batch)_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_community <- count(community_all, COG_category)
write.csv(categories_community, file = "COG_community(batch).csv")

KB1_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `KB1_1`, `KB1_2`, `KB1_3`) %>% 
  filter((!is.na(`KB1_1`) & !is.na(`KB1_2`)) | (!is.na(`KB1_1`) & !is.na(`KB1_3`)) | (!is.na(`KB1_2`) & !is.na(`KB1_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_KB1 <- count(KB1_all, COG_category)
write.csv(categories_KB1, file = "COG_KB1.csv")

YL2_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `YL2_1`, `YL2_2`, `YL2_3`) %>% 
  filter((!is.na(`YL2_1`) & !is.na(`YL2_2`)) | (!is.na(`YL2_1`) & !is.na(`YL2_3`)) | (!is.na(`YL2_2`) & !is.na(`YL2_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_YL2 <- count(YL2_all, COG_category)
write.csv(categories_YL2, file = "COG_YL2.csv")

KB18_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `KB18_1`, `KB18_2`, `KB18_3`) %>% 
  filter((!is.na(`KB18_1`) & !is.na(`KB18_2`)) | (!is.na(`KB18_1`) & !is.na(`KB18_3`)) | (!is.na(`KB18_2`) & !is.na(`KB18_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_KB18 <- count(KB18_all, COG_category)
write.csv(categories_KB18, file = "COG_KB18.csv")

YL27_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `YL27_1`, `YL27_2`, `YL27_3`) %>% 
  filter((!is.na(`YL27_1`) & !is.na(`YL27_2`)) | (!is.na(`YL27_1`) & !is.na(`YL27_3`)) | (!is.na(`YL27_2`) & !is.na(`YL27_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_YL27 <- count(YL27_all, COG_category)
write.csv(categories_YL27, file = "COG_YL27.csv")

YL31_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `YL31_1`, `YL31_2`, `YL31_3`) %>% 
  filter((!is.na(`YL31_1`) & !is.na(`YL31_2`)) | (!is.na(`YL31_1`) & !is.na(`YL31_3`)) | (!is.na(`YL31_2`) & !is.na(`YL31_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_YL31 <- count(YL31_all, COG_category)
write.csv(categories_YL31, file = "COG_YL31.csv")

YL32_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `YL32_1`, `YL32_2`, `YL32_3`) %>% 
  filter((!is.na(`YL32_1`) & !is.na(`YL32_2`)) | (!is.na(`YL32_1`) & !is.na(`YL32_3`)) | (!is.na(`YL32_2`) & !is.na(`YL32_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_YL32 <- count(YL32_all, COG_category)
write.csv(categories_YL32, file = "COG_YL32.csv")

YL44_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `YL44_1`, `YL44_2`, `YL44_3`) %>% 
  filter((!is.na(`YL44_1`) & !is.na(`YL44_2`)) | (!is.na(`YL44_1`) & !is.na(`YL44_3`)) | (!is.na(`YL44_2`) & !is.na(`YL44_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_YL44 <- count(YL44_all, COG_category)
write.csv(categories_YL44, file = "COG_YL44.csv")

YL45_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `YL45_1`, `YL45_2`, `YL45_3`) %>% 
  filter((!is.na(`YL45_1`) & !is.na(`YL45_2`)) | (!is.na(`YL45_1`) & !is.na(`YL45_3`)) | (!is.na(`YL45_2`) & !is.na(`YL45_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_YL45 <- count(YL45_all, COG_category)
write.csv(categories_YL45, file = "COG_YL45.csv")

I46_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `I46_1`, `I46_2`, `I46_3`) %>% 
  filter((!is.na(`I46_1`) & !is.na(`I46_2`)) | (!is.na(`I46_1`) & !is.na(`I46_3`)) | (!is.na(`I46_2`) & !is.na(`I46_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_I46 <- count(I46_all, COG_category)
write.csv(categories_I46, file = "COG_I46.csv")

I49_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `I49_1`, `I49_2`, `I49_3`) %>% 
  filter((!is.na(`I49_1`) & !is.na(`I49_2`)) | (!is.na(`I49_1`) & !is.na(`I49_3`)) | (!is.na(`I49_2`) & !is.na(`I49_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_I49 <- count(I49_all, COG_category)
write.csv(categories_I49, file = "COG_I49.csv")

YL58_all <- single_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `YL58_1`, `YL58_2`, `YL58_3`) %>% 
  filter((!is.na(`YL58_1`) & !is.na(`YL58_2`)) | (!is.na(`YL58_1`) & !is.na(`YL58_3`)) | (!is.na(`YL58_2`) & !is.na(`YL58_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_YL58 <- count(YL58_all, COG_category)
write.csv(categories_YL58, file = "COG_YL58.csv")

##for I48

#load data
I48_vs_community <- read_xlsx("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/reanalysed_single_vs_community.xlsx", sheet = "single vs community_I48")

#add further annotations
I48_vs_community <- rename(I48_vs_community, "locus_tag" = "locus tag_old" )
I48_vs_community <- I48_vs_community %>% left_join(annotations, by = "locus_tag") 
I48_vs_community_unique <- distinct(I48_vs_community, locus_tag, .keep_all = TRUE) 
I48_vs_community_short <- select(I48_vs_community_unique, "locus_tag", "species", "refseq annotation", "koeggNOG", "prodigal", "cds", "rxn", "ec", "pathway", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "COG_category", "community(batch)_1", "community(batch)_2", "community(batch)_3", "I48_1", "I48_2", "I48_3")

I48_all <- I48_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `I48_1`, `I48_2`, `I48_3`) %>% 
  filter((!is.na(`I48_1`) & !is.na(`I48_2`)) | (!is.na(`I48_1`) & !is.na(`I48_3`)) | (!is.na(`I48_2`) & !is.na(`I48_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)
categories_I48 <- count(I48_all, COG_category)
write.csv(categories_I48, file = "COG_I48.csv")

I48_community <- I48_vs_community_short %>% 
  select(locus_tag, species, koeggNOG, COG_category, `community(batch)_1`, `community(batch)_2`, `community(batch)_3`) %>% 
  filter((!is.na(`community(batch)_1`) & !is.na(`community(batch)_2`)) | (!is.na(`community(batch)_1`) & !is.na(`community(batch)_3`)) | (!is.na(`community(batch)_2`) & !is.na(`community(batch)_3`))) %>%
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)

#by replicate, community data
species_11 <- c("E. faecalis KB1", "B. animalis YL2", "A. muris KB18", "M. intestinale YL27", "F. plautii YL31", 
                "E. clostridioformis YL32", "A. muciniphila YL44", "T. muris YL45", "C. innocuum I46", 
                "L. reuteri I49", "B. coccoides YL58")

subset_1 <- filter(I48_vs_community_short, species == "B. caecimuris I48")%>%select(COG_category, `community(batch)_1`)%>%filter(!is.na(`community(batch)_1`))
count_1 <- count(subset_1, COG_category)
colnames(count_1) <- c("COG_category", paste("B. caecimuris I48", "_1"))
assign (paste("comm", "B. caecimuris I48", "1", sep = "_"), count_1)
comm_all <- count_1
subset_2 <- filter(I48_vs_community_short, species == "B. caecimuris I48")%>%select(COG_category, `community(batch)_2`)%>%filter(!is.na(`community(batch)_2`))
count_2 <- count(subset_2, COG_category)
colnames(count_2) <- c("COG_category", paste("B. caecimuris I48", "_2"))
assign (paste("comm", "B. caecimuris I48", "1", sep = "_"), count_2)
comm_all <- left_join(comm_all, count_2)
subset_3 <- filter(I48_vs_community_short, species == "B. caecimuris I48")%>%select(COG_category, `community(batch)_3`)%>%filter(!is.na(`community(batch)_3`))
count_3 <- count(subset_3, COG_category)
colnames(count_3) <- c("COG_category", "B. caecimuris I48", "_3")
assign (paste("comm", "B. caecimuris I48", "1", sep = "_"), count_3)
comm_all <- left_join(comm_all, count_3)
 for (i in 1:11){
   subset_1 <- filter(single_vs_community_short, species == species_11[i])%>%select(COG_category, `community(batch)_1`)%>%filter(!is.na(`community(batch)_1`))
   count_1 <- count(subset_1, COG_category)
   colnames(count_1) <- c("COG_category", paste(species_11[i], "_1"))
   assign (paste("comm", species_11[i], "1", sep = "_"), count_1)
   comm_all <- left_join(comm_all, count_1)
   subset_2 <- filter(single_vs_community_short, species == species_11[i])%>%select(COG_category, `community(batch)_2`)%>%filter(!is.na(`community(batch)_2`))
   count_2 <- count(subset_2, COG_category)
   colnames(count_2) <- c("COG_category", paste(species_11[i], "_2"))
   assign (paste("comm", species_11[i], "1", sep = "_"), count_2)
   comm_all <- left_join(comm_all, count_2)
   subset_3 <- filter(single_vs_community_short, species == species_11[i])%>%select(COG_category, `community(batch)_3`)%>%filter(!is.na(`community(batch)_3`))
   count_3 <- count(subset_3, COG_category)
   colnames(count_3) <- c("COG_category", paste(species_11[i], "_3"))
   assign (paste("comm", species_11[i], "1", sep = "_"), count_3)
   comm_all <- left_join(comm_all, count_3)
 }

rm(subset, count)
comm_all <- comm_all %>% replace(is.na(.), 0)
write.csv(comm_all, file = "comm_all.csv")

#by replicate, single data

I48_vs_community_short <- I48_vs_community_short %>% 
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)


subset_1 <- filter(I48_vs_community_short, species == "B. caecimuris I48")%>%select(13,17)%>%filter(!is.na(2))
count_1 <- count(subset_1, COG_category)
colnames(count_1) <- c("COG_category", "I48_1")
single_all <- count_1
subset_2 <- filter(I48_vs_community_short, species == "B. caecimuris I48")%>%select(13,18)%>%filter(!is.na(2))
count_2 <- count(subset_2, COG_category)
colnames(count_2) <- c("COG_category", "I48_2")
single_all <- left_join(single_all, count_2, by = c("COG_category"))
subset_3 <- filter(I48_vs_community_short, species == "B. caecimuris I48")%>%select(13,19)%>%filter(!is.na(2))
count_3 <- count(subset_3, COG_category)
colnames(count_3) <- c("COG_category", "I48_3")
single_all <- left_join(single_all, count_3)

single_vs_community_short <- single_vs_community_short %>% 
  separate_rows(COG_category, sep="") %>% filter(COG_category != "") %>%
  mutate_at('COG_category', as.factor)

for (i in 17:49){
  subset <- single_vs_community_short%>%select(2,13,i)
  dataset <- colnames(subset)[3]
  strain <- str_replace(dataset, "_1", "")
  strain <- str_replace(strain, "_2", "")
  strain <- str_replace(strain, "_3", "")
  subset <- subset %>% filter(grepl(strain, species)) %>%drop_na(3)
  count <- count(subset, COG_category)
  colnames(count) <- c("COG_category", dataset)
  if (nrow(count)==0) {
    count <- data.frame(COG_category=single_all$COG_category, new = NA)
    colnames(count) <- c("COG_category", dataset) }
  single_all <- left_join(single_all, count, by = c("COG_category"))
}

single_all <- single_all %>% replace(is.na(.), 0)
write.csv(single_all, file = "single_all.csv")
