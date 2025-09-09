library(tidyverse)
library(dplyr)
library(readxl)
library(ggplot2)
setwd("~/Anna/20231208")
# read data from working directory
annotations <- read.delim("protein_annotation.tsv")
P200_07_expression <- read_excel("reanalysed_mouse_vs_invitro.xlsx", 
                                 sheet = "Tabelle1") %>% rename(locus_tag = "locus tag_old")

#join annotations and proteome data
prot_ann <- left_join(P200_07_expression, annotations, by = "locus_tag")

#filter out useful columns
prot_short <- prot_ann %>% select(locus_tag, org, COG_category, AAM_1, AAM_2, AAM_3, M10_1, M10_2, M10_3, mouse_1, mouse_2, mouse_3)
prot_short$COG_category[is.na(prot_short$COG_category)] <- " "                  

#separate rows with more than one COG category, keep only first one
prot_unique_COG <- prot_short %>% separate_longer_position(COG_category, 1, keep_empty = TRUE) %>% distinct(locus_tag, .keep_all = TRUE)

#filter by conditions
unique_COG_AAM <- prot_unique_COG %>% select (locus_tag, org, COG_category, AAM_1, AAM_2, AAM_3) %>% filter((AAM_1!="NA" & AAM_2!="NA") | (AAM_1!="NA" & AAM_3!="NA") | (AAM_2!="NA" & AAM_3!="NA")) 
unique_COG_M10 <- prot_unique_COG %>% select (locus_tag, org, COG_category, M10_1, M10_2, M10_3) %>% filter((M10_1!="NA" & M10_2!="NA") | (M10_1!="NA" & M10_3!="NA") | (M10_2!="NA" & M10_3!="NA"))
unique_COG_mouse <- prot_unique_COG %>% select (locus_tag, org, COG_category, mouse_1, mouse_2, mouse_3) %>% filter((mouse_1!="NA" & mouse_2!="NA") | (mouse_1!="NA" & mouse_3!="NA") | (mouse_2!="NA" & mouse_3!="NA"))
unique_COG_AAM$COG_category <- as.factor(unique_COG_AAM$COG_category)
unique_COG_M10$COG_category <- as.factor(unique_COG_M10$COG_category)
unique_COG_mouse$COG_category <- as.factor(unique_COG_mouse$COG_category)
categories_AAM <- count(unique_COG_AAM, COG_category)
categories_M10 <- count(unique_COG_M10, COG_category)
categories_mouse <- count(unique_COG_mouse, COG_category)


#plot distribution
ggplot(categories_AAM, aes(x="", y = n, fill=COG_category)) + 
  geom_col(position = "stack", color = "black") + 
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_viridis_d() + 
  geom_text (aes(label = COG_category), position=position_stack((vjust=0.5))) + 
  ggtitle("AAM")

ggplot(categories_M10, aes(x="", y = n, fill=COG_category)) + 
  geom_col(position = "stack", color = "black") + 
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_viridis_d() + 
  geom_text (aes(label = COG_category), position=position_stack((vjust=0.5))) + 
  ggtitle("M10")

ggplot(categories_mouse, aes(x="", y = n, fill=COG_category)) + 
  geom_col(position = "stack", color = "black") + 
  coord_polar(theta = "y") + 
  theme_void() +
  scale_fill_viridis_d() + 
  geom_text (aes(label = COG_category), position=position_stack((vjust=0.5))) + 
  ggtitle("mouse")
                         
write.csv(categories_AAM, file = "COG_AAM.csv")     
write.csv(categories_M10, file = "COG_M10.csv")  
write.csv(categories_mouse, file = "COG_mouse.csv")  

##count separate replicates 

COGs_AAM_1 <- unique_COG_AAM %>%
              select(locus_tag, org, COG_category, AAM_1)%>%
              filter(AAM_1 != "NA")
categories_AAM_1 <- count(COGs_AAM_1, COG_category) %>% rename("AF_1" = "n")

COGs_AAM_2 <- unique_COG_AAM %>%
  select(locus_tag, org, COG_category, AAM_2)%>%
  filter(AAM_2 != "NA")
categories_AAM_2 <- count(COGs_AAM_2, COG_category) %>% rename("AF_2" = "n")


COGs_AAM_3 <- unique_COG_AAM %>%
  select(locus_tag, org, COG_category, AAM_3)%>%
  filter(AAM_3 != "NA")
categories_AAM_3 <- count(COGs_AAM_3, COG_category) %>% rename("AF_3" = "n")

COGs_M10_1 <- unique_COG_M10 %>%
  select(locus_tag, org, COG_category, M10_1)%>%
  filter(M10_1 != "NA")
categories_M10_1 <- count(COGs_M10_1, COG_category) %>% rename("APF_1" = "n")

COGs_M10_2 <- unique_COG_M10 %>%
  select(locus_tag, org, COG_category, M10_2)%>%
  filter(M10_2 != "NA")
categories_M10_2 <- count(COGs_M10_2, COG_category) %>% rename("APF_2" = "n")


COGs_M10_3 <- unique_COG_M10 %>%
  select(locus_tag, org, COG_category, M10_3)%>%
  filter(M10_3 != "NA")
categories_M10_3 <- count(COGs_M10_3, COG_category) %>% rename("APF_3" = "n")

COGs_mouse_1 <- unique_COG_mouse %>%
  select(locus_tag, org, COG_category, mouse_1)%>%
  filter(mouse_1 != "NA")
categories_mouse_1 <- count(COGs_mouse_1, COG_category) %>% rename("mouse_1" = "n")

COGs_mouse_2 <- unique_COG_mouse %>%
  select(locus_tag, org, COG_category, mouse_2)%>%
  filter(mouse_2 != "NA")
categories_mouse_2 <- count(COGs_mouse_2, COG_category) %>% rename("mouse_2" = "n")


COGs_mouse_3 <- unique_COG_mouse %>%
  select(locus_tag, org, COG_category, mouse_3)%>%
  filter(mouse_3 != "NA")
categories_mouse_3 <- count(COGs_mouse_3, COG_category) %>% rename("mouse_3" = "n")

categories <- list(categories_AAM_1, categories_AAM_2, categories_AAM_3, categories_M10_1, categories_M10_2, 
                   categories_M10_3, categories_mouse_1, categories_mouse_2, categories_mouse_3) %>%
              reduce(full_join, by = "COG_category")
write.csv(categories, file = "categories_replicates.csv")

#species contribution per COG category
unique_COG_AAM$org <- as.factor(unique_COG_AAM$org)
unique_COG_M10$org <- as.factor(unique_COG_M10$org)
unique_COG_mouse$org <- as.factor(unique_COG_mouse$org)

COGs <- levels(categories_mouse$COG_category)

#AF
list_AF <- list()
for (i in 1:length(COGs)){

x <- unique_COG_AAM %>% select(2,3,4)%>%
                                dplyr::filter(AAM_1 != "NA")%>%
                                dplyr::filter(COG_category == COGs[i])
count_x<- count(x, org)
count_x<- rename(count_x, setNames("n", paste0("AF_1_", COGs[i])))

y <- unique_COG_AAM %>% select(2,3,5)%>%
  dplyr::filter(AAM_2 != "NA")%>%
  dplyr::filter(COG_category == COGs[i])
count_y<- count(y, org)
count_y<- rename(count_y, setNames("n", paste0("AF_2_", COGs[i])))

z <- unique_COG_AAM %>% select(2,3,6)%>%
  dplyr::filter(AAM_3 != "NA")%>%
  dplyr::filter(COG_category == COGs[i])
count_z<- count(z, org)
count_z<- rename(count_z, setNames("n", paste0("AF_3_", COGs[i])))

res <- full_join(count_x, count_y)
res <- full_join(res, count_z)
res <- pivot_longer(res, cols= !org, names_to = c("med", "rep", ".cat"), names_sep = "_", values_to = "count")%>%
  select(org, rep, .cat, count)%>%
  unite(org_rep, org, rep, sep = "_")%>%
  pivot_wider(names_from = .cat, values_from = count)
list_AF[[i]] <- res
rm(x, y, z, count_x, count_y, count_z, res)}

AF_species_contrib <- reduce(list_AF, full_join, by = "org_rep")%>%column_to_rownames("org_rep")%>%t()
  

#APF
list_APF <- list()
for (i in 1:length(COGs)){
  
  x <- unique_COG_M10 %>% select(2,3,4)%>%
    dplyr::filter(M10_1 != "NA")%>%
    dplyr::filter(COG_category == COGs[i])
  count_x<- count(x, org)
  count_x<- rename(count_x, setNames("n", paste0("APF_1_", COGs[i])))
  
  y <- unique_COG_M10 %>% select(2,3,5)%>%
    dplyr::filter(M10_2 != "NA")%>%
    dplyr::filter(COG_category == COGs[i])
  count_y<- count(y, org)
  count_y<- rename(count_y, setNames("n", paste0("APF_2_", COGs[i])))
  
  z <- unique_COG_M10 %>% select(2,3,6)%>%
    dplyr::filter(M10_3 != "NA")%>%
    dplyr::filter(COG_category == COGs[i])
  count_z<- count(z, org)
  count_z<- rename(count_z, setNames("n", paste0("APF_3_", COGs[i])))
  
  res <- full_join(count_x, count_y)
  res <- full_join(res, count_z)
  res <- pivot_longer(res, cols= !org, names_to = c("med", "rep", ".cat"), names_sep = "_", values_to = "count")%>%
    select(org, rep, .cat, count)%>%
    unite(org_rep, org, rep, sep = "_")%>%
    pivot_wider(names_from = .cat, values_from = count)
  
  list_APF[[i]] <- res
  rm(x, y, z, count_x, count_y, count_z, res)}

APF_species_contrib <- reduce(list_APF, full_join, by = "org_rep")%>%column_to_rownames("org_rep")%>%t()


#mouse
list_mouse <- list()
for (i in 1:length(COGs)){
  
  x <- unique_COG_mouse %>% select(2,3,4)%>%
    dplyr::filter(mouse_1 != "NA")%>%
    dplyr::filter(COG_category == COGs[i])
  count_x<- count(x, org)
  count_x<- rename(count_x, setNames("n", paste0("mouse_1_", COGs[i])))
  
  y <- unique_COG_mouse %>% select(2,3,5)%>%
    dplyr::filter(mouse_2 != "NA")%>%
    dplyr::filter(COG_category == COGs[i])
  count_y<- count(y, org)
  count_y<- rename(count_y, setNames("n", paste0("mouse_2_", COGs[i])))
  
  z <- unique_COG_mouse %>% select(2,3,6)%>%
    dplyr::filter(mouse_3 != "NA")%>%
    dplyr::filter(COG_category == COGs[i])
  count_z<- count(z, org)
  count_z<- rename(count_z, setNames("n", paste0("mouse_3_", COGs[i])))
  
  res <- full_join(count_x, count_y)
  res <- full_join(res, count_z)
  res <- pivot_longer(res, cols= !org, names_to = c("med", "rep", ".cat"), names_sep = "_", values_to = "count")%>%
    select(org, rep, .cat, count)%>%
    unite(org_rep, org, rep, sep = "_")%>%
    pivot_wider(names_from = .cat, values_from = count)
  
  list_mouse[[i]] <- res
  rm(x, y, z, count_x, count_y, count_z, res)}

mouse_species_contrib <- reduce(list_mouse, full_join, by = "org_rep")%>%column_to_rownames("org_rep")%>%t()


write.csv(AF_species_contrib, file = "AF_species_contrib.csv")
write.csv(APF_species_contrib, file = "APF_species_contrib.csv")
write.csv(mouse_species_contrib, file = "mouse_species_contrib.csv")

