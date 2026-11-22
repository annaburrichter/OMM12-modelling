library(tidyverse)
library(readxl)

P200_01_reanalysis <- read_excel("P200_01_reanalysis.xlsx", sheet = "Expression")
P200_03_reanalysis <- read_excel("P200_03_reanalysis.xlsx", sheet = "Expression")
P200_07_reanalysis <- read_excel("P200_07_reanalysis_12-23_ohne germfree.xlsx", sheet = "Expression")

P200_03_counts <- list()
P200_03_intensities <- list()
for (i in 2:13){
sample<- select(P200_03_reanalysis, ID, all_of(i)) %>% 
                filter(!if_any(everything(), is.na))%>%
                separate_wider_delim(ID, delim="_", names = c("X", "locus_strain", "locus_gene"), too_many = "drop") %>%
                separate_wider_delim(locus_gene, delim="..", names = c("locus_gene"), too_many = "drop")
colname <- names(sample)[4]
count <- count(sample, locus_strain)
intensity <- sample %>% group_by(locus_strain) %>%
  summarise(sum_value = sum(10^.data[[colname]], na.rm = TRUE),
  .groups = "drop")
name<- names(P200_03_reanalysis)[i]
P200_03_counts[[i]] <- count
names(P200_03_counts)[i] <- paste(name, "count", sep="_")
P200_03_intensities[[i]] <- intensity
names(P200_03_intensities)[i] <- paste(name, "intensity", sep="_")
}
rm(sample, count, name, i, colname, intensity)

P200_07_counts <- list()
P200_07_intensities <- list()
for (i in 2:13){
  sample<- select(P200_07_reanalysis, ID, all_of(i)) %>% 
    filter(!if_any(everything(), is.na))%>%
    separate_wider_delim(ID, delim="_", names = c("X", "locus_strain", "locus_gene"), too_many = "drop", too_few = "align_end") %>%
    separate_wider_delim(locus_gene, delim="..", names = c("locus_gene"), too_many = "drop")
  colname <- names(sample)[4]
  count <- count(sample, locus_strain)
  intensity <- sample %>% group_by(locus_strain) %>%
    summarise(sum_value = sum(10^.data[[colname]], na.rm = TRUE),
              .groups = "drop")
  name<- names(P200_07_reanalysis)[i]
  P200_07_counts[[i]] <- count
  names(P200_07_counts)[i] <- paste(name, "count", sep="_")
  P200_07_intensities[[i]] <- intensity
  names(P200_07_intensities)[i] <- paste(name, "intensity", sep="_")
}
rm(sample, count, name, i, colname, intensity)

P200_01_counts <- list()
P200_01_intensities <- list()
for (i in 2:43){
  sample<- select(P200_01_reanalysis, ID, all_of(i)) %>% 
    filter(!if_any(everything(), is.na))%>%
    separate_wider_delim(ID, delim="_", names = c("X", "locus_strain", "locus_gene"), too_many = "drop") %>%
    separate_wider_delim(locus_gene, delim="..", names = c("locus_gene"), too_many = "drop")
  colname <- names(sample)[4]
  count <- count(sample, locus_strain)
  intensity <- sample %>% group_by(locus_strain) %>%
    summarise(sum_value = sum(10^.data[[colname]], na.rm = TRUE),
              .groups = "drop")
  name<- names(P200_01_reanalysis)[i]
  P200_01_counts[[i]] <- count
  names(P200_01_counts)[i] <- paste(name, "count", sep="_")
  P200_01_intensities[[i]] <- intensity
  names(P200_01_intensities)[i] <- paste(name, "intensity", sep="_")
}
rm(sample, count, name, i, colname, intensity)

P200_03_counts <- P200_03_counts[-1]
P200_03_count <- reduce(P200_03_counts, full_join, by ='locus_strain')
colnames(P200_03_count) <- c("locus_strain", names(P200_03_counts))
P200_03_intensities <- P200_03_intensities[-1]
P200_03_intensity <- reduce(P200_03_intensities, full_join, by ='locus_strain')
colnames(P200_03_intensity) <- c("locus_strain", names(P200_03_intensities))

P200_07_counts <- P200_07_counts[-1]
P200_07_count <- reduce(P200_07_counts, full_join, by ='locus_strain')
colnames(P200_07_count) <- c("locus_strain", names(P200_07_counts))
P200_07_count <- filter(P200_07_count, str_detect(locus_strain, "^I5Q"))
P200_07_intensities <- P200_07_intensities[-1]
P200_07_intensity <- reduce(P200_07_intensities, full_join, by ='locus_strain')
colnames(P200_07_intensity) <- c("locus_strain", names(P200_07_intensities))
P200_07_intensity <- filter(P200_07_intensity, str_detect(locus_strain, "^I5Q"))

P200_01_counts <- P200_01_counts[-1]
P200_01_count <- reduce(P200_01_counts, full_join, by ='locus_strain')
colnames(P200_01_count) <- c("locus_strain", names(P200_01_counts))
P200_01_intensities <- P200_01_intensities[-1]
P200_01_intensity <- reduce(P200_01_intensities, full_join, by ='locus_strain')
colnames(P200_01_intensity) <- c("locus_strain", names(P200_01_intensities))

write.csv(P200_07_count, "P200_07_count.csv")
                                     