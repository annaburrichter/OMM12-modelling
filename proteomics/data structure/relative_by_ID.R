library(tidyverse)
library(data.table)

relative_ID <- select(P200_01_count, locus_strain:batch_3_count) %>% 
  mutate(strain = c("I48", "I49", "YL45", "KB18", "YL32", "YL31", "YL44", "YL58", "I46", "YL2", "KB1", "YL27"))%>%
  bind_rows(summarise(., across(where(is.numeric), sum),across(where(is.character), ~'Total')))%>%
  select(batch_1_count:strain)
Sample_Name <- colnames(relative_ID)

relative_ID <- data.frame(relative_ID) %>%
  data.table::transpose()
colnames(relative_ID) <- relative_ID[4,]
relative_ID <- relative_ID %>% mutate(Sample_Name = Sample_Name) %>% slice(1:3)%>% 
  select(Sample_Name, KB1, YL2, KB18, YL27, YL31, YL32, YL44, YL45, I46, I48, I49, YL58)
relative_ID$Sample_Name <- gsub("batch", "community", relative_ID$Sample_Name)
relative_ID$Sample_Name <- gsub("_count", "", relative_ID$Sample_Name)

relative_ID_vivo <- select(P200_07_count, locus_strain:mouse_3_count) %>% slice(1:11)%>%
  mutate(strain = c("I48", "I49", "YL45", "KB18", "YL32", "YL31", "YL44", "YL58", "I46", "KB1", "YL27"))%>%
  bind_rows(summarise(., across(where(is.numeric), sum),across(where(is.character), ~'Total')))%>%
  select(AAM_1_count:strain)
Sample_Name <- colnames(relative_ID_vivo)

relative_ID_vivo <- data.frame(relative_ID_vivo) %>%
  data.table::transpose()
colnames(relative_ID_vivo) <- relative_ID_vivo[10,]
relative_ID_vivo <- relative_ID_vivo %>% mutate(Sample_Name = Sample_Name) %>% slice(1:9)%>%
  select(Sample_Name, KB1,  KB18, YL27, YL31, YL32, YL44, YL45, I46, I48, I49, YL58)
relative_ID_vivo$Sample_Name <- gsub("AAM", "AF", relative_ID_vivo$Sample_Name)
relative_ID_vivo$Sample_Name <- gsub("M10", "APF", relative_ID_vivo$Sample_Name)
relative_ID_vivo$Sample_Name <- gsub("_count", "", relative_ID_vivo$Sample_Name)

rel_ID <- relative_ID %>% mutate(across(c(KB1:YL58), as.numeric))%>%
                          rowwise()%>%
                          mutate(total=sum(c_across(KB1:YL58)),across(KB1:YL58, ~.x/total))%>% select(-total)

rel_ID_vivo <- relative_ID_vivo %>% mutate(across(c(KB1:YL58), as.numeric))%>%
  rowwise()%>%
  mutate(total=sum(c_across(KB1:YL58), na.rm=TRUE),across(KB1:YL58, ~.x/total))%>% select(-total)

  
relative_qPCR_single <- read.delim("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/data structure/qPCR_single_vs_comm.txt")
relative_qPCR_single$Sample_Name <- gsub("average 24h", "community", relative_qPCR_single$Sample_Name)
relative_qPCR <- read.delim("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/qPCR/20220725_data_relative_relevant.txt")


##plot all 

load("L:/AG_Stecher/!Anna Burrichter/R/all_probes.rda")
load("L:/AG_Stecher/!Anna Burrichter/R/Colorscheme_AnnaW_ordered.rda")

#single vs. comm by protein ID
samples <- rel_ID$Sample_Name
sample <- c("strain", samples)
row.names(rel_ID) <- rel_ID$Sample_Name
rel_ID <- select(rel_ID,-Sample_Name)
rel_ID_transposed <- as.data.frame(t(rel_ID))
rel_ID_transposed <-rownames_to_column(rel_ID_transposed)
colnames(rel_ID_transposed)[colnames(rel_ID_transposed) == "rowname"] <- "strain"
data_tibble <- as_tibble(rel_ID_transposed)
colnames(data_tibble) <- sample
data_tibble_long_rel_ID <- data_tibble %>% pivot_longer(2:length(data_tibble), names_to = "sample", values_to = "abundance")


#plot data as stacked columns, colour scheme by Anna W., strains ordered numerical
ggplot(data_tibble_long_rel_ID, aes(x = sample, y = abundance, fill=factor(strain, levels=all_probes))) + 
  geom_col() + scale_fill_manual(values = col_annaw_ordered) + theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", axis.text.x = element_text(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  labs(title = "relative abundance by protein count")

#single vs. comm by qPCR
samples <- relative_qPCR_single$Sample_Name
sample <- c("strain", samples)
row.names(relative_qPCR_single) <- relative_qPCR_single$Sample_Name
relative_qPCR_single <- select(relative_qPCR_single,-Sample_Name)
relative_qPCR_single_transposed <- as.data.frame(t(relative_qPCR_single))
relative_qPCR_single_transposed <-rownames_to_column(relative_qPCR_single_transposed)
colnames(relative_qPCR_single_transposed)[colnames(relative_qPCR_single_transposed) == "rowname"] <- "strain"
data_tibble <- as_tibble(relative_qPCR_single_transposed)
colnames(data_tibble) <- sample
data_tibble_long_relative_qPCR_single <- data_tibble %>% pivot_longer(2:length(data_tibble), names_to = "sample", values_to = "abundance")


#plot data as stacked columns, colour scheme by Anna W., strains ordered numerical
ggplot(data_tibble_long_relative_qPCR_single, aes(x = sample, y = abundance, fill=factor(strain, levels=all_probes))) + 
  geom_col() + scale_fill_manual(values = col_annaw_ordered) + theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", axis.text.x = element_text(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  labs(title = "relative abundance by qPCR")

#vitro vs. vivo by protein ID
rel_ID_vivo <- rel_ID_vivo %>% mutate(YL2 = 0)%>% relocate(YL2, .before = KB18)
samples <- rel_ID_vivo$Sample_Name
sample <- c("strain", samples)
row.names(rel_ID_vivo) <- rel_ID_vivo$Sample_Name
rel_ID_vivo <- select(rel_ID_vivo,-Sample_Name)
rel_ID_vivo_transposed <- as.data.frame(t(rel_ID_vivo))
rel_ID_vivo_transposed <-rownames_to_column(rel_ID_vivo_transposed)
colnames(rel_ID_vivo_transposed)[colnames(rel_ID_vivo_transposed) == "rowname"] <- "strain"
data_tibble <- as_tibble(rel_ID_vivo_transposed)
colnames(data_tibble) <- sample
data_tibble_long_rel_ID_vivo <- data_tibble %>% pivot_longer(2:length(data_tibble), names_to = "sample", values_to = "abundance")


#plot data as stacked columns, colour scheme by Anna W., strains ordered numerical
ggplot(data_tibble_long_rel_ID_vivo, aes(x = sample, y = abundance, fill=factor(strain, levels=all_probes))) + 
  geom_col() + scale_fill_manual(values = col_annaw_ordered) + theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", axis.text.x = element_text(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  labs(title = "relative abundance by protein count")

#vitro vs. vivo by qPCR
samples <- relative_qPCR$Sample_Name
sample <- c("strain", samples)
row.names(relative_qPCR) <- relative_qPCR$Sample_Name
relative_qPCR <- select(relative_qPCR,-Sample_Name)
relative_qPCR_transposed <- as.data.frame(t(relative_qPCR))
relative_qPCR_transposed <-rownames_to_column(relative_qPCR_transposed)
colnames(relative_qPCR_transposed)[colnames(relative_qPCR_transposed) == "rowname"] <- "strain"
data_tibble <- as_tibble(relative_qPCR_transposed)
colnames(data_tibble) <- sample
data_tibble_long_relative_qPCR <- data_tibble %>% pivot_longer(2:length(data_tibble), names_to = "sample", values_to = "abundance")


#plot data as stacked columns, colour scheme by Anna W., strains ordered numerical
ggplot(data_tibble_long_relative_qPCR, aes(x = sample, y = abundance, fill=factor(strain, levels=all_probes))) + 
  geom_col() + scale_fill_manual(values = col_annaw_ordered) + theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", axis.text.x = element_text(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  labs(title = "relative abundance by qPCR")


##correlation protein count and qPCR

#community vs. single
corr_single <- left_join(data_tibble_long_rel_ID, data_tibble_long_relative_qPCR_single, by =c('strain' = 'strain', 'sample' = 'sample'))
ggplot(corr_single, aes(x = abundance.x, y = abundance.y, color =factor(strain, levels=all_probes))) + 
  geom_point() + scale_color_manual(values = col_annaw_ordered) + 
  xlim(0, 0.65) + ylim(0,0.65)+
  geom_abline(slope=1, intercept=0, linetype = 'dashed', colour = 'darkgrey')+
  theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "community") + xlab("relative abundance by protein count") + ylab("relative abundance by qPCR")

#in vitro vs. in vivo
corr_vivo <- right_join(data_tibble_long_rel_ID_vivo, data_tibble_long_relative_qPCR, by =c('strain' = 'strain', 'sample' = 'sample'))
corr_AF <- filter(corr_vivo, str_detect(sample, "^AF"))
corr_APF <- filter(corr_vivo, str_detect(sample, "^APF"))
corr_mouse <- filter(corr_vivo, str_detect(sample, "^mouse"))

ggplot(corr_AF, aes(x = abundance.x, y = abundance.y, color =factor(strain, levels=all_probes))) + 
  geom_point() + scale_color_manual(values = col_annaw_ordered) + 
  xlim(0, 0.65) + ylim(0,0.65)+
  geom_abline(slope=1, intercept=0, linetype = 'dashed', colour = 'darkgrey')+
  theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "community in AF medium") + xlab("relative abundance by protein count") + ylab("relative abundance by qPCR")

ggplot(corr_APF, aes(x = abundance.x, y = abundance.y, color =factor(strain, levels=all_probes))) + 
  geom_point() + scale_color_manual(values = col_annaw_ordered) + 
  xlim(0, 0.65) + ylim(0,0.65)+
  geom_abline(slope=1, intercept=0, linetype = 'dashed', colour = 'darkgrey')+
  theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "community in APF medium") + xlab("relative abundance by protein count") + ylab("relative abundance by qPCR")


ggplot(corr_mouse, aes(x = abundance.x, y = abundance.y, color =factor(strain, levels=all_probes))) + 
  geom_point() + scale_color_manual(values = col_annaw_ordered) + 
  xlim(0, 0.65) + ylim(0,0.65)+
  geom_abline(slope=1, intercept=0, linetype = 'dashed', colour = 'darkgrey')+
  theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "community in mouse cecum") + xlab("relative abundance by protein count") + ylab("relative abundance by qPCR")

