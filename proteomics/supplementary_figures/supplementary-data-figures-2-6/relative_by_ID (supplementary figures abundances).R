library(tidyverse)
library(data.table)
locus_start <- read_delim("locustag_start.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)


##relative community composition by counts
counts_community <- select(P200_01_count, locus_strain:batch_3_count) %>% 
  right_join(locus_start, by = join_by(locus_strain==locus_tag))%>%
  bind_rows(summarise(., across(where(is.numeric), sum),across(where(is.character), ~'Total')))%>%
  select(batch_1_count:strain)
Sample_Name <- colnames(counts_community)

counts_community <- data.frame(counts_community) %>%
  data.table::transpose()
colnames(counts_community) <- counts_community[4,]
counts_community <- counts_community %>% mutate(Sample_Name = Sample_Name) %>% slice(1:3)%>% 
  select(Sample_Name, KB1, YL2, KB18, YL27, YL31, YL32, YL44, YL45, I46, I48, I49, YL58)
counts_community$Sample_Name <- gsub("batch", "community", counts_community$Sample_Name)
counts_community$Sample_Name <- gsub("_count", "", counts_community$Sample_Name)

counts_community_vivo <- select(P200_07_count, locus_strain:mouse_3_count) %>% slice(1:11)%>%
  right_join(locus_start, by = join_by(locus_strain==locus_tag))%>%
  bind_rows(summarise(., across(where(is.numeric), sum),across(where(is.character), ~'Total')))%>%
  select(AAM_1_count:strain)
Sample_Name <- colnames(counts_community_vivo)

counts_community_vivo <- data.frame(counts_community_vivo) %>%
  data.table::transpose()
colnames(counts_community_vivo) <- counts_community_vivo[10,]
counts_community_vivo <- counts_community_vivo %>% mutate(Sample_Name = Sample_Name) %>% slice(1:9)%>%
  select(Sample_Name, KB1,  KB18, YL27, YL31, YL32, YL44, YL45, I46, I48, I49, YL58)
counts_community_vivo$Sample_Name <- gsub("AAM", "AF", counts_community_vivo$Sample_Name)
counts_community_vivo$Sample_Name <- gsub("M10", "APF", counts_community_vivo$Sample_Name)
counts_community_vivo$Sample_Name <- gsub("_count", "", counts_community_vivo$Sample_Name)

relative_counts <- counts_community %>% mutate(across(c(KB1:YL58), as.numeric))%>%
                          rowwise()%>%
                          mutate(total=sum(c_across(KB1:YL58)),across(KB1:YL58, ~.x/total))%>% select(-total)

relative_counts_vivo <- counts_community_vivo %>% mutate(across(c(KB1:YL58), as.numeric))%>%
  rowwise()%>%
  mutate(total=sum(c_across(KB1:YL58), na.rm=TRUE),across(KB1:YL58, ~.x/total))%>% select(-total)


##relative community composition by intensities

intensity_community <- select(P200_01_intensity, locus_strain:batch_3_intensity) %>% 
  right_join(locus_start, by = join_by(locus_strain==locus_tag))%>%
  bind_rows(summarise(., across(where(is.numeric), sum),across(where(is.character), ~'Total')))%>%
  select(batch_1_intensity:strain)
Sample_Name <- colnames(intensity_community)

intensity_community <- data.frame(intensity_community) %>%
  data.table::transpose()
colnames(intensity_community) <- intensity_community[4,]
intensity_community <- intensity_community %>% mutate(Sample_Name = Sample_Name) %>% slice(1:3)%>% 
  select(Sample_Name, KB1, YL2, KB18, YL27, YL31, YL32, YL44, YL45, I46, I48, I49, YL58)
intensity_community$Sample_Name <- gsub("batch", "community", intensity_community$Sample_Name)
intensity_community$Sample_Name <- gsub("_intensity", "", intensity_community$Sample_Name)

intensity_community_vivo <- select(P200_07_intensity, locus_strain:mouse_3_intensity) %>% slice(1:11)%>%
  right_join(locus_start, by = join_by(locus_strain==locus_tag))%>%
  bind_rows(summarise(., across(where(is.numeric), sum),across(where(is.character), ~'Total')))%>%
  select(AAM_1_intensity:strain)
Sample_Name <- colnames(intensity_community_vivo)

intensity_community_vivo <- data.frame(intensity_community_vivo) %>%
  data.table::transpose()
colnames(intensity_community_vivo) <- intensity_community_vivo[10,]
intensity_community_vivo <- intensity_community_vivo %>% mutate(Sample_Name = Sample_Name) %>% slice(1:9)%>%
  select(Sample_Name, KB1,  KB18, YL27, YL31, YL32, YL44, YL45, I46, I48, I49, YL58)
intensity_community_vivo$Sample_Name <- gsub("AAM", "AF", intensity_community_vivo$Sample_Name)
intensity_community_vivo$Sample_Name <- gsub("M10", "APF", intensity_community_vivo$Sample_Name)
intensity_community_vivo$Sample_Name <- gsub("_intensity", "", intensity_community_vivo$Sample_Name)

relative_intensity <- intensity_community %>% mutate(across(c(KB1:YL58), as.numeric))%>%
  rowwise()%>%
  mutate(total=sum(c_across(KB1:YL58)),across(KB1:YL58, ~.x/total))%>% select(-total)

relative_intensity_vivo <- intensity_community_vivo %>% mutate(across(c(KB1:YL58), as.numeric))%>%
  rowwise()%>%
  mutate(total=sum(c_across(KB1:YL58), na.rm=TRUE),across(KB1:YL58, ~.x/total))%>% select(-total)

##relative community composition by intensities
  
relative_qPCR_single <- read.delim("qPCR_single_vs_comm.txt")
relative_qPCR_single$Sample_Name <- gsub("average 24h", "community", relative_qPCR_single$Sample_Name)
relative_qPCR <- read.delim("qPCR/20220725_data_relative_relevant.txt")


##plot all 

load("L:/Members/!Anna Burrichter/R/all_probes.rda")
load("L:/Members/!Anna Burrichter/R/Colorscheme_AnnaW_ordered.rda")

#single vs. comm by protein count
samples <- relative_counts$Sample_Name
sample <- c("strain", samples)
row.names(relative_counts) <- relative_counts$Sample_Name
relative_counts <- select(relative_counts,-Sample_Name)
relative_counts_transposed <- as.data.frame(t(relative_counts))
relative_counts_transposed <-rownames_to_column(relative_counts_transposed)
colnames(relative_counts_transposed)[colnames(relative_counts_transposed) == "rowname"] <- "strain"
data_tibble <- as_tibble(relative_counts_transposed)
colnames(data_tibble) <- sample
data_tibble_long_relative_counts <- data_tibble %>% pivot_longer(2:length(data_tibble), names_to = "sample", values_to = "abundance")


#plot data as stacked columns, colour scheme by Anna W., strains ordered numerical
ggplot(data_tibble_long_relative_counts, aes(x = sample, y = abundance, fill=factor(strain, levels=all_probes))) + 
  geom_col() + scale_fill_manual(values = col_annaw_ordered) + theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", axis.text.x = element_text(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  labs(title = "relative abundance by protein count")

#single vs. comm by protein intensity
samples <- relative_intensity$Sample_Name
sample <- c("strain", samples)
row.names(relative_intensity) <- relative_intensity$Sample_Name
relative_intensity <- select(relative_intensity,-Sample_Name)
relative_intensity_transposed <- as.data.frame(t(relative_intensity))
relative_intensity_transposed <-rownames_to_column(relative_intensity_transposed)
colnames(relative_intensity_transposed)[colnames(relative_intensity_transposed) == "rowname"] <- "strain"
data_tibble <- as_tibble(relative_intensity_transposed)
colnames(data_tibble) <- sample
data_tibble_long_relative_intensity <- data_tibble %>% pivot_longer(2:length(data_tibble), names_to = "sample", values_to = "abundance")


#plot data as stacked columns, colour scheme by Anna W., strains ordered numerical
ggplot(data_tibble_long_relative_intensity, aes(x = sample, y = abundance, fill=factor(strain, levels=all_probes))) + 
  geom_col() + scale_fill_manual(values = col_annaw_ordered) + theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", axis.text.x = element_text(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  labs(title = "relative abundance by protein intensity")

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
  theme(legend.title=element_blank(), legend.position = "left", axis.text.x = element_text(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  labs(title = "relative abundance by qPCR")

#vitro vs. vivo by protein intensity
relative_counts_vivo <- relative_counts_vivo %>% mutate(YL2 = 0)%>% relocate(YL2, .before = KB18)
samples <- relative_counts_vivo$Sample_Name
sample <- c("strain", samples)
row.names(relative_counts_vivo) <- relative_counts_vivo$Sample_Name
relative_counts_vivo <- select(relative_counts_vivo,-Sample_Name)
relative_counts_vivo_transposed <- as.data.frame(t(relative_counts_vivo))
relative_counts_vivo_transposed <-rownames_to_column(relative_counts_vivo_transposed)
colnames(relative_counts_vivo_transposed)[colnames(relative_counts_vivo_transposed) == "rowname"] <- "strain"
data_tibble <- as_tibble(relative_counts_vivo_transposed)
colnames(data_tibble) <- sample
data_tibble_long_relative_counts_vivo <- data_tibble %>% pivot_longer(2:length(data_tibble), names_to = "sample", values_to = "abundance")


#plot data as stacked columns, colour scheme by Anna W., strains ordered numerical
ggplot(data_tibble_long_relative_counts_vivo, aes(x = sample, y = abundance, fill=factor(strain, levels=all_probes))) + 
  geom_col() + scale_fill_manual(values = col_annaw_ordered) + theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", axis.text.x = element_text(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  labs(title = "relative abundance by protein intensity")

#vitro vs. vivo by protein ID
relative_intensity_vivo <- relative_intensity_vivo %>% mutate(YL2 = 0)%>% relocate(YL2, .before = KB18)
samples <- relative_intensity_vivo$Sample_Name
sample <- c("strain", samples)
row.names(relative_intensity_vivo) <- relative_intensity_vivo$Sample_Name
relative_intensity_vivo <- select(relative_intensity_vivo,-Sample_Name)
relative_intensity_vivo_transposed <- as.data.frame(t(relative_intensity_vivo))
relative_intensity_vivo_transposed <-rownames_to_column(relative_intensity_vivo_transposed)
colnames(relative_intensity_vivo_transposed)[colnames(relative_intensity_vivo_transposed) == "rowname"] <- "strain"
data_tibble <- as_tibble(relative_intensity_vivo_transposed)
colnames(data_tibble) <- sample
data_tibble_long_relative_intensity_vivo <- data_tibble %>% pivot_longer(2:length(data_tibble), names_to = "sample", values_to = "abundance")


#plot data as stacked columns, colour scheme by Anna W., strains ordered numerical
ggplot(data_tibble_long_relative_intensity_vivo, aes(x = sample, y = abundance, fill=factor(strain, levels=all_probes))) + 
  geom_col() + scale_fill_manual(values = col_annaw_ordered) + theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", axis.text.x = element_text(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
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
  theme(legend.title=element_blank(), legend.position = "left", axis.text.x = element_text(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  labs(title = "relative abundance by qPCR")


##correlation protein count and qPCR

#community vs. single
corr_single <- left_join(data_tibble_long_relative_counts, data_tibble_long_relative_qPCR_single, by =c('strain' = 'strain', 'sample' = 'sample'))
ggplot(corr_single, aes(x = abundance.x, y = abundance.y, color =factor(strain, levels=all_probes))) + 
  geom_point() + scale_color_manual(values = col_annaw_ordered) + 
  xlim(0, 0.65) + ylim(0,0.65)+
  geom_abline(slope=1, intercept=0, linetype = 'dashed', colour = 'darkgrey')+
  theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "community") + xlab("relative abundance by protein count") + ylab("relative abundance by qPCR")

#in vitro vs. in vivo
corr_vivo <- right_join(data_tibble_long_relative_counts_vivo, data_tibble_long_relative_qPCR, by =c('strain' = 'strain', 'sample' = 'sample'))
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

##correlation protein intensity and qPCR

#community vs. single
corr_single <- left_join(data_tibble_long_relative_intensity, data_tibble_long_relative_qPCR_single, by =c('strain' = 'strain', 'sample' = 'sample'))
ggplot(corr_single, aes(x = abundance.x, y = abundance.y, color =factor(strain, levels=all_probes))) + 
  geom_point() + scale_color_manual(values = col_annaw_ordered) + 
  xlim(0, 0.65) + ylim(0,0.65)+
  geom_abline(slope=1, intercept=0, linetype = 'dashed', colour = 'darkgrey')+
  theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "community") + xlab("relative abundance by protein intensity") + ylab("relative abundance by qPCR")

#in vitro vs. in vivo
corr_vivo <- right_join(data_tibble_long_relative_intensity_vivo, data_tibble_long_relative_qPCR, by =c('strain' = 'strain', 'sample' = 'sample'))
corr_AF <- filter(corr_vivo, str_detect(sample, "^AF"))
corr_APF <- filter(corr_vivo, str_detect(sample, "^APF"))
corr_mouse <- filter(corr_vivo, str_detect(sample, "^mouse"))

ggplot(corr_AF, aes(x = abundance.x, y = abundance.y, color =factor(strain, levels=all_probes))) + 
  geom_point() + scale_color_manual(values = col_annaw_ordered) + 
  xlim(0, 0.65) + ylim(0,0.65)+
  geom_abline(slope=1, intercept=0, linetype = 'dashed', colour = 'darkgrey')+
  theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "community in AF medium") + xlab("relative abundance by protein intensity") + ylab("relative abundance by qPCR")

ggplot(corr_APF, aes(x = abundance.x, y = abundance.y, color =factor(strain, levels=all_probes))) + 
  geom_point() + scale_color_manual(values = col_annaw_ordered) + 
  xlim(0, 0.65) + ylim(0,0.65)+
  geom_abline(slope=1, intercept=0, linetype = 'dashed', colour = 'darkgrey')+
  theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "community in APF medium") + xlab("relative abundance by protein intensity") + ylab("relative abundance by qPCR")


ggplot(corr_mouse, aes(x = abundance.x, y = abundance.y, color =factor(strain, levels=all_probes))) + 
  geom_point() + scale_color_manual(values = col_annaw_ordered) + 
  xlim(0, 0.65) + ylim(0,0.65)+
  geom_abline(slope=1, intercept=0, linetype = 'dashed', colour = 'darkgrey')+
  theme_bw()+
  theme(legend.title=element_blank(), legend.position = "left", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "community in mouse cecum") + xlab("relative abundance by protein intensity") + ylab("relative abundance by qPCR")

