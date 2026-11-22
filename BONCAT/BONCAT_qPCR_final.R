# load packages
library(tidyverse)
library(ggpubr)
library(rstatix)
#load presaved file of detection limits 
load("detection_limits_2025.rda")
#set working directory (ignore if working from a project)

# data cleanup ------------------------------------------------------------


#search all exported qPCR files (ending in Abs Quant)
files <- list.files(pattern = "Quant.txt")

# Extract data from .csv--------------------------------------------------

qPCR <- tibble()

for (i in 1:length(files)) {
  
file <- read_tsv(files[i]) %>% mutate(across(everything(), as.character))
qPCR <- bind_rows(qPCR, file)
}
rm(file)
names(qPCR) <- str_replace_all(names(qPCR), c(" " = "_")) #remove any spaces to avoid issues later on

# change column types for later analysis
qPCR$Gene_Name <- as.factor(qPCR$Gene_Name)
qPCR$Concentration <- as.numeric(qPCR$Concentration)
qPCR$Concentration_Mean <- as.numeric(qPCR$Concentration_Mean)
#remove standard, negative control data points marked as excluded ! this does not adjust the mean! check again when necessary for a sample!
qPCR <- qPCR %>% filter(Excluded == "Unchecked") %>% filter(Sample_Type == "Unknown") 
qPCR$Sample_Name <- as.factor(qPCR$Sample_Name)
#remove duplicate rows and keep concentration mean, sample name and gene name
qPCR <- qPCR %>% distinct(Sample_Name, Gene_Name, .keep_all =TRUE) %>% select(Sample_Name, Gene_Name, Concentration_Mean)

#calculate total
qPCR_total <- qPCR %>% group_by(Sample_Name)%>%
  summarise(Concentration_Mean = sum(Concentration_Mean, na.rm = TRUE))%>%
  mutate(Gene_Name = replicate (length(levels(Sample_Name)), "total"), .before = Concentration_Mean)
#add LOD samples
LOD <- tibble(
  x = c("LOD", "LOD", "LOD", "LOD", "LOD", "LOD", "LOD", "LOD", "LOD", "LOD", "LOD", "LOD"),
  y = detection_limits_2025$detection_limits_strains,
  z = detection_limits_2025$limits_2025)
names(LOD) <- c("Sample_Name", "Gene_Name", "Concentration_Mean")

qPCR <- bind_rows(qPCR, qPCR_total, LOD)


# relative plots (stacked bars) -------------------------------------------

#create dataframe of relative values
qPCR_wide <- pivot_wider(qPCR, names_from = Gene_Name, values_from = Concentration_Mean)
qPCR_relative <- data.frame(qPCR_wide$Sample_Name, qPCR_wide[-1]/qPCR_wide$total)
colnames(qPCR_relative)[colnames(qPCR_relative) == "qPCR_wide.Sample_Name"] <- "Sample_Name"
qPCR_relative <- qPCR_relative%>% select(Sample_Name, KB1, YL2, KB18, YL27, YL31, YL32, YL44, YL45, I46, I48, I49, YL58)%>% 
  filter(Sample_Name != "LOD") 

#load all strain designations and Anna w.'s colour scheme
load("all_probes.rda")
col_annaw <- load("Colorscheme_AnnaW_ordered.rda")

#transpose data and transform to tibble
qPCR_relative <- column_to_rownames(qPCR_relative, var = "Sample_Name")
qPCR_relative_long <- as.data.frame(t(qPCR_relative)) %>% rownames_to_column(var = "strain")%>% 
  as_tibble()%>% pivot_longer(!strain, names_to = "sample", values_to = "abundance")

#plot all samples

ggplot(qPCR_relative_long, aes(x = sample, y = abundance, fill=factor(strain, levels=all_probes))) + 
  geom_col() + scale_fill_manual(values = col_annaw_ordered) +
  ylab("relative abundance") + 
  theme(legend.title=element_blank(), axis.title.x =element_blank(), legend.position = "left", axis.text.x = element_text(angle=90))


# plot absolute values ----------------------------------------------------


#read DNA concentrations (same sample names as qPCR, including one row "LOD" with the lowest measured concentration and highest sample amount, 
#and a column "sample_amount" with the amount of sample in ml, mg or other reference unit)

DNA <- read_tsv("DNA_BONCAT.txt")
qPCR_absolute <- full_join(qPCR_wide, DNA, by ="Sample_Name")%>%
  mutate(across(KB1:total, ~.x*qPCR_absolute$`total amount of DNA extracted (ng)`*100000/(qPCR_absolute$`final amount used in qPCR reaction`*qPCR_absolute$`number of cells`)))

#absolute plot
  
#separate sample names to be able to make groups (change sample name components how it fits your samples)
qPCR_absolute_long <- qPCR_absolute %>% 
  select(Sample_Name, KB1, YL2, KB18, YL27, YL31, YL32, YL44, YL45, I46, I48, I49, YL58, total)%>%
  pivot_longer(!Sample_Name, names_to = "probe", values_to = "abundance")%>%
  separate_wider_delim("Sample_Name", delim ="_", names=c("community", "BONCAT", "carbon_source"), too_few = "align_start")%>%
  separate_wider_regex(community, patterns = c(community = ".*(?=.$)", replicate = ".$"))

#reorder for plot
qPCR_absolute_long$probe <- factor(qPCR_absolute_long$probe,
                                   levels=c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58", "total"))

ttests <- qPCR_absolute_long %>%
  pivot_wider(names_from = BONCAT,values_from = abundance) %>%
  group_by(community, carbon_source, probe) %>%
  group_modify(~{dat <- .x %>%
      filter(!is.na(`B+`),!is.na(`B-`))
    # Need at least 2 complete pairs for a paired t-test
    if (nrow(dat) < 2) {return(tibble(n_pairs = nrow(dat),p = NA_real_))}
    tibble( n_pairs = nrow(dat),p = t.test(log10(dat$`B+`),log10(dat$`B-`),paired = TRUE)$p.value)}) %>%
  ungroup() %>%mutate(
    p.signif = case_when(
      is.na(p) ~ NA_character_,
      p <= 0.001 ~ "***",
      p <= 0.01  ~ "**",
      p <= 0.05  ~ "*",
      TRUE       ~ "ns" ))

label_pos <- qPCR_absolute_long %>%
  group_by(community, carbon_source, probe) %>%
  summarise(y.position = max(abundance, na.rm = TRUE) * 10000,.groups = "drop")

ttests <- left_join(ttests, label_pos,by = c("community", "carbon_source", "probe"))

#plot with statistics(change axis limits to 1e+6, 1e+11 to zoom in) (choose which data to plot using the subset function in the second line)


plot_tre <- 
  ggplot(qPCR_absolute_long, aes(y = abundance, x =community, color=BONCAT)) +
  geom_boxplot(data=subset(qPCR_absolute_long, community!="LO" & carbon_source == "Tre"),aes(), outlier.shape=NA, alpha=0.5) +
  geom_hline(data=subset(qPCR_absolute_long, community=="LO"), aes(yintercept = abundance), linetype =2) +
  geom_text(data = subset(ttests, carbon_source == "Tre"),aes(x = community,y = y.position,label = p.signif),inherit.aes = FALSE,color = "black",size = 4)+
  facet_wrap(~probe)+ 
  theme_bw()+
  scale_y_log10(limits=c(0.1,1e+16))+
  theme(axis.text.x = element_text(angle = 90))+ylab("copies/100 000 cells")+ggtitle("Trehalose")
print(plot_tre)

plot_pyr <- 
  ggplot(qPCR_absolute_long, aes(y = abundance, x =community, color=BONCAT)) +
  geom_boxplot(data=subset(qPCR_absolute_long, community!="LO" & carbon_source == "Pyr"),aes(), outlier.shape=NA, alpha=0.5) +
  geom_hline(data=subset(qPCR_absolute_long, community=="LO"), aes(yintercept = abundance), linetype =2) +
  geom_text(data = subset(ttests, carbon_source == "Pyr"),aes(x = community,y = y.position,label = p.signif),inherit.aes = FALSE,color = "black",size = 4)+
  facet_wrap(~probe)+ 
  theme_bw()+
  scale_y_log10(limits=c(0.1,1e+16))+
  theme(axis.text.x = element_text(angle = 90))+ylab("copies/100 000 cells")+ggtitle("Pyruvate")
print(plot_pyr)

plot_control <- 
  ggplot(qPCR_absolute_long, aes(y = abundance, x =community, color=BONCAT)) +
  geom_boxplot(data=subset(qPCR_absolute_long, community!="LO" & carbon_source == "-C"),aes(), outlier.shape=NA, alpha=0.5) +
  geom_hline(data=subset(qPCR_absolute_long, community=="LO"), aes(yintercept = abundance), linetype =2) +
  geom_text(data = subset(ttests, carbon_source == "-C"),aes(x = community,y = y.position,label = p.signif),inherit.aes = FALSE,color = "black",size = 4)+
  facet_wrap(~probe)+ 
  theme_bw()+
  scale_y_log10(limits=c(0.1,1e+16))+
  theme(axis.text.x = element_text(angle = 90))+ylab("copies/100 000 cells")+ggtitle("no carbon source")
print(plot_control)


# save data ---------------------------------------------------------------

#write relative and absolute data to csv
write.table(qPCR, file = "qPCR_base.csv", row.names= FALSE)
write.table(qPCR_absolute, file = "qPCR_absolute.csv", row.names= FALSE)
write.table(qPCR_relative, file = "qPCR_relative.csv", row.names= FALSE)


