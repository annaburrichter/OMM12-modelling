setwd("//nas.ads.mwn.de/lmmepki$/AG_Stecher/!Anna Burrichter/metabolic models AG Kaleta/proteomics data/qPCR")
library(tidyverse)
#import textfile corrected for sampled volume and DNA concentration so that output data is reads/ml culture (correction done in Excel)
proteomics_absolute <- read_tsv("absolute_ml.txt")
proteomics_absolute_rep <- as_tibble(proteomics_absolute)

#separate sample names to be able to make groups
rep_long <- proteomics_absolute_rep %>% pivot_longer(2:length(proteomics_absolute_rep), names_to = "probe", values_to = "abundance")
rep_long_sep <- separate(rep_long, "sample", sep="_", into=c("condition", "replicate"))

#remove LOD (also use to separate for separate plots)
LOD <- filter(rep_long_sep, condition=="LOD")
LOD$probe <- factor(LOD$probe,
                    levels=c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58", "E.coli", "total"))
#filter out
relevant <- filter(rep_long_sep, condition=="AF" | condition=="APF" |condition=="mouse")
relevant$probe <- factor(relevant$probe,
                    levels=c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58", "E.coli", "total"))
germfree <- filter(rep_long_sep, condition=="germfree")
germfree$probe <- factor(germfree$probe,
                         levels=c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58", "E.coli", "total"))

relevant_GF <- filter(rep_long_sep, condition=="AF" | condition=="APF" |condition=="mouse"|condition=="germfree")
relevant_GF$probe <- factor(relevant_GF$probe,
                         levels=c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58", "E.coli", "total"))
relevant_GF$condition <- factor(relevant_GF$condition,
                         levels=c("AF", "APF", "mouse", "germfree"))
relevant_GF$replicate <- as.numeric(relevant_GF$replicate)
relevant_GF$replicate <- replace_na(relevant_GF$replicate, 1)


#plot 
ggplot(relevant, aes(y = abundance, x =condition, color=condition)) +
  geom_boxplot(data=relevant,aes(), outlier.shape=NA, alpha=0.5) +
  scale_color_manual(values = c("#FFA040", "#FF6000", "#800080")) + 
  geom_hline(data=LOD, aes(yintercept = abundance), linetype =2) +
  theme_bw() + 
  #geom_jitter(aes(col=sample), size=1.5, width=0.2) +    
  facet_wrap(~probe, nrow = 3, ncol = 5)+ 
  scale_y_log10(limits=c(0.1,1e+16))+
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")+ylab("copies/1ml culture")


#plot germfree
ggplot(germfree, aes(y = abundance, x =condition, color=condition)) +
  geom_boxplot(data=germfree,aes(), outlier.shape=NA, alpha=0.5) +
  geom_hline(data=LOD, aes(yintercept = abundance), linetype =2) +
  #geom_jitter(aes(col=sample), size=1.5, width=0.2) +    
  theme_bw() + 
  facet_wrap(~probe, nrow = 3, ncol = 5)+ 
  scale_y_log10(limits=c(0.1,1e+16))+
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")+ylab("copies/1ml culture")

#plot all
ggplot(relevant_GF, aes(y = abundance, x =condition, color=condition)) +
  geom_boxplot(data=relevant_GF,aes(), outlier.shape=NA, alpha=0.5) +
  scale_color_manual(values = c("#FFA040", "#FF6000", "#800080", "steelblue")) + 
  geom_hline(data=LOD, aes(yintercept = abundance), linetype =2) +
  theme_bw() + 
  #geom_jitter(aes(col=sample), size=1.5, width=0.2) +    
  facet_wrap(~probe, nrow = 3, ncol = 5)+ 
  scale_y_log10(limits=c(0.1,1e+16))+
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")+ylab("copies/1ml culture")

