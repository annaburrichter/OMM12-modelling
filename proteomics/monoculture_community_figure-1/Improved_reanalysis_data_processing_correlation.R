library(openxlsx)
library(tidyverse)
library(readxl)

annotations <- read.delim2("protein_annotation.tsv")

##read data, replace I48 dataset with separate one (contamination of original)

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

#list of species

species <- c("E. faecalis KB1", "B. animalis YL2", "A. muris KB18", "M. intestinale YL27", "F. plautii YL31",
             "E. clostridioformis YL32", "A. muciniphila YL44", "T. muris YL45", "C. innocuum I46", 
             "B. caecimuris I48", "L. reuteri I49", "B. coccoides YL58")
species_short <- c("KB1", "YL2", "KB18","YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58")

##for each strain, filter out samples with only one datapoint (make those all NA), 
##then rank, compute mean rank and percentile rank, plot correlation,
##save table to csv and into list, Spearman's rho to correlations table
single_community_list_ranks <- list()
correlations <- data.frame(strain = character(), comparison = character(), rho = numeric(),p = numeric())


for (i in 1:12){
  
  table <- single_vs_community %>% 
    filter(strain == species_short[[i]]) %>% 
    select(`locus tag_old`,genus, species, strain, `refseq annotation`, koeggNOG,
           `community(batch)_1`, `community(batch)_2`, `community(batch)_3`, starts_with(species_short[[i]]))
  #if there is only one value for a protein, set this to NA also  
  cols <- c(7, 8, 9)  # community set
  table <- table %>% mutate(across(all_of(cols),
        ~ replace(.,rowSums(!is.na(pick(all_of(cols)))) < 2,NA)))
  cols <- c(10, 11, 12)  # monoculture set
  table <- table %>% mutate(across(all_of(cols),
         ~ replace(.,rowSums(!is.na(pick(all_of(cols)))) < 2,NA)))
  rm(cols)
  
  table <- table %>% 
    mutate(rank_community_1 = rank(-.[[7]], na.last = "keep"),
           rank_community_2 = rank(-.[[8]], na.last = "keep"),
           rank_community_3 = rank(-.[[9]], na.last = "keep")) %>%
    mutate(rank_single_1 = rank(-.[[10]], na.last = "keep"),
           rank_single_2 = rank(-.[[11]], na.last = "keep"),
           rank_single_3 = rank(-.[[12]], na.last = "keep")) %>%
    mutate(mean_rank_community = rowMeans(select(., rank_community_1, rank_community_2, rank_community_3), na.rm = TRUE),
           mean_rank_single = rowMeans(select(., rank_single_1, rank_single_2, rank_single_3), na.rm = TRUE))%>%
    mutate(across(starts_with("mean_"), ~replace_na(., NA)))%>%
    mutate(percentile_rank_single = mean_rank_single/max(mean_rank_single, na.rm = TRUE),
          percentile_rank_comm = mean_rank_community/max(mean_rank_community, na.rm = TRUE))
  
  single_community_list_ranks[[species_short[i]]] <- table
  
  #plot correlation
  axis_lim = max(c(table$mean_rank_single, table$mean_rank_community), na.rm = TRUE)*1.1
  p<- ggplot(table, aes(x = mean_rank_single, y = mean_rank_community)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
    scale_y_continuous(limits = c(-1, axis_lim), expand = c(0, 0)) +
    scale_x_continuous(limits = c(-1, axis_lim), expand = c(0, 0)) +
    labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
         x = "Mean Rank in monoculture",
         y = "Mean Rank in community (AF)") +
    theme_minimal()
  
  ggsave(  filename = paste0(paste(species_short[i], "single", "community", sep = "_"), ".png"),
           plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300, 
           path =  "./correlation plots single vs. community")
  ggsave(  filename = paste0(paste(species_short[i], "single", "community", sep = "_"), ".pdf"),
           plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",
           path =  "./correlation plots single vs. community")
  
  #calculate Spearman's rho and save to correlations table
  tmp <- select(table, mean_rank_single, mean_rank_community)
  colnames(tmp) <- c("x", "y")  
  tmp <- na.omit(tmp[, c("x", "y")])
  
  if (nrow(tmp) < 3) {
    rho <- NA
    pval <- NA} else {
    test <- cor.test(tmp$x,tmp$y,method = "spearman",exact = FALSE)
    rho <- unname(test$estimate)
    pval <- test$p.value}
  
  correlations <- rbind(correlations,
                        data.frame(strain = species_short[i],comparison = "single_community", rho = rho, p = pval))
  rm(tmp, rho, pval, test)
  
  #save ranks to csv
  filename <- paste(species_short[[i]], "ranks", "csv", sep = ".")
  write.csv(table, filename, row.names = FALSE)
  rm(table, filename, axis_lim, p)
}


##plot correlations
correlations$strain <- factor(correlations$strain, levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46","I48", "I49", "YL58"))

ggplot(correlations, aes(x = strain, y = rho, fill = strain))+
  geom_col()+
  scale_fill_manual(values = c("#0e8b18","#0000cc","#97b497", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc", "#871900", "#365934","#ff9b05", "#004000","#0b6760"))+
  scale_x_discrete(labels = species)+
  labs(title = "Spearman Correlation Coefficients for Rank Comparison between Monoculture and Community",
       x = "OMM12 species",
       y = "Spearman's rho for correlation between single and community") +
  scale_y_continuous(limits = c(-0.5, 1), expand = c(0, 0)) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.position = "none")

##plot correlations without KB18, YL2, YL44, I49)

species_filtered <- c("E. faecalis KB1",  "M. intestinale YL27", "F. plautii YL31",
             "E. clostridioformis YL32", "T. muris YL45", "C. innocuum I46", 
             "B. caecimuris I48", "B. coccoides YL58")


ggplot(filter(correlations, strain != "KB18" & strain != "YL2" & strain != "YL44" & strain != "I49"), 
       aes(x = strain, y = rho, fill = strain))+
  geom_col()+
  scale_fill_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#871900", "#365934","#ff9b05","#0b6760"))+
  scale_x_discrete(labels = species_filtered)+
  labs(title = "Spearman Correlation Coefficients for Rank Comparison between Monoculture and Community",
       x = "OMM12 species",
       y = "Spearman's rho for correlation between single and community") +
  scale_y_continuous(limits = c(-0.5, 1), expand = c(0, 0)) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.position = "none")
