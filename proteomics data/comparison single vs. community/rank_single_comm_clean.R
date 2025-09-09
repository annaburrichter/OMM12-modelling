library(openxlsx)
library(tidyverse)
library(readxl)
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
library(clusterProfiler)

annotations <- read.delim2("protein_annotation.tsv")
single_vs_community <- read_excel("reanalysed_single_vs_community.xlsx", 
                                             sheet = "single vs community",
                                  col_types = c("skip", "text", "text", "skip", "skip", "text", "text", "skip", 
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))


i48_vs_community <- read_excel("reanalysed_single_vs_community.xlsx", 
                                  sheet = "single vs community_I48", 
                               col_types = c("skip", "text", "text", "skip", "skip", "text", "text", "skip", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))

species <- c("E. faecalis KB1", "B. animalis YL2", "A. muris KB18", "M. intestinale YL27", "F. plautii YL31",
             "E. clostridioformis YL32", "A. muciniphila YL44", "T. muris YL45", "C. innocuum I46", 
             "B. caecimuris I48", "L. reuteri I49", "B. coccoides YL58")
specieswoI48 <- c("E. faecalis KB1", "B. animalis YL2", "A. muris KB18", "M. intestinale YL27", "F. plautii YL31",
             "E. clostridioformis YL32", "A. muciniphila YL44", "T. muris YL45", "C. innocuum I46", 
             "L. reuteri I49", "B. coccoides YL58")
species_short <- c("KB1_", "YL2_", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I49", "YL58")

# test statistics without mean
I48_direct <- i48_vs_community %>% 
  filter(species == "B. caecimuris I48") %>% 
  mutate(rank_I48_1 = rank(-I48_1, na.last = "keep")) %>%
  mutate(rank_I48_2 = rank(-I48_2, na.last = "keep")) %>%
  mutate(rank_I48_3 = rank(-I48_3, na.last = "keep")) %>%
  mutate(rank_comm_1 = rank(-`community(batch)_1`, na.last = "keep")) %>%
  mutate(rank_comm_2 = rank(-`community(batch)_2`, na.last = "keep")) %>%
  mutate(rank_comm_3 = rank(-`community(batch)_3`, na.last = "keep")) %>%
  left_join(annotations, by = c("locus tag_old" = "locus_tag")) %>%
  select(`locus tag_old`, species, `refseq annotation`, koeggNOG, KEGG_ko, rank_I48_1, rank_I48_2, rank_I48_3,rank_comm_1,rank_comm_2, rank_comm_3)
write.csv(I48_direct, "I48_direct.csv")

# repeat for all other strains

for (i in 1:length(specieswoI48)){

table <- single_vs_community %>% 
  filter(species == specieswoI48[[i]]) %>% 
  left_join(annotations, by = c("locus tag_old" = "locus_tag")) %>%
  select(`locus tag_old`, species, `refseq annotation`, koeggNOG, KEGG_ko, `community(batch)_1`, `community(batch)_2`, `community(batch)_3`, starts_with(species_short[[i]]))%>%
  mutate(rank_single_1 = rank(-.[[9]], na.last = "keep")) %>%
  mutate(rank_single_2 = rank(-.[[10]], na.last = "keep")) %>%
  mutate(rank_single_3 = rank(-.[[11]], na.last = "keep")) %>%
  mutate(rank_comm_1 = rank(-`community(batch)_1`, na.last = "keep")) %>%
  mutate(rank_comm_2 = rank(-`community(batch)_2`, na.last = "keep")) %>%
  mutate(rank_comm_3 = rank(-`community(batch)_3`, na.last = "keep"))
table_short <- table %>% select(`locus tag_old`, rank_single_1, rank_single_2, rank_single_3, rank_comm_1, rank_comm_2, rank_comm_3)
filename_short <- paste(species_short[[i]], "direct", "csv", sep = ".")
filename <- paste(species_short[[i]], "ranks", "csv", sep = ".")
write.csv(table_short, filename_short, row.names = FALSE)
write.csv(table, filename)
assign(paste(species_short[[i]], "ranks", sep = "."), table)
rm(table, filename, table_short, filename_short)}


##t-tests in Prism, export discoveries and import into R

sheets <- excel_sheets("discoveries.xlsx")
list_all_up <- data.frame(matrix(nrow=0, ncol=2))
list_all_down <- data.frame(matrix(nrow=0, ncol=2))

# enrich by species (difference cutoff 100)
for (i in 1:12){
  
  disc <- read_excel("discoveries.xlsx", sheet = i, col_names = TRUE) 
  list <- disc %>%
    left_join(annotations, relationship = "many-to-many") %>%
    select(KEGG_ko, Difference) %>%
    filter(!is.na(KEGG_ko)) %>%  filter(KEGG_ko != "-") %>% filter(KEGG_ko != "") %>% 
    separate_longer_delim(KEGG_ko, ",") %>% distinct(KEGG_ko, .keep_all = TRUE)
    list$KEGG_ko <- str_replace_all(list$KEGG_ko, c("ko:" = ""))
  list_comm_up <- filter(list, Difference > 100)
  list_all_up <- bind_rows(list_all_up, list_comm_up)
  list_comm_up <- list_comm_up$KEGG_ko
  if (length(list_comm_up) >0) {
    enrichment_up <- enrichKEGG(list_comm_up, organism = "ko")
  plot_up <- enrichment_up@result %>% 
    separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)
  plot_up$genes <- as.numeric(unlist(plot_up$genes))
  plot_up$all_genes <- as.numeric(unlist(plot_up$all_genes))
  plot_up <- mutate(plot_up, gene_ratio = genes/all_genes)
  plot_up <-  dplyr::filter(plot_up, p.adjust < 0.05)
  plot_up <-  dplyr::slice(plot_up, 1:10)}
  else {
    plot_up <- data.frame(matrix(nrow=0, ncol=17))
  }
  
  list_comm_down <- filter(list, Difference < -100)
  list_all_down <- bind_rows(list_all_down, list_comm_down)
  list_comm_down <- list_comm_down$KEGG_ko
  if (length(list_comm_down) >0) {
  enrichment_down <- enrichKEGG(list_comm_down, organism = "ko")
  plot_down <- enrichment_down@result %>% 
    separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)
  plot_down$genes <- as.numeric(unlist(plot_down$genes))
  plot_down$all_genes <- as.numeric(unlist(plot_down$all_genes))
  plot_down <- mutate(plot_down, gene_ratio = genes/all_genes)
  plot_down$gene_ratio <- plot_down$gene_ratio*(-1)
  plot_down <-  dplyr::filter(plot_down, p.adjust < 0.05)
  plot_down <-  dplyr::slice(plot_down, 1:10)}
  else {
    plot_down <- data.frame(matrix(nrow=0, ncol=17))
  }
  
  plot <- bind_rows(plot_up, plot_down) %>% arrange (desc(gene_ratio))
 
  ggplot(plot, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
      geom_point(shape=16) + 
      scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
      scale_size_continuous(limits = c(0,101)) +
      coord_cartesian(xlim = c(-0.5, 0.5)) + 
      geom_vline(xintercept = 0, color = "grey") + 
      theme_bw(base_size = 12) + 
      labs(x = "gene ratio", y = "", title = paste(sheets[[i]], "monoculture vs. OMM12 community", sep = " " ), caption = c("up in community", "up in monoculture")) + 
      theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
            plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 138, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))
  
  ggsave(filename = paste(sheets[[i]], "enrichment monoculture vs. OMM12 community", "pdf", sep = "." ) , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm", width = 8, height = 6.60, unit = "in")
  ggsave(filename = paste(sheets[[i]], "enrichment monoculture vs. OMM12 community", "png", sep = "." ) , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm", width = 9.86, height = 6.25, unit = "in")
  
  
  assign(sheets[[i]], disc)
  assign(paste(sheets[[i]], "enrichment_up", sep = "_"), enrichment_up)
  assign(paste(sheets[[i]], "enrichment_down", sep = "_"), enrichment_down)
  }

rm(disc, enrichment_down, enrichment_up, list_comm_down, list_comm_up, list, plot_up, plot_down, plot)

#enrich and plot over all species
# filter with difference cutoff 100

list_all_up <-  dplyr::filter(list_all_up, Difference > 100)
list_all_down <-  dplyr::filter(list_all_down, Difference < -100)


list_all_up_ko <- list_all_up$KEGG_ko
enrichment_all_up <- enrichKEGG(list_all_up_ko, organism = "ko")
plot_all_up <- enrichment_all_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)
plot_all_up$genes <- as.numeric(unlist(plot_all_up$genes))
plot_all_up$all_genes <- as.numeric(unlist(plot_all_up$all_genes))
plot_all_up <- mutate(plot_all_up, gene_ratio = genes/all_genes)
plot_all_up <-  dplyr::filter(plot_all_up, p.adjust < 0.05)

list_all_down_ko <- list_all_down$KEGG_ko
enrichment_all_down <- enrichKEGG(list_all_down_ko, organism = "ko")
plot_all_down <- enrichment_all_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)
plot_all_down$genes <- as.numeric(unlist(plot_all_down$genes))
plot_all_down$all_genes <- as.numeric(unlist(plot_all_down$all_genes))
plot_all_down <- mutate(plot_all_down, gene_ratio = genes/all_genes)
plot_all_down$gene_ratio <- plot_all_down$gene_ratio*(-1)
plot_all_down <-  dplyr::filter(plot_all_down, p.adjust < 0.05)

plot_all <- bind_rows(plot_all_up, plot_all_down) %>% arrange (desc(gene_ratio))

ggplot(plot_all, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,101)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "monoculture vs. OMM12 community - all strains", caption = c("up in community", "up in monoculture")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 12), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 230, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment monoculture vs. community - all strains - rank100.pdf" , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm", width = 11.7, height = 10.50, unit = "in")
ggsave(filename = "enrichment monoculture vs. community - all strains - rank100.png" , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm", width = 11.7, height = 9.9, unit = "in")


#export enrichment data to Excel file as supplementary data

# Create a new workbook
wb <- createWorkbook()

# Add worksheets and write data frames to them
addWorksheet(wb, "all species up")
writeData(wb, "all species up", enrichment_all_up@result)
addWorksheet(wb, "all species down")
writeData(wb, "all species down", enrichment_all_down@result)

addWorksheet(wb, "KB1 up")
writeData(wb, "KB1 up", KB1_discoveries_enrichment_up@result)
addWorksheet(wb, "KB1 down")
writeData(wb, "KB1 down", KB1_discoveries_enrichment_down@result)

addWorksheet(wb, "YL2 up")
writeData(wb, "YL2 up", YL2_discoveries_enrichment_up@result)
addWorksheet(wb, "YL2 down")
writeData(wb, "YL2 down", YL2_discoveries_enrichment_down@result)

addWorksheet(wb, "KB18 up")
writeData(wb, "KB18 up", KB18_discoveries_enrichment_up@result)
addWorksheet(wb, "KB18 down")
writeData(wb, "KB18 down", KB18_discoveries_enrichment_down@result)

addWorksheet(wb, "YL27 up")
writeData(wb, "YL27 up", YL27_discoveries_enrichment_up@result)
addWorksheet(wb, "YL27 down")
writeData(wb, "YL27 down", YL27_discoveries_enrichment_down@result)

addWorksheet(wb, "YL31 up")
writeData(wb, "YL31 up", YL31_discoveries_enrichment_up@result)
addWorksheet(wb, "YL31 down")
writeData(wb, "YL31 down", YL31_discoveries_enrichment_down@result)

addWorksheet(wb, "YL32 up")
writeData(wb, "YL32 up", YL32_discoveries_enrichment_up@result)
addWorksheet(wb, "YL32 down")
writeData(wb, "YL32 down", YL32_discoveries_enrichment_down@result)

addWorksheet(wb, "YL44 up")
writeData(wb, "YL44 up", YL44_discoveries_enrichment_up@result)
addWorksheet(wb, "YL44 down")
writeData(wb, "YL44 down", YL44_discoveries_enrichment_down@result)

addWorksheet(wb, "YL45 up")
writeData(wb, "YL45 up", YL45_discoveries_enrichment_up@result)
addWorksheet(wb, "YL45 down")
writeData(wb, "YL45 down", YL45_discoveries_enrichment_down@result)

addWorksheet(wb, "I46 up")
writeData(wb, "I46 up", I46_discoveries_enrichment_up@result)
addWorksheet(wb, "I46 down")
writeData(wb, "I46 down", I46_discoveries_enrichment_down@result)

addWorksheet(wb, "I48 up")
writeData(wb, "I48 up", I48_discoveries_enrichment_up@result)
addWorksheet(wb, "I48 down")
writeData(wb, "I48 down", I48_discoveries_enrichment_down@result)

addWorksheet(wb, "I49 up")
writeData(wb, "I49 up", I49_discoveries_enrichment_up@result)
addWorksheet(wb, "I49 down")
writeData(wb, "I49 down", I49_discoveries_enrichment_down@result)

addWorksheet(wb, "YL58 up")
writeData(wb, "YL58 up", YL58_discoveries_enrichment_up@result)
addWorksheet(wb, "YL58 down")
writeData(wb, "YL58 down", YL58_discoveries_enrichment_down@result)

# Save the workbook to an Excel file
saveWorkbook(wb, "Table S3 full datasets for enrichment monoculture - community.xlsx", overwrite = TRUE)


#volcano plot of all single vs community - cutoff 100

library(prism2R)

prism_file <- get_analysis_tables("single vs community rank t-tests.prism")
ann_short <- select(vitro_vs_vivo, `locus tag_old`, species, `refseq annotation`, koeggNOG)
ann_short <- rename(ann_short,locus_tag = `locus tag_old`)
colors <- c("#ff9b05","#365934","#004000","#0000cc","#cc6a00","#7bc57b","#11a196","#6b00bc","#871900","#0b6760", "#0e8b18", "#97b497")

for (i in 1:12){

prism <- prism_file[[i]]
prism <- rename(prism,locus_tag = X1)
prism$`q value`<- as.numeric(prism$`q value`)
prism <- prism %>% full_join(ann_short, relationship = "many-to-many")
prism <- mutate(prism, log_q = -log(`q value`))
prism$colorhelp <- ifelse(prism$log_q > 1.301029995664e0, "sig", "ns")
prism$color <- ifelse(prism$Difference < -100|prism$Difference > 100, prism$colorhelp, "ns")


ggplot(prism, aes(x=Difference, y=log_q, color = color)) + 
  geom_point(shape=16) + 
  scale_color_manual(values=c("lightgrey", colors[i]))+
  geom_vline(xintercept = 0, color = "black") + 
  geom_vline(xintercept = 100, color = "grey") + 
  geom_vline(xintercept = -100, color = "grey") + 
  geom_hline(yintercept = 1.301029995664e0, color = "grey") +
  coord_cartesian(xlim = c(-1200, 1200),ylim = c(0, 15)) + 
  theme_classic(base_size=14) + 
  labs(x = "Rank Difference monoculture vs. community", y = "-log(q value)", title = paste("single vs. community", names(prism_file)[i], sep = ""), caption = c("up in community", "up in monoculture")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 5), axis.text = element_text(size = 12), plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 130, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"), legend.position = "none")
  

ggsave(filename = paste("single vs. community", "rank100", names(prism_file)[i], "pdf", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")
ggsave(filename = paste("single vs. community", "rank100",names(prism_file)[i], "png", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")


assign(names(prism_file)[i], p)
}

mark_YL27 <- c("I5Q90_01190", "I5Q90_07635", "I5Q90_05770", "I5Q90_05760")
mark_YL31 <- c("I5Q84_01935", "I5Q84_01940", "I5Q84_01950", "I5Q84_01955", "I5Q84_01970", "I5Q84_01980")
mark_YL32 <- c("I5Q83_03300", "I5Q83_03305", "I5Q83_03225", "I5Q83_03315", "I5Q83_03295",
               "I5Q83_03230", "I5Q83_03310", "I5Q83_03240", "I5Q83_27725", "I5Q83_19950")
mark_I46 <- c("I5Q87_12730", "I5Q87_16370", "I5Q87_09085", "I5Q87_09075","I5Q87_12070")
mark_I48 <- c("I5Q79_08170", "I5Q79_10475", "I5Q79_08180", "I5Q79_08185", "I5Q79_08165", "I5Q79_10470", 
              "I5Q79_07450", "I5Q79_07445", "I5Q79_08670", "I5Q79_12300", "I5Q79_00510", "I5Q79_12305", 
              "I5Q79_07420", "I5Q79_08665", "I5Q79_10485")
mark_YL58 <- c("I5Q86_17800", "I5Q86_15835", "I5Q86_11945", "I5Q86_20645", "I5Q86_03035", "I5Q86_06130", 
               "I5Q86_21630", "I5Q86_21635", "I5Q86_21640")

# marked plot YL27
i <- 5
prism <- prism_file[[5]]
prism <- rename(prism,locus_tag = X1)
prism$`q value`<- as.numeric(prism$`q value`)
prism <- prism %>% full_join(ann_short, relationship = "many-to-many")
prism <- mutate(prism, log_q = -log(`q value`))
prism$colorhelp1 <- ifelse(prism$log_q > 1.301029995664e0, "sig", "ns")
prism$colorhelp2 <- ifelse(prism$Difference < -100|prism$Difference > 100, prism$colorhelp1, "ns")
prism$color<- ifelse(prism$locus_tag %in% mark_YL27, "special", prism$colorhelp2)


ggplot(prism, aes(x=Difference, y=log_q, color = color)) + 
  geom_point(shape=16) + 
  scale_color_manual(values=c("lightgrey", colors[i], "coral2"))+
  geom_vline(xintercept = 0, color = "black") + 
  geom_vline(xintercept = 100, color = "grey") + 
  geom_vline(xintercept = -100, color = "grey") + 
  geom_hline(yintercept = 1.301029995664e0, color = "grey") +
  coord_cartesian(xlim = c(-1200, 1200),ylim = c(0, 15)) + 
  theme_classic(base_size=14) + 
  labs(x = "Rank Difference monoculture vs. community", y = "-log(q value)", title = paste("single vs. community", names(prism_file)[i], sep = ""), caption = c("up in community", "up in monoculture")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 5), axis.text = element_text(size = 12), plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 130, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"), legend.position = "none")


ggsave(filename = paste("single vs. community","rank100", names(prism_file)[i], "marked", "pdf", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")
ggsave(filename = paste("single vs. community","rank100", names(prism_file)[i], "marked", "png", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")


# marked plot YL31
i <- 6
prism <- prism_file[[6]]
prism <- rename(prism,locus_tag = X1)
prism$`q value`<- as.numeric(prism$`q value`)
prism <- prism %>% full_join(ann_short, relationship = "many-to-many")
prism <- mutate(prism, log_q = -log(`q value`))
prism$colorhelp1 <- ifelse(prism$log_q > 1.301029995664e0, "sig", "ns")
prism$colorhelp2 <- ifelse(prism$Difference < -100|prism$Difference > 100, prism$colorhelp1, "ns")
prism$color<- ifelse(prism$locus_tag %in% mark_YL31, "special", prism$colorhelp2)

ggplot(prism, aes(x=Difference, y=log_q, color = color)) + 
  geom_point(shape=16) + 
  scale_color_manual(values=c("lightgrey", colors[i], "coral2"))+
  geom_vline(xintercept = 0, color = "black") + 
  geom_vline(xintercept = 100, color = "grey") + 
  geom_vline(xintercept = -100, color = "grey") + 
  geom_hline(yintercept = 1.301029995664e0, color = "grey") +
  coord_cartesian(xlim = c(-1200, 1200),ylim = c(0, 15)) + 
  theme_classic(base_size=14) + 
  labs(x = "Rank Difference monoculture vs. community", y = "-log(q value)", title = paste("single vs. community", names(prism_file)[i], sep = ""), caption = c("up in community", "up in monoculture")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 5), axis.text = element_text(size = 12), plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 130, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"), legend.position = "none")


ggsave(filename = paste("single vs. community","rank100", names(prism_file)[i], "marked", "pdf", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")
ggsave(filename = paste("single vs. community", "rank100",names(prism_file)[i], "marked", "png", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")


# marked plot YL32
i<- 7
prism <- prism_file[[7]]
prism <- rename(prism,locus_tag = X1)
prism$`q value`<- as.numeric(prism$`q value`)
prism <- prism %>% full_join(ann_short, relationship = "many-to-many")
prism <- mutate(prism, log_q = -log(`q value`))
prism$colorhelp1 <- ifelse(prism$log_q > 1.301029995664e0, "sig", "ns")
prism$colorhelp2 <- ifelse(prism$Difference < -100|prism$Difference > 100, prism$colorhelp1, "ns")
prism$color<- ifelse(prism$locus_tag %in% mark_YL32, "special", prism$colorhelp2)


ggplot(prism, aes(x=Difference, y=log_q, color = color)) + 
  geom_point(shape=16) + 
  scale_color_manual(values=c("lightgrey", colors[i], "coral2"))+
  geom_vline(xintercept = 0, color = "black") + 
  geom_vline(xintercept = 100, color = "grey") + 
  geom_vline(xintercept = -100, color = "grey") + 
  geom_hline(yintercept = 1.301029995664e0, color = "grey") +
  coord_cartesian(xlim = c(-1200, 1200),ylim = c(0, 15)) + 
  theme_classic(base_size=14) + 
  labs(x = "Rank Difference monoculture vs. community", y = "-log(q value)", title = paste("single vs. community", names(prism_file)[i], sep = ""), caption = c("up in community", "up in monoculture")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 5), axis.text = element_text(size = 12), plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 130, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"), legend.position = "none")


ggsave(filename = paste("single vs. community", "rank100",names(prism_file)[i], "marked", "pdf", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")
ggsave(filename = paste("single vs. community", "rank100",names(prism_file)[i], "marked", "png", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")


# marked plot I48
i <- 1
prism <- prism_file[[1]]
prism <- rename(prism,locus_tag = X1)
prism$`q value`<- as.numeric(prism$`q value`)
prism <- prism %>% full_join(ann_short, relationship = "many-to-many")
prism <- mutate(prism, log_q = -log(`q value`))
prism$colorhelp1 <- ifelse(prism$log_q > 1.301029995664e0, "sig", "ns")
prism$colorhelp2 <- ifelse(prism$Difference < -100|prism$Difference > 100, prism$colorhelp1, "ns")
prism$color<- ifelse(prism$locus_tag %in% mark_I48, "special", prism$colorhelp2)


ggplot(prism, aes(x=Difference, y=log_q, color = color)) + 
  geom_point(shape=16) + 
  scale_color_manual(values=c("lightgrey", colors[i], "coral2"))+
  geom_vline(xintercept = 0, color = "black") + 
  geom_vline(xintercept = 100, color = "grey") + 
  geom_vline(xintercept = -100, color = "grey") + 
  geom_hline(yintercept = 1.301029995664e0, color = "grey") +
  coord_cartesian(xlim = c(-1200, 1200),ylim = c(0, 15)) + 
  theme_classic(base_size=14) + 
  labs(x = "Rank Difference monoculture vs. community", y = "-log(q value)", title = paste("single vs. community", names(prism_file)[i], sep = ""), caption = c("up in community", "up in monoculture")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 5), axis.text = element_text(size = 12), plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 130, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"), legend.position = "none")


ggsave(filename = paste("single vs. community", "rank100",names(prism_file)[i], "marked", "pdf", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")
ggsave(filename = paste("single vs. community", "rank100",names(prism_file)[i], "marked", "png", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")


# marked plot I46
i<- 2
prism <- prism_file[[2]]
prism <- rename(prism,locus_tag = X1)
prism$`q value`<- as.numeric(prism$`q value`)
prism <- prism %>% full_join(ann_short, relationship = "many-to-many")
prism <- mutate(prism, log_q = -log(`q value`))
prism$colorhelp1 <- ifelse(prism$log_q > 1.301029995664e0, "sig", "ns")
prism$colorhelp2 <- ifelse(prism$Difference < -100|prism$Difference > 100, prism$colorhelp1, "ns")
prism$color<- ifelse(prism$locus_tag %in% mark_I46, "special", prism$colorhelp2)


ggplot(prism, aes(x=Difference, y=log_q, color = color)) + 
  geom_point(shape=16) + 
  scale_color_manual(values=c("lightgrey", colors[i], "coral2"))+
  geom_vline(xintercept = 0, color = "black") + 
  geom_vline(xintercept = 100, color = "grey") + 
  geom_vline(xintercept = -100, color = "grey") + 
  geom_hline(yintercept = 1.301029995664e0, color = "grey") +
  coord_cartesian(xlim = c(-1200, 1200),ylim = c(0, 15)) + 
  theme_classic(base_size=14) + 
  labs(x = "Rank Difference monoculture vs. community", y = "-log(q value)", title = paste("single vs. community", names(prism_file)[i], sep = ""), caption = c("up in community", "up in monoculture")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 5), axis.text = element_text(size = 12), plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 130, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"), legend.position = "none")


ggsave(filename = paste("single vs. community", "rank100",names(prism_file)[i], "marked", "pdf", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")
ggsave(filename = paste("single vs. community", "rank100",names(prism_file)[i], "marked", "png", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")


# marked plot YL58
i<- 10
prism <- prism_file[[10]]
prism <- rename(prism,locus_tag = X1)
prism$`q value`<- as.numeric(prism$`q value`)
prism <- prism %>% full_join(ann_short, relationship = "many-to-many")
prism <- mutate(prism, log_q = -log(`q value`))
prism$colorhelp1 <- ifelse(prism$log_q > 1.301029995664e0, "sig", "ns")
prism$colorhelp2 <- ifelse(prism$Difference < -100|prism$Difference > 100, prism$colorhelp1, "ns")
prism$color<- ifelse(prism$locus_tag %in% mark_YL58, "special", prism$colorhelp2)


ggplot(prism, aes(x=Difference, y=log_q, color = color)) + 
  geom_point(shape=16) + 
  scale_color_manual(values=c("lightgrey", colors[i], "coral2"))+
  geom_vline(xintercept = 0, color = "black") + 
  geom_vline(xintercept = 100, color = "grey") + 
  geom_vline(xintercept = -100, color = "grey") + 
  geom_hline(yintercept = 1.301029995664e0, color = "grey") +
  coord_cartesian(xlim = c(-1200, 1200),ylim = c(0, 15)) + 
  theme_classic(base_size=14) + 
  labs(x = "Rank Difference monoculture vs. community", y = "-log(q value)", title = paste("single vs. community", names(prism_file)[i], sep = ""), caption = c("up in community", "up in monoculture")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 5), axis.text = element_text(size = 12), plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 130, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"), legend.position = "none")


ggsave(filename = paste("single vs. community", "rank100",names(prism_file)[i], "marked", "pdf", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")
ggsave(filename = paste("single vs. community", "rank100",names(prism_file)[i], "marked", "png", sep = "." ), path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")


## count proteins with significant rank change 

all_discoveries <- data.frame(matrix(nrow=0, ncol=0))

for (i in 1:12){
  
  disc <- read_excel("discoveries.xlsx", sheet = i, col_names = TRUE, 
                     col_types = c("text", "text", "numeric", "numeric",
                                   "numeric", "numeric","numeric", "numeric", "text"))
  disc <- left_join(disc, annotations, relationship = "many-to-many")
  all_discoveries <- bind_rows(all_discoveries, disc)
  rm(disc)}

all_discoveries_filtered_50 <- all_discoveries %>%
                                dplyr::filter(Difference > 50|Difference < -50) %>%
                                dplyr::filter(!is.na(org)) %>%
                                dplyr::group_by(org)
cutoff_50 <- tally(all_discoveries_filtered_50)
cutoff_50$org <- factor(cutoff_50$org, levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58"))


all_discoveries_filtered_100 <- all_discoveries %>%
  dplyr::filter(Difference > 100|Difference < -100) %>%
  dplyr::filter(!is.na(org)) %>%
  dplyr::group_by(org)
cutoff_100 <- tally(all_discoveries_filtered_100)
cutoff_100$org <- factor(cutoff_100$org, levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58"))
species_it <- expression(italic("E. faecalis")~"KB1", italic("B.animalis")~"YL2",
                         italic("A. muris")~"KB18", italic("M. intestinale")~"YL27",
                         italic("F. plautii")~"YL31", italic("E. clostridioformis")~"YL32",
                         italic("A. muciniphila")~"YL44", italic("T. muris")~"YL45",
                         italic("C. innocuum")~"I46", italic("B.caecimuris")~"I48",
                         italic("L. reuteri")~"I49", italic("B.coccoides")~"YL58") 

colors <-  c("#0e8b18","#0000cc","#97b497", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc", "#871900", "#365934","#ff9b05", "#004000","#0b6760")

ggplot(cutoff_50, aes(x=n, y = factor(org, levels = c("YL58", "I49", "I48", "I46", "YL45", "YL44","YL32" , "YL31", "YL27", "KB18","YL2" ,"KB1")), fill = org)) + 
  geom_col() + 
  scale_fill_manual(values = colors) + 
  scale_y_discrete(labels = rev(species_it)) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = 'none', axis.text = element_text(size = 14)) + 
  labs(y="", x = "number of proteins", title = "Number of proteins with significant rank change > 50")

ggsave(filename = "rank change 50.pdf" , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")
ggsave(filename = "rank change 50.png" , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")


ggplot(cutoff_100, aes(x=n, y = factor(org, levels = c("YL58", "I49", "I48", "I46", "YL45", "YL44","YL32" , "YL31", "YL27", "KB18","YL2" ,"KB1")), fill = org)) + 
  geom_col() + 
  scale_fill_manual(values = colors) + 
  scale_y_discrete(labels = rev(species_it)) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = 'none', axis.text = element_text(size = 14)) + 
  labs(y="", x = "number of proteins", title = "Number of proteins with significant rank change > 100")

ggsave(filename = "rank change 100.pdf" , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")
ggsave(filename = "rank change 100.png" , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")


#cutoff 100, two-sided

all_disc_100_up <- dplyr::filter(all_discoveries_filtered_100, Difference >100) %>% dplyr::group_by(org)
cutoff_100_up <- tally(all_disc_100_up)
all_disc_100_down <- dplyr::filter(all_discoveries_filtered_100, Difference < (-100)) %>% dplyr::group_by(org)
cutoff_100_down <- tally(all_disc_100_down) %>% mutate(n = (-1)*n) %>% rename(m = n)
cutoff_100_two_sided <- left_join(cutoff_100_up, cutoff_100_down)%>%pivot_longer(-org)


ggplot(cutoff_100_two_sided, aes(x = factor(org, levels = c("YL58", "I49", "I48", "I46", "YL45", "YL44","YL32" , "YL31", "YL27", "KB18","YL2" ,"KB1")), y = value, 
                                 fill = factor(org, levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32","YL44" , "YL45", "I46", "I48","I49" ,"YL58")))) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip(ylim = c(-200, 200)) + # Flip coordinates for horizontal bars
  scale_y_continuous(labels = abs) + # Show positive labels on both sides
  geom_hline(yintercept = 0, color = "black")+ 
  labs(
    title = "Number of proteins with significant rank change > 100",
    x = "species",
    y = "number of proteins"
  ) +
  theme_bw(base_size =14) + theme(legend.title=element_blank())+
  scale_fill_manual(values = colors)

ggsave(filename = "rank change 100.pdf" , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")
ggsave(filename = "rank change 100.png" , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")

##differently regulated proteins as percentage of all proteins detected in community

prism <- prism_file[[1]]
proteins_I48 <- sum(!is.na(prism[,5]))
prism <- prism_file[[11]]
proteins_KB1 <- sum(!is.na(prism[,5]))
prism <- prism_file[[5]]
proteins_YL27 <- sum(!is.na(prism[,5]))
prism <- prism_file[[6]]
proteins_YL31 <- sum(!is.na(prism[,5]))
prism <- prism_file[[7]]
proteins_YL32 <- sum(!is.na(prism[,5]))
prism <- prism_file[[9]]
proteins_YL45 <- sum(!is.na(prism[,5]))
prism <- prism_file[[10]]
proteins_YL58 <- sum(!is.na(prism[,5]))

protein_numbers <- data.frame(org = c("I48", "KB1", "YL27", "YL31", "YL32", "YL45", "YL58"), 
                              proteins = c(proteins_I48, proteins_KB1, proteins_YL27, proteins_YL31, proteins_YL32, proteins_YL45, proteins_YL58))

perc_down <- left_join(cutoff_100_down, protein_numbers)%>%mutate(down=m/proteins)
perc_up <- left_join(cutoff_100_up, protein_numbers)%>%mutate(up=n/proteins)

perc_two_sided <- left_join(perc_up, perc_down)%>% select(org, up, down)%>%pivot_longer(-org)


ggplot(perc_two_sided, aes(x = factor(org, levels = c("YL58", "I49", "I48", "I46", "YL45", "YL44","YL32" , "YL31", "YL27", "KB18","YL2" ,"KB1")), y = value, 
                                 fill = factor(org, levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32","YL44" , "YL45", "I46", "I48","I49" ,"YL58")))) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip(ylim = c(-0.5, 0.5)) + # Flip coordinates for horizontal bars
  scale_y_continuous(labels = abs) + # Show positive labels on both sides
  geom_hline(yintercept = 0, color = "black")+ 
  labs(
    title = "Fraction of proteins with significant rank change > 100",
    x = "species",
    y = "fraction of regulated proteins"
  ) +
  theme_bw(base_size =14) + theme(legend.title=element_blank())+
  scale_fill_manual(values = colors)

ggsave(filename = "fraction rank change 100.pdf" , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")
ggsave(filename = "fraction rank change 100.png" , path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/plots_single_comm")

