library(openxlsx)
library(tidyverse)
library(readxl)
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
library(clusterProfiler)

annotations <- read.delim2("protein_annotation.tsv")



# enrich and plot by species (rank cutoff 50) -----------------------------

#create list to collect data over all species
list_all_up <- tibble()
list_all_down <- tibble()
all_regulated_static <- list()
enrichments_list_static <- list()

for (i in 1:12){
  
  data <- single_community_list_ranks[[i]]%>%
    mutate(rank_difference = mean_rank_single-mean_rank_community)%>%
    rename_with(.cols = 10, ~"single_1")%>%
    rename_with(.cols = 11, ~"single_2")%>%
    rename_with(.cols = 12, ~"single_3")%>%
    left_join(annotations,by = c("locus tag_old" = "locus_tag"), relationship = "many-to-many")
  sig_up <- data %>% filter(rank_difference >50)
  exclusive_up <- data %>% filter(is.na(rank_difference) & percentile_rank_comm <0.1)
  all_regulated_up <- bind_rows(sig_up, exclusive_up)
  all_regulated_static[[paste(species_short[i], "_up")]] <- all_regulated_up
    list_up <- bind_rows(sig_up, exclusive_up)%>%
    select(KEGG_ko) %>%
    filter(!is.na(KEGG_ko)) %>%  filter(KEGG_ko != "-") %>% filter(KEGG_ko != "") %>% 
    separate_longer_delim(KEGG_ko, ",") %>% distinct(KEGG_ko, .keep_all = TRUE)%>%
    mutate(KEGG_ko = str_replace_all(KEGG_ko, "ko:", ""))
  list_all_up <- bind_rows(list_all_up, list_up)
  list_up <- list_up$KEGG_ko
  
  if (length(list_up) >0) {
    enrichment_up <- enrichKEGG(list_up, organism = "ko")
    plot_up <- enrichment_up@result %>% 
      separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)
    plot_up$genes <- as.numeric(unlist(plot_up$genes))
    plot_up$all_genes <- as.numeric(unlist(plot_up$all_genes))
    plot_up <- mutate(plot_up, gene_ratio = genes/all_genes)%>% arrange(desc(gene_ratio))
    plot_up <-  dplyr::filter(plot_up, p.adjust < 0.05)
    plot_up <-  dplyr::slice(plot_up, 1:10)} else {
    plot_up <- data.frame(matrix(nrow=0, ncol=17))
  }
  
  sig_down <- data %>% filter(rank_difference < -50)
  exclusive_down <- data %>% filter(is.na(rank_difference) & percentile_rank_single <0.1)
  all_regulated_down <- bind_rows(sig_down, exclusive_down)
  all_regulated_static[[paste0(species_short[i], "_down")]] <- all_regulated_down
  list_down <- bind_rows(sig_down, exclusive_down)%>%
    select(KEGG_ko) %>%
    filter(!is.na(KEGG_ko)) %>%  filter(KEGG_ko != "-") %>% filter(KEGG_ko != "") %>% 
    separate_longer_delim(KEGG_ko, ",") %>% distinct(KEGG_ko, .keep_all = TRUE)%>%
    mutate(KEGG_ko = str_replace_all(KEGG_ko, "ko:", ""))
  list_all_down <- bind_rows(list_all_down, list_down)
  list_down <- list_down$KEGG_ko
  if (length(list_down) >0) {
    enrichment_down <- enrichKEGG(list_down, organism = "ko")
    plot_down <- enrichment_down@result %>% 
      separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)
    plot_down$genes <- as.numeric(unlist(plot_down$genes))
    plot_down$all_genes <- as.numeric(unlist(plot_down$all_genes))
    plot_down <- mutate(plot_down, gene_ratio = genes/all_genes)%>% arrange(desc(gene_ratio))
    plot_down$gene_ratio <- plot_down$gene_ratio*(-1)
    plot_down <-  dplyr::filter(plot_down, p.adjust < 0.05)
    plot_down <-  dplyr::slice(plot_down, 1:10)} else {
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
    labs(x = "gene ratio", y = "", title = paste(species[[i]], "monoculture vs. OMM12 community", sep = " " ), caption = c("up in community", "up in monoculture")) + 
    theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
          plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 138, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))
  
  ggsave(filename = paste(species[[i]], "enrichment monoculture vs. OMM12 community (static cutoff 50)", "pdf", sep = "." ) , path = "./enrichment_plots_static", width = 8, height = 6.60, unit = "in")
  ggsave(filename = paste(species[[i]], "enrichment monoculture vs. OMM12 community (static cutoff 50)", "png", sep = "." ) , path = "./enrichment_plots_static", width = 9.86, height = 6.25, unit = "in")

  enrichments_list_static[[paste0(species_short[i], "enrichment_down")]] <- enrichment_down
  enrichments_list_static[[paste0(species_short[i], "enrichment_up")]] <- enrichment_up
  rm(data,all_regulated_down, all_regulated_up,enrichment_down, enrichment_up, list_down, list_up, sig_up, sig_down, exclusive_up, exclusive_down, plot_up, plot_down, plot)

}



# enrich and plot by species (rank cutoff 50), plotting all categories--------

for (i in 1:12){
  
  data <- single_community_list_ranks[[i]]%>%
    mutate(rank_difference = mean_rank_single-mean_rank_community)%>%
    rename_with(.cols = 10, ~"single_1")%>%
    rename_with(.cols = 11, ~"single_2")%>%
    rename_with(.cols = 12, ~"single_3")%>%
    left_join(annotations,by = c("locus tag_old" = "locus_tag"), relationship = "many-to-many")
  sig_up <- data %>% filter(rank_difference >50)
  exclusive_up <- data %>% filter(is.na(rank_difference) & percentile_rank_comm <0.1)
  all_regulated_up <- bind_rows(sig_up, exclusive_up)
  list_up <- bind_rows(sig_up, exclusive_up)%>%
    select(KEGG_ko) %>%
    filter(!is.na(KEGG_ko)) %>%  filter(KEGG_ko != "-") %>% filter(KEGG_ko != "") %>% 
    separate_longer_delim(KEGG_ko, ",") %>% distinct(KEGG_ko, .keep_all = TRUE)%>%
    mutate(KEGG_ko = str_replace_all(KEGG_ko, "ko:", ""))
  list_up <- list_up$KEGG_ko
  
  if (length(list_up) >0) {
    enrichment_up <- enrichKEGG(list_up, organism = "ko")
    plot_up <- enrichment_up@result %>% 
      separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)
    plot_up$genes <- as.numeric(unlist(plot_up$genes))
    plot_up$all_genes <- as.numeric(unlist(plot_up$all_genes))
    plot_up <- mutate(plot_up, gene_ratio = genes/all_genes)%>% arrange(desc(gene_ratio))
    plot_up <-  dplyr::filter(plot_up, p.adjust < 0.05)} else {
    plot_up <- data.frame(matrix(nrow=0, ncol=17))
    }
  
  sig_down <- data %>% filter(rank_difference < -50)
  exclusive_down <- data %>% filter(is.na(rank_difference) & percentile_rank_single <0.1)
  all_regulated_down <- bind_rows(sig_down, exclusive_down)
  list_down <- bind_rows(sig_down, exclusive_down)%>%
    select(KEGG_ko) %>%
    filter(!is.na(KEGG_ko)) %>%  filter(KEGG_ko != "-") %>% filter(KEGG_ko != "") %>% 
    separate_longer_delim(KEGG_ko, ",") %>% distinct(KEGG_ko, .keep_all = TRUE)%>%
    mutate(KEGG_ko = str_replace_all(KEGG_ko, "ko:", ""))
  list_down <- list_down$KEGG_ko
  if (length(list_down) >0) {
    enrichment_down <- enrichKEGG(list_down, organism = "ko")
    plot_down <- enrichment_down@result %>% 
      separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)
    plot_down$genes <- as.numeric(unlist(plot_down$genes))
    plot_down$all_genes <- as.numeric(unlist(plot_down$all_genes))
    plot_down <- mutate(plot_down, gene_ratio = genes/all_genes)%>% arrange(desc(gene_ratio))
    plot_down$gene_ratio <- plot_down$gene_ratio*(-1)
    plot_down <-  dplyr::filter(plot_down, p.adjust < 0.05)} else {
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
    labs(x = "gene ratio", y = "", title = paste(species[[i]], "monoculture vs. OMM12 community", sep = " " ), caption = c("up in community", "up in monoculture")) + 
    theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
          plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 138, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))
  
  plot_length <- n_distinct(plot$Description)
  ggsave(filename = paste(species[[i]], "enrichment monoculture vs. OMM12 community (static cutoff 50) long", "pdf", sep = "." ) , path = "./enrichment_plots_static_long", width = 8, height = 0.6 + 0.375 * plot_length, unit = "in")
  ggsave(filename = paste(species[[i]], "enrichment monoculture vs. OMM12 community (static cutoff 50) long", "png", sep = "." ) , path = "./enrichment_plots_static_long", width = 9.86, height = 6.25, unit = "in")
  
  rm(data,all_regulated_down, all_regulated_up,enrichment_down, enrichment_up, list_down, list_up, sig_up, sig_down, exclusive_up, exclusive_down, plot_up, plot_down, plot)
  
}

# enrich and plot over all species ----------------------------------------

list_all_up<- list_all_up%>% distinct(KEGG_ko, .keep_all = TRUE)
list_all_up_ko <- list_all_up$KEGG_ko

enrichment_all_up <- enrichKEGG(list_all_up_ko, organism = "ko")
plot_all_up <- enrichment_all_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)
plot_all_up$genes <- as.numeric(unlist(plot_all_up$genes))
plot_all_up$all_genes <- as.numeric(unlist(plot_all_up$all_genes))
plot_all_up <- mutate(plot_all_up, gene_ratio = genes/all_genes)
plot_all_up <-  dplyr::filter(plot_all_up, p.adjust < 0.05)

list_all_down<- list_all_down%>% distinct(KEGG_ko, .keep_all = TRUE)
list_all_down_ko <- list_all_down$KEGG_ko

enrichment_all_down <- enrichKEGG(list_all_down_ko, organism = "ko")
plot_all_down <- enrichment_all_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)
plot_all_down$genes <- as.numeric(unlist(plot_all_down$genes))
plot_all_down$all_genes <- as.numeric(unlist(plot_all_down$all_genes))
plot_all_down <- mutate(plot_all_down, gene_ratio = genes/all_genes)
plot_all_down$gene_ratio <- plot_all_down$gene_ratio*(-1)
plot_all_down <-  dplyr::filter(plot_all_down, p.adjust < 0.05)

enrichments_list_static[["all_enrichment_down"]] <- enrichment_all_down
enrichments_list_static[["all_enrichment_up"]] <- enrichment_all_up

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

plot_length <- n_distinct(plot_all$Description)
ggsave(filename = "enrichment monoculture vs. community - all strains (static cutoff 50).pdf" , 
       path =  "./enrichment_plots_static", width = 11.7, height = 0.6 + 0.375 * plot_length, unit = "in")
ggsave(filename = "enrichment monoculture vs. community - all strains (static cutoff 50).png" , 
       path =  "./enrichment_plots_static", width = 11.7, height = 9.9, unit = "in")

rm(enrichment_all_down, enrichment_all_up, list_all_down_ko, list_all_up_ko, list_all_down, list_all_up, plot_all_down, plot_all_up, plot_all)

#export enrichments results to Excel

wb <- createWorkbook()

for (i in seq_along(enrichments_list_static)) {
  sheet_name <- names(enrichments_list_static)[i]
  result_df <- enrichments_list_static[[i]]@result
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = result_df)
}

saveWorkbook(wb, "enrichment_results.xlsx", overwrite = TRUE)

# count proteins with significant rank change+exclusive -------------------
#remove KB18 columns from KB1 table and YL27 columns from YL2 table

all_regulated_static[[1]] <- all_regulated_static[[1]]%>%select(-KB18_1, -KB18_2, -KB18_3)
all_regulated_static[[2]] <- all_regulated_static[[2]]%>%select(-KB18_1, -KB18_2, -KB18_3)
all_regulated_static[[3]] <- all_regulated_static[[3]]%>%select(-YL27_1, -YL27_2, -YL27_3)
all_regulated_static[[4]] <- all_regulated_static[[4]]%>%select(-YL27_1, -YL27_2, -YL27_3)

all_discoveries_up <- all_regulated_static[grepl("_up$", names(all_regulated_static))]%>%
  bind_rows()
all_discoveries_down <- all_regulated_static[grepl("_down$", names(all_regulated_static))]%>%
  bind_rows()

#cutoff 50, two-sided

cutoff_50_up <- all_discoveries_up %>%dplyr::group_by(org)%>%tally()
cutoff_50_down <- all_discoveries_down %>%dplyr::group_by(org)%>%tally()%>% mutate(n = (-1)*n) %>% rename(m = n)

cutoff_50_two_sided <- full_join(cutoff_50_up, cutoff_50_down)%>%pivot_longer(-org)%>%
  filter(!is.na(org))

colors <-  c("#0e8b18","#0000cc","#97b497", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc", "#871900", "#365934","#ff9b05", "#004000","#0b6760")

ggplot(cutoff_50_two_sided, aes(x = factor(org, levels = c("YL58", "I49", "I48", "I46", "YL45", "YL44","YL32" , "YL31", "YL27", "KB18","YL2" ,"KB1")), y = value, 
                                fill = factor(org, levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32","YL44" , "YL45", "I46", "I48","I49" ,"YL58")))) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip(ylim = c(-400, 400)) + # Flip coordinates for horizontal bars
  scale_y_continuous(labels = abs) + # Show positive labels on both sides
  geom_hline(yintercept = 0, color = "black")+ 
  labs(
    title = "Number of proteins with rank change > 50 + exclusive proteins",
    x = "species",
    y = "number of regulated proteins",
    caption = c("up in community", "up in monoculture")
  ) +
  theme_bw(base_size =14) + theme(legend.title=element_blank(), 
                                  plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 115, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"),
                                  plot.title = element_text(vjust = 10, size = 14))+
  scale_fill_manual(values = colors)

ggsave(filename = "rank change 50 + exclusive.pdf" , path =   "./enrichment_plots_static",
       width = 8, height = 6, units = "in")
ggsave(filename = "rank change 50 + exclusive.png" , path =   "./enrichment_plots_static",
       width = 8, height = 6, units = "in")



# differently regulated proteins as percentage of all proteins det --------

protein_numbers <- purrr::imap_dfr(single_community_list_ranks,~ tibble(org = .y, n_non_na = sum(!is.na(.x$mean_rank_community)|!is.na(.x$mean_rank_single))))%>%
  mutate(org = str_replace_all(org, ".ranks", ""))

perc_two_sided <- cutoff_50_two_sided%>%
  pivot_wider(names_from = name, values_from = value)%>%
  left_join(protein_numbers)%>%
  mutate(down=m/n_non_na, up=n/n_non_na)%>%
  select(org, up, down)%>%pivot_longer(-org)


ggplot(perc_two_sided, aes(x = factor(org, levels = c("YL58", "I49", "I48", "I46", "YL45", "YL44","YL32" , "YL31", "YL27", "KB18","YL2" ,"KB1")), y = value, 
                           fill = factor(org, levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32","YL44" , "YL45", "I46", "I48","I49" ,"YL58")))) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip(ylim = c(-0.2, 0.2)) + # Flip coordinates for horizontal bars
  scale_y_continuous(labels = abs) + # Show positive labels on both sides
  geom_hline(yintercept = 0, color = "black")+ 
  labs(
    title = "Fraction of regulated proteins with significant rank change > 50 + exclusive proteins",
    x = "species",
    y = "fraction of regulated proteins in all detected proteins",
    caption = c("up in community", "up in monoculture")
  ) +
  theme_bw(base_size =14) + theme(legend.title=element_blank(), 
                                  plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 100, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"),
                                  plot.title = element_text(vjust = 10, size = 14))+
  scale_fill_manual(values = colors)

ggsave(filename = "fraction rank change 50 + exclusive_slim.pdf" , path =  "./enrichment_plots_static",
       width = 5, height = 6, units = "in")
ggsave(filename = "fraction rank change 50 + exclusive_slim.png" , path =  "./enrichment_plots_static",
       width = 5, height = 6, units = "in")

rm(all_discoveries_down, all_discoveries_up, cutoff_50_down, cutoff_50_up, cutoff_50_two_sided, perc_two_sided, protein_numbers)
