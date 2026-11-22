###run enrichment over all species first to be able to reuse filtered tables
library(openxlsx)
###enrichment mouse vs. APF per species (no data for KB18, YL2, I49)

species_filtered <- c("KB1","YL27", "YL31", "YL32", "YL44", "YL45", "I46" , "I48" , "YL58")
wb <- createWorkbook()

for (i in 1:length(species_filtered)){
#up
mouse_APF_spec_up <- mouse_APF_up %>%
  dplyr::filter(org == species_filtered[i])%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_spec_up <- mouse_APF_spec_up$KEGG_ko
mouse_APF_spec_up <- str_replace_all(mouse_APF_spec_up, c("ko:" = ""))

if (length(mouse_APF_spec_up > 0)){
           
enrichment_mouse_APF_spec_up <- enrichKEGG(mouse_APF_spec_up, organism = "ko")

plot_mouse_APF_spec_up <- enrichment_mouse_APF_spec_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)
}else{
  plot_mouse_APF_spec_up <- data.frame(category = character(), subcategory = character(), ID = character(),
                                       Description = character(), genes = double(), all_genes = double(), 
                                       GeneRatio = character(), BgRatio = character(), RichFactor = double(),
                                       FoldEnrichment = double(), zScore = double(), pvalue = double(), 
                                       `p.adjust` = double(), qvalue = double(), geneID = character(), 
                                       Count = integer(), gene_ratio = double())
}
  
#down

mouse_APF_spec_down <- mouse_APF_down %>%
  dplyr::filter(org == species_filtered[i])%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_spec_down <- mouse_APF_spec_down$KEGG_ko
mouse_APF_spec_down <- str_replace_all(mouse_APF_spec_down, c("ko:" = ""))

if (length(mouse_APF_spec_down > 0)){
enrichment_mouse_APF_spec_down <- enrichKEGG(mouse_APF_spec_down, organism = "ko")

plot_mouse_APF_spec_down <- enrichment_mouse_APF_spec_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)    
plot_mouse_APF_spec_down$gene_ratio <-plot_mouse_APF_spec_down$gene_ratio*(-1)
}else{
  plot_mouse_APF_spec_down <- data.frame(category = character(), subcategory = character(), ID = character(),
                                       Description = character(), genes = double(), all_genes = double(), 
                                       GeneRatio = character(), BgRatio = character(), RichFactor = double(),
                                       FoldEnrichment = double(), zScore = double(), pvalue = double(), 
                                       `p.adjust` = double(), qvalue = double(), geneID = character(), 
                                       Count = integer(), gene_ratio = double())
}

plot_mouse_APF_spec <- dplyr::bind_rows(plot_mouse_APF_spec_up, plot_mouse_APF_spec_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF_spec, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = paste("mouse vs. APF for ", species_filtered[i]), caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = paste("enrichment mouse vs. APF + exclusive - ", species_filtered[i], ".pdf"), path = "./enrichment by species + exclusive", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = paste("enrichment mouse vs. APF + exclusive - ", species_filtered[i], ".png"), path = "./enrichment by species + exclusive", width = 9.86, height = 6.15, unit = "in")

if (length(mouse_APF_spec_up > 0)){
  addWorksheet(wb, paste(species_filtered[i], " APF - mouse up"))
writeData(wb, paste(species_filtered[i], " APF - mouse up"), enrichment_mouse_APF_spec_up@result)}
if (length(mouse_APF_spec_down > 0)){
addWorksheet(wb, paste(species_filtered[i], " APF - mouse down"))
writeData(wb, paste(species_filtered[i], " APF - mouse down"), enrichment_mouse_APF_spec_down@result)}


rm(mouse_APF_spec_down, mouse_APF_spec_up, plot_mouse_APF_spec, plot_mouse_APF_spec_down, plot_mouse_APF_spec_up, enrichment_mouse_APF_spec_down, enrichment_mouse_APF_spec_up)
}


# Save the workbook to an Excel file
saveWorkbook(wb, "Table S4 full datasets for enrichment in vitro - in vivo.xlsx", overwrite = TRUE)
