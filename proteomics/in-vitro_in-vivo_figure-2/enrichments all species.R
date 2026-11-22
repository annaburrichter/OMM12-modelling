library(tidyverse)
library(clusterProfiler)

prism_file <- get_analysis_tables("all_normalised.prism")

t_APF_AF <- prism_file[[1]]%>%
  rename(locus_tag = X1)%>%
  left_join(annotations, relationship = "many-to-many")%>%
  left_join(data,relationship = "many-to-many")%>%
  select(locus_tag, org, KEGG_ko, koeggNOG, `Discovery?`, `P value`, `Mean of M10_1`, `Mean of AAM_1`, Difference,
         `SE of difference`, `q value`)%>%
  dplyr::filter(!is.na(org))

t_mouse_AF <- prism_file[[2]]%>%
rename(locus_tag = X1)%>%
  left_join(annotations, relationship = "many-to-many")%>%
  left_join(data,relationship = "many-to-many")%>%
  select(locus_tag, org, KEGG_ko, koeggNOG, `Discovery?`, `P value`, `Mean of mouse_1`, `Mean of AAM_1`, Difference,
         `SE of difference`, `q value`)%>%
  dplyr::filter(!is.na(org))

t_mouse_APF <- prism_file[[3]]%>%
rename(locus_tag = X1)%>%
  left_join(annotations, relationship = "many-to-many")%>%
  left_join(data,relationship = "many-to-many")%>%
  select(locus_tag, org, KEGG_ko, koeggNOG, `Discovery?`, `P value`, `Mean of mouse_1`, `Mean of M10_1`, Difference,
         `SE of difference`, `q value`)%>%
  dplyr::filter(!is.na(org))



##enrichment mouse vs. APF
enrich_mouse_APF <- dplyr::filter(t_mouse_APF, `q value` <= 0.05)
exclusive_mouse <- dplyr::filter(t_mouse_APF, `Mean of mouse_1` >= 0.0001 & is.na(`Mean of M10_1`))
exclusive_APF <- dplyr::filter(t_mouse_APF, is.na(`Mean of mouse_1`) & `Mean of M10_1` >= 0.0001)

#up
mouse_APF_up <- enrich_mouse_APF %>%
  mutate(fold_change = `Mean of mouse_1`/`Mean of M10_1`)%>%
  dplyr::filter(fold_change >= 2)%>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_up <- bind_rows(mouse_APF_up, exclusive_mouse)
mouse_APF_up_ko <- mouse_APF_up$KEGG_ko
mouse_APF_up_ko <- str_replace_all(mouse_APF_up_ko, c("ko:" = ""))

enrichment_mouse_APF_up <- enrichKEGG(mouse_APF_up_ko, organism = "ko")

plot_mouse_APF_up <- enrichment_mouse_APF_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::filter(Description != "Carbon metabolism")%>%
  dplyr::filter(Description != "Ribosome")%>%
  dplyr::slice(1:10)

#down

mouse_APF_down <- enrich_mouse_APF %>%
  mutate(fold_change = `Mean of mouse_1`/`Mean of M10_1`)%>%
  dplyr::filter(fold_change <= 0.5)%>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_down <- bind_rows(mouse_APF_down, exclusive_APF)
mouse_APF_down_ko <- mouse_APF_down$KEGG_ko
mouse_APF_down_ko <- str_replace_all(mouse_APF_down_ko, c("ko:" = ""))

enrichment_mouse_APF_down <- enrichKEGG(mouse_APF_down_ko, organism = "ko")

plot_mouse_APF_down <- enrichment_mouse_APF_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::filter(Description != "Carbon metabolism")%>%
  dplyr::filter(Description != "Ribosome")%>%
  dplyr::slice(1:10)    
plot_mouse_APF_down$gene_ratio <-plot_mouse_APF_down$gene_ratio*(-1)

plot_mouse_APF <- dplyr::bind_rows(plot_mouse_APF_up, plot_mouse_APF_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,75)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. APF", caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. APF + exclusive.pdf", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. APF + exclusive.png", width = 9.86, height = 6.15, unit = "in")

##enrichment mouse vs. AF
enrich_mouse_AF <- dplyr::filter(t_mouse_AF, `q value` <= 0.05)
exclusive_mouse <- dplyr::filter(t_mouse_AF, `Mean of mouse_1` >= 0.0001 & is.na(`Mean of AAM_1`))
exclusive_AF <- dplyr::filter(t_mouse_AF, is.na(`Mean of mouse_1`) & `Mean of AAM_1` >= 0.0001)

#up
mouse_AF_up <- enrich_mouse_AF %>%
  mutate(fold_change = `Mean of mouse_1`/`Mean of AAM_1`)%>%
  dplyr::filter(fold_change >= 2)%>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_AF_up <- bind_rows(mouse_AF_up, exclusive_mouse)
mouse_AF_up_ko <- mouse_AF_up$KEGG_ko
mouse_AF_up_ko <- str_replace_all(mouse_AF_up_ko, c("ko:" = ""))

enrichment_mouse_AF_up <- enrichKEGG(mouse_AF_up_ko, organism = "ko")

plot_mouse_AF_up <- enrichment_mouse_AF_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::filter(Description != "Carbon metabolism")%>%
  dplyr::filter(Description != "Ribosome")%>%
  dplyr::slice(1:10)

#down

mouse_AF_down <- enrich_mouse_AF %>%
  mutate(fold_change = `Mean of mouse_1`/`Mean of AAM_1`)%>%
  dplyr::filter(fold_change <= 0.5)%>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_AF_down <- bind_rows(mouse_AF_down, exclusive_AF)
mouse_AF_down_ko <- mouse_AF_down$KEGG_ko
mouse_AF_down_ko <- str_replace_all(mouse_AF_down_ko, c("ko:" = ""))

enrichment_mouse_AF_down <- enrichKEGG(mouse_AF_down_ko, organism = "ko")

plot_mouse_AF_down <- enrichment_mouse_AF_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::filter(Description != "Carbon metabolism")%>%
  dplyr::filter(Description != "Ribosome")%>%
  dplyr::slice(1:10)    
plot_mouse_AF_down$gene_ratio <-plot_mouse_AF_down$gene_ratio*(-1)

plot_mouse_AF <- dplyr::bind_rows(plot_mouse_AF_up, plot_mouse_AF_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_AF, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,75)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. AF", caption = c("up in mouse", "up in AF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. AF + exclusive.pdf", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. AF + exclusive.png", width = 9.86, height = 6.15, unit = "in")

##enrichment APF vs. AF
enrich_APF_AF <- dplyr::filter(t_APF_AF, `q value` <= 0.05)
exclusive_APF <- dplyr::filter(t_APF_AF, `Mean of M10_1` >= 0.0001 & is.na(`Mean of AAM_1`))
exclusive_AF <- dplyr::filter(t_APF_AF, is.na(`Mean of M10_1`) & `Mean of AAM_1` >= 0.0001)

#up
APF_AF_up <- enrich_APF_AF %>%
  mutate(fold_change = `Mean of M10_1`/`Mean of AAM_1`)%>%
  dplyr::filter(fold_change >= 2)%>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
APF_AF_up <- bind_rows(APF_AF_up, exclusive_APF)
APF_AF_up_ko <- APF_AF_up$KEGG_ko
APF_AF_up_ko <- str_replace_all(APF_AF_up_ko, c("ko:" = ""))

enrichment_APF_AF_up <- enrichKEGG(APF_AF_up_ko, organism = "ko")

plot_APF_AF_up <- enrichment_APF_AF_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::filter(Description != "Carbon metabolism")%>%
  dplyr::filter(Description != "Ribosome")%>%
  dplyr::slice(1:10)

#down

APF_AF_down <- enrich_APF_AF %>%
  mutate(fold_change = `Mean of M10_1`/`Mean of AAM_1`)%>%
  dplyr::filter(fold_change <= 0.5)%>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
APF_AF_down <- bind_rows(APF_AF_down, exclusive_AF)
APF_AF_down_ko <- APF_AF_down$KEGG_ko
APF_AF_down_ko <- str_replace_all(APF_AF_down_ko, c("ko:" = ""))

enrichment_APF_AF_down <- enrichKEGG(APF_AF_down_ko, organism = "ko")

plot_APF_AF_down <- enrichment_APF_AF_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::filter(Description != "Carbon metabolism")%>%
  dplyr::filter(Description != "Ribosome")%>%
  dplyr::slice(1:10)    
plot_APF_AF_down$gene_ratio <-plot_APF_AF_down$gene_ratio*(-1)

plot_APF_AF <- dplyr::bind_rows(plot_APF_AF_up, plot_APF_AF_down) %>% arrange (desc(gene_ratio))

ggplot(plot_APF_AF, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,75)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "APF vs. AF", caption = c("up in APF", "up in AF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment APF vs. AF + exclusive.pdf", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment APF vs. AF + exclusive.png", width = 9.86, height = 6.15, unit = "in")

write.csv(enrichment_mouse_APF_up@result, file = "enrichment_mouse_APF_exclusive_up.csv")
write.csv(enrichment_mouse_APF_down@result, file = "enrichment_mouse_APF_exclusive_down.csv")
write.csv(enrichment_mouse_AF_up@result, file = "enrichment_mouse_AF_exclusive_up.csv")
write.csv(enrichment_mouse_AF_down@result, file = "enrichment_mouse_AF_exclusive_down.csv")
write.csv(enrichment_APF_AF_up@result, file = "enrichment_APF_AF_exclusive_up.csv")
write.csv(enrichment_APF_AF_down@result, file = "enrichment_APF_AF_exclusive_down.csv")

rm(plot_APF_AF, plot_APF_AF_down, plot_APF_AF_up, plot_mouse_AF, plot_mouse_AF_down, plot_mouse_AF_up,
   plot_mouse_APF, plot_mouse_APF_down, plot_mouse_APF_up, enrich_APF_AF, enrich_mouse_AF, enrich_mouse_APF,
   enrichment_APF_AF_down, enrichment_APF_AF_up, enrichment_mouse_AF_down, enrichment_mouse_AF_up,
   enrichment_mouse_APF_down, enrichment_mouse_APF_up, exclusive_AF, exclusive_APF, exclusive_mouse)
