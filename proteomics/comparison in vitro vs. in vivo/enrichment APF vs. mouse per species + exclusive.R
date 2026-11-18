###enrichment mouse vs. APF per species (no data for KB18, YL2, I49)

##YL32
#up
mouse_APF_YL32_up <- mouse_APF_up %>%
  dplyr::filter(org == "YL32")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL32_up <- mouse_APF_YL32_up$KEGG_ko
mouse_APF_YL32_up <- str_replace_all(mouse_APF_YL32_up, c("ko:" = ""))

enrichment_mouse_APF_YL32_up <- enrichKEGG(mouse_APF_YL32_up, organism = "ko")

plot_mouse_APF_YL32_up <- enrichment_mouse_APF_YL32_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)

#down

mouse_APF_YL32_down <- mouse_APF_down %>%
  dplyr::filter(org == "YL32")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL32_down <- mouse_APF_YL32_down$KEGG_ko
mouse_APF_YL32_down <- str_replace_all(mouse_APF_YL32_down, c("ko:" = ""))

enrichment_mouse_APF_YL32_down <- enrichKEGG(mouse_APF_YL32_down, organism = "ko")

plot_mouse_APF_YL32_down <- enrichment_mouse_APF_YL32_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)    
plot_mouse_APF_YL32_down$gene_ratio <-plot_mouse_APF_YL32_down$gene_ratio*(-1)

plot_mouse_APF_YL32 <- dplyr::bind_rows(plot_mouse_APF_YL32_up, plot_mouse_APF_YL32_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF_YL32, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. APF for YL32", caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. APF + exclusive - YL32.pdf", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. APF + exclusive - YL32.png", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.15, unit = "in")


##YL31
#up
mouse_APF_YL31_up <- mouse_APF_up %>%
  dplyr::filter(org == "YL31")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL31_up <- mouse_APF_YL31_up$KEGG_ko
mouse_APF_YL31_up <- str_replace_all(mouse_APF_YL31_up, c("ko:" = ""))

enrichment_mouse_APF_YL31_up <- enrichKEGG(mouse_APF_YL31_up, organism = "ko")

plot_mouse_APF_YL31_up <- enrichment_mouse_APF_YL31_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)

#down

mouse_APF_YL31_down <- mouse_APF_down %>%
  dplyr::filter(org == "YL31")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL31_down <- mouse_APF_YL31_down$KEGG_ko
mouse_APF_YL31_down <- str_replace_all(mouse_APF_YL31_down, c("ko:" = ""))

enrichment_mouse_APF_YL31_down <- enrichKEGG(mouse_APF_YL31_down, organism = "ko")

plot_mouse_APF_YL31_down <- enrichment_mouse_APF_YL31_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)    
plot_mouse_APF_YL31_down$gene_ratio <-plot_mouse_APF_YL31_down$gene_ratio*(-1)

plot_mouse_APF_YL31 <- dplyr::bind_rows(plot_mouse_APF_YL31_up, plot_mouse_APF_YL31_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF_YL31, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. APF for YL31", caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. APF + exclusive - YL31.pdf", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. APF + exclusive - YL31.png", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.15, unit = "in")


##YL44
#up
mouse_APF_YL44_up <- mouse_APF_up %>%
  dplyr::filter(org == "YL44")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL44_up <- mouse_APF_YL44_up$KEGG_ko
mouse_APF_YL44_up <- str_replace_all(mouse_APF_YL44_up, c("ko:" = ""))

enrichment_mouse_APF_YL44_up <- enrichKEGG(mouse_APF_YL44_up, organism = "ko")

plot_mouse_APF_YL44_up <- enrichment_mouse_APF_YL44_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)

#down

mouse_APF_YL44_down <- mouse_APF_down %>%
  dplyr::filter(org == "YL44")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL44_down <- mouse_APF_YL44_down$KEGG_ko
mouse_APF_YL44_down <- str_replace_all(mouse_APF_YL44_down, c("ko:" = ""))

enrichment_mouse_APF_YL44_down <- enrichKEGG(mouse_APF_YL44_down, organism = "ko")

plot_mouse_APF_YL44_down <- enrichment_mouse_APF_YL44_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)    
plot_mouse_APF_YL44_down$gene_ratio <-plot_mouse_APF_YL44_down$gene_ratio*(-1)

plot_mouse_APF_YL44 <- dplyr::bind_rows(plot_mouse_APF_YL44_up, plot_mouse_APF_YL44_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF_YL44, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. APF for YL44", caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. APF + exclusive - YL44.pdf", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. APF + exclusive - YL44.png", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.15, unit = "in")


##YL45
#up
mouse_APF_YL45_up <- mouse_APF_up %>%
  dplyr::filter(org == "YL45")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL45_up <- mouse_APF_YL45_up$KEGG_ko
mouse_APF_YL45_up <- str_replace_all(mouse_APF_YL45_up, c("ko:" = ""))

enrichment_mouse_APF_YL45_up <- enrichKEGG(mouse_APF_YL45_up, organism = "ko")

plot_mouse_APF_YL45_up <- enrichment_mouse_APF_YL45_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)

#down

mouse_APF_YL45_down <- mouse_APF_down %>%
  dplyr::filter(org == "YL45")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL45_down <- mouse_APF_YL45_down$KEGG_ko
mouse_APF_YL45_down <- str_replace_all(mouse_APF_YL45_down, c("ko:" = ""))

enrichment_mouse_APF_YL45_down <- enrichKEGG(mouse_APF_YL45_down, organism = "ko")

plot_mouse_APF_YL45_down <- enrichment_mouse_APF_YL45_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)    
plot_mouse_APF_YL45_down$gene_ratio <-plot_mouse_APF_YL45_down$gene_ratio*(-1)

plot_mouse_APF_YL45 <- dplyr::bind_rows(plot_mouse_APF_YL45_up, plot_mouse_APF_YL45_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF_YL45, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. APF for YL45", caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. APF + exclusive - YL45.pdf", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. APF + exclusive - YL45.png", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.15, unit = "in")


##YL58
#up
mouse_APF_YL58_up <- mouse_APF_up %>%
  dplyr::filter(org == "YL58")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL58_up <- mouse_APF_YL58_up$KEGG_ko
mouse_APF_YL58_up <- str_replace_all(mouse_APF_YL58_up, c("ko:" = ""))

enrichment_mouse_APF_YL58_up <- enrichKEGG(mouse_APF_YL58_up, organism = "ko")

plot_mouse_APF_YL58_up <- enrichment_mouse_APF_YL58_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)

#down

mouse_APF_YL58_down <- mouse_APF_down %>%
  dplyr::filter(org == "YL58")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL58_down <- mouse_APF_YL58_down$KEGG_ko
mouse_APF_YL58_down <- str_replace_all(mouse_APF_YL58_down, c("ko:" = ""))

enrichment_mouse_APF_YL58_down <- enrichKEGG(mouse_APF_YL58_down, organism = "ko")

plot_mouse_APF_YL58_down <- enrichment_mouse_APF_YL58_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)    
plot_mouse_APF_YL58_down$gene_ratio <-plot_mouse_APF_YL58_down$gene_ratio*(-1)

plot_mouse_APF_YL58 <- dplyr::bind_rows(plot_mouse_APF_YL58_up, plot_mouse_APF_YL58_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF_YL58, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. APF for YL58", caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. APF + exclusive - YL58.pdf", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. APF + exclusive - YL58.png", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.15, unit = "in")


##YL27
mouse_APF_YL27_up <- mouse_APF_up %>%
  dplyr::filter(org == "YL27")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL27_up <- mouse_APF_YL27_up$KEGG_ko
mouse_APF_YL27_up <- str_replace_all(mouse_APF_YL27_up, c("ko:" = ""))

enrichment_mouse_APF_YL27_up <- enrichKEGG(mouse_APF_YL27_up, organism = "ko")

plot_mouse_APF_YL27_up <- enrichment_mouse_APF_YL27_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)

#down

mouse_APF_YL27_down <- enrich_mouse_APF %>%
  mutate(fold_change = `Mean of mouse_1`/`Mean of M10_1`)%>%
  dplyr::filter(fold_change <= 0.5)%>%
  dplyr::filter(org == "YL27")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_YL27_down <- mouse_APF_YL27_down$KEGG_ko
mouse_APF_YL27_down <- str_replace_all(mouse_APF_YL27_down, c("ko:" = ""))

enrichment_mouse_APF_YL27_down <- enrichKEGG(mouse_APF_YL27_down, organism = "ko")

plot_mouse_APF_YL27_down <- enrichment_mouse_APF_YL27_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)    
plot_mouse_APF_YL27_down$gene_ratio <-plot_mouse_APF_YL27_down$gene_ratio*(-1)

plot_mouse_APF_YL27 <- dplyr::bind_rows(plot_mouse_APF_YL27_up, plot_mouse_APF_YL27_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF_YL27, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. APF for YL27", caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. APF + exclusive - YL27.pdf", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. APF + exclusive - YL27.png", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.15, unit = "in")


##KB1
#up
mouse_APF_KB1_up <- mouse_APF_up %>%
  dplyr::filter(org == "KB1")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_KB1_up <- mouse_APF_KB1_up$KEGG_ko
mouse_APF_KB1_up <- str_replace_all(mouse_APF_KB1_up, c("ko:" = ""))

enrichment_mouse_APF_KB1_up <- enrichKEGG(mouse_APF_KB1_up, organism = "ko")

plot_mouse_APF_KB1_up <- enrichment_mouse_APF_KB1_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)

#down - no data

plot_mouse_APF_KB1_down <- data.frame(matrix(nrow=0, ncol=17))


plot_mouse_APF_KB1 <- dplyr::bind_rows(plot_mouse_APF_KB1_up, plot_mouse_APF_KB1_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF_KB1, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. APF for KB1", caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. APF + exclusive - KB1.pdf", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. APF + exclusive - KB1.png", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.15, unit = "in")


##I46
#up
mouse_APF_I46_up <- mouse_APF_up %>%
  dplyr::filter(org == "I46")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_I46_up <- mouse_APF_I46_up$KEGG_ko
mouse_APF_I46_up <- str_replace_all(mouse_APF_I46_up, c("ko:" = ""))

enrichment_mouse_APF_I46_up <- enrichKEGG(mouse_APF_I46_up, organism = "ko")

plot_mouse_APF_I46_up <- enrichment_mouse_APF_I46_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)

#down

mouse_APF_I46_down <- mouse_APF_down %>%
  dplyr::filter(org == "I46")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_I46_down <- mouse_APF_I46_down$KEGG_ko
mouse_APF_I46_down <- str_replace_all(mouse_APF_I46_down, c("ko:" = ""))

enrichment_mouse_APF_I46_down <- enrichKEGG(mouse_APF_I46_down, organism = "ko")

plot_mouse_APF_I46_down <- enrichment_mouse_APF_I46_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)    
plot_mouse_APF_I46_down$gene_ratio <-plot_mouse_APF_I46_down$gene_ratio*(-1)

plot_mouse_APF_I46 <- dplyr::bind_rows(plot_mouse_APF_I46_up, plot_mouse_APF_I46_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF_I46, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. APF for I46", caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. APF + exclusive - I46.pdf", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. APF + exclusive - I46.png", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.15, unit = "in")


##I48
#up
mouse_APF_I48_up <- mouse_APF_up %>%
  dplyr::filter(org == "I48")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_I48_up <- mouse_APF_I48_up$KEGG_ko
mouse_APF_I48_up <- str_replace_all(mouse_APF_I48_up, c("ko:" = ""))

enrichment_mouse_APF_I48_up <- enrichKEGG(mouse_APF_I48_up, organism = "ko")

plot_mouse_APF_I48_up <- enrichment_mouse_APF_I48_up@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)

#down

mouse_APF_I48_down <- mouse_APF_down %>%
  dplyr::filter(org == "I48")%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_I48_down <- mouse_APF_I48_down$KEGG_ko
mouse_APF_I48_down <- str_replace_all(mouse_APF_I48_down, c("ko:" = ""))

enrichment_mouse_APF_I48_down <- enrichKEGG(mouse_APF_I48_down, organism = "ko")

plot_mouse_APF_I48_down <- enrichment_mouse_APF_I48_down@result %>% 
  separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
  mutate(across(genes:all_genes, as.numeric))%>%
  mutate(gene_ratio = genes/all_genes)%>%
  dplyr::filter(p.adjust < 0.05)%>%
  dplyr::slice(1:10)    
plot_mouse_APF_I48_down$gene_ratio <-plot_mouse_APF_I48_down$gene_ratio*(-1)

plot_mouse_APF_I48 <- dplyr::bind_rows(plot_mouse_APF_I48_up, plot_mouse_APF_I48_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF_I48, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. APF for I48", caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. APF + exclusive - I48.pdf", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. APF + exclusive - I48.png", path = "C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/enrichment by species + exclusive", width = 9.86, height = 6.15, unit = "in")

write.csv(enrichment_mouse_APF_YL32_up@result, file = "enrichment_mouse_APF_YL32_exclusive_up.csv")
write.csv(enrichment_mouse_APF_YL32_down@result, file = "enrichment_mouse_APF_YL32_exclusive_down.csv")

write.csv(enrichment_mouse_APF_YL31_up@result, file = "enrichment_mouse_APF_YL31_exclusive_up.csv")
write.csv(enrichment_mouse_APF_YL31_down@result, file = "enrichment_mouse_APF_YL31_exclusive_down.csv")

write.csv(enrichment_mouse_APF_YL44_up@result, file = "enrichment_mouse_APF_YL44_exclusive_up.csv")
write.csv(enrichment_mouse_APF_YL44_down@result, file = "enrichment_mouse_APF_YL44_exclusive_down.csv")

write.csv(enrichment_mouse_APF_YL45_up@result, file = "enrichment_mouse_APF_YL45_exclusive_up.csv")
write.csv(enrichment_mouse_APF_YL45_down@result, file = "enrichment_mouse_APF_YL45_exclusive_down.csv")

write.csv(enrichment_mouse_APF_YL58_up@result, file = "enrichment_mouse_APF_YL58_exclusive_up.csv")
write.csv(enrichment_mouse_APF_YL58_down@result, file = "enrichment_mouse_APF_YL58_exclusive_down.csv")

write.csv(enrichment_mouse_APF_YL27_up@result, file = "enrichment_mouse_APF_YL27_exclusive_up.csv")
write.csv(enrichment_mouse_APF_YL27_down@result, file = "enrichment_mouse_APF_YL27_exclusive_down.csv")

write.csv(enrichment_mouse_APF_KB1_up@result, file = "enrichment_mouse_APF_KB1_exclusive_up.csv")

write.csv(enrichment_mouse_APF_I46_up@result, file = "enrichment_mouse_APF_I46_exclusive_up.csv")
write.csv(enrichment_mouse_APF_I46_down@result, file = "enrichment_mouse_APF_I46_exclusive_down.csv")

write.csv(enrichment_mouse_APF_I48_up@result, file = "enrichment_mouse_APF_I48_exclusive_up.csv")
write.csv(enrichment_mouse_APF_I48_down@result, file = "enrichment_mouse_APF_I48_exclusive_down.csv")


#export enrichment data to Excel file as supplementary data

# Create a new workbook
wb <- createWorkbook()

# Add worksheets and write data frames to them
addWorksheet(wb, "mouse - AF all species up")
writeData(wb, "mouse - AF all species up", enrichment_mouse_AF_up@result)
addWorksheet(wb, "mouse - AF all species down")
writeData(wb, "mouse - AF all species down", enrichment_mouse_AF_down@result)

addWorksheet(wb, "APF - AF all species up")
writeData(wb, "APF - AF all species up", enrichment_APF_AF_up@result)
addWorksheet(wb, "APF - AF all species down")
writeData(wb, "APF - AF all species down", enrichment_APF_AF_down@result)

addWorksheet(wb, "mouse - APF all species up")
writeData(wb, "mouse - APF all species up", enrichment_mouse_APF_up@result)
addWorksheet(wb, "mouse - APF all species down")
writeData(wb, "mouse - APF all species down", enrichment_mouse_APF_down@result)

addWorksheet(wb, "KB1 APF - mouse up")
writeData(wb, "KB1 APF - mouse up", enrichment_mouse_APF_KB1_up@result)
addWorksheet(wb, "KB1 APF - mouse down")
writeData(wb, "KB1 APF - mouse down", enrichment_mouse_APF_KB1_down)

addWorksheet(wb, "YL27 APF - mouse up")
writeData(wb, "YL27 APF - mouse up", enrichment_mouse_APF_YL27_up)
addWorksheet(wb, "YL27 APF - mouse down")
writeData(wb, "YL27 APF - mouse down", enrichment_mouse_APF_YL27_down@result)

addWorksheet(wb, "YL31 APF - mouse up")
writeData(wb, "YL31 APF - mouse up", enrichment_mouse_APF_YL31_up@result)
addWorksheet(wb, "YL31 APF - mouse down")
writeData(wb, "YL31 APF - mouse down", enrichment_mouse_APF_YL31_down@result)

addWorksheet(wb, "YL32 APF - mouse up")
writeData(wb, "YL32 APF - mouse up", enrichment_mouse_APF_YL32_up@result)
addWorksheet(wb, "YL32 APF - mouse down")
writeData(wb, "YL32 APF - mouse down", enrichment_mouse_APF_YL32_down@result)

addWorksheet(wb, "YL44 APF - mouse up")
writeData(wb, "YL44 APF - mouse up", enrichment_mouse_APF_YL44_up@result)
addWorksheet(wb, "YL44 APF - mouse down")
writeData(wb, "YL44 APF - mouse down", enrichment_mouse_APF_YL44_down@result)

addWorksheet(wb, "YL45 APF - mouse up")
writeData(wb, "YL45 APF - mouse up", enrichment_mouse_APF_YL45_up@result)
addWorksheet(wb, "YL45 APF - mouse down")
writeData(wb, "YL45 APF - mouse down", enrichment_mouse_APF_YL45_down@result)

addWorksheet(wb, "I46 APF - mouse up")
writeData(wb, "I46 APF - mouse up", enrichment_mouse_APF_I46_up@result)
addWorksheet(wb, "I46 APF - mouse down")
writeData(wb, "I46 APF - mouse down", enrichment_mouse_APF_I46_down@result)

addWorksheet(wb, "I48 APF - mouse up")
writeData(wb, "I48 APF - mouse up", enrichment_mouse_APF_I48_up@result)
addWorksheet(wb, "I48 APF - mouse down")
writeData(wb, "I48 APF - mouse down", enrichment_mouse_APF_I48_down@result)

addWorksheet(wb, "YL58 APF - mouse up")
writeData(wb, "YL58 APF - mouse up", enrichment_mouse_APF_YL58_up@result)
addWorksheet(wb, "YL58 APF - mouse down")
writeData(wb, "YL58 APF - mouse down", enrichment_mouse_APF_YL58_down@result)

# Save the workbook to an Excel file
saveWorkbook(wb, "Table S4 full datasets for enrichment in vitro - in vivo.xlsx", overwrite = TRUE)
