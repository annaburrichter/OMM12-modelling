library(tidyverse)
library(readxl)
library(prism2R)
library(clusterProfiler)

## import analysis from PRISM

prism_file <- get_analysis_tables("all_normalised.prism")
t_APF_AF <- prism_file[[1]]
t_mouse_AF <- prism_file[[2]]
t_mouse_APF <- prism_file[[3]]

colors <- c("#ff9b05","#365934","#004000","#0000cc","#cc6a00","#7bc57b","#11a196","#6b00bc","#871900","#0b6760", "#0e8b18", "#97b497")


###enrichment and volcano plot of APF vs. AF
t_APF_AF <- t_APF_AF %>%
            rename(locus_tag = X1)%>%
            left_join(annotations, relationship = "many-to-many")%>%
            left_join(data,relationship = "many-to-many")%>%
            select(locus_tag, org, KEGG_ko, koeggNOG, `Discovery?`, `P value`, `Mean of M10_1`, `Mean of AAM_1`, Difference,
                   `SE of difference`, `q value`)%>%
            dplyr::filter(!is.na(org))

##volcano plot
volcano_APF_AF <- t_APF_AF %>%
                  mutate(log_q = -log(`q value`))%>%
                  mutate(log_fold_change = log(`Mean of M10_1`/`Mean of AAM_1`))
  
volcano_APF_AF$colorhelp <- ifelse(volcano_APF_AF$log_q > 1.301029995664e0, volcano_APF_AF$org, "ns")
volcano_APF_AF$color <- ifelse(volcano_APF_AF$log_fold_change < -0.30102999566398119521373889472449|volcano_APF_AF$log_fold_change > 0.30102999566398119521373889472449, volcano_APF_AF$colorhelp, "ns")
volcano_APF_AF$color <- factor(volcano_APF_AF$color, levels = c("KB1", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58", "ns"))
  

colors = c("#0e8b18", "#97b497", "#cc6a00", "#7bc57b","#11a196", "#6b00bc", "#871900","#365934", "#ff9b05", "#0b6760","lightgrey")
labs <- expression(italic("E. faecalis")~"KB1", italic("A. muris")~"KB18", italic("M. intestinale")~"YL27",
                                  italic("F. plautii")~"YL31", italic("E. clostridioformis")~"YL32", italic("A. muciniphila")~"YL44", italic("T. muris")~"YL45",
                                  italic("C. innocuum")~"I46", italic("B.caecimuris")~"I48", italic("B.coccoides")~"YL58", "ns", "ns")


ggplot(volcano_APF_AF, aes(x=log_fold_change, y=log_q, color = color)) + 
  geom_point(shape=16) + 
  scale_color_manual(values= colors , labels = labs) +  
  geom_vline(xintercept = 0, color = "black") + 
  geom_vline(xintercept = 0.30102999566398119521373889472449, color = "grey") + 
  geom_vline(xintercept = -0.30102999566398119521373889472449, color = "grey") + 
  geom_hline(yintercept = 1.301029995664e0, color = "grey") +
  coord_cartesian(xlim = c(-5, 5)) + 
  theme_classic(base_size = 14) + 
  theme(axis.text = element_text(size=12)) + theme(axis.title.x = element_text(size = 12, colour = "grey30")) +
  labs(x = "log fold change", y = "-log(q value)", title = "APF vs AF - all strains", caption = c("up in APF", "up in AF"),color = "") + 
  theme(plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14))

ggsave(filename = "volcano_all_APF_AF.pdf", width = 9.86, height = 6.60, unit = "in")
ggsave(filename = "volcano_all_APF_AF.png", width = 9.86, height = 6.35, unit = "in")

##enrichment APF vs. AF
enrich_APF_AF <- dplyr::filter(t_APF_AF, `q value` <= 0.05)
#up
APF_AF_up <- enrich_APF_AF %>%
  mutate(fold_change = `Mean of M10_1`/`Mean of AAM_1`)%>%
  dplyr::filter(fold_change >= 2)%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
APF_AF_up <- APF_AF_up$KEGG_ko
APF_AF_up <- str_replace_all(APF_AF_up, c("ko:" = ""))

enrichment_APF_AF_up <- enrichKEGG(APF_AF_up, organism = "ko")

plot_APF_AF_up <- enrichment_APF_AF_up@result %>% 
    separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
    mutate(across(genes:all_genes, as.numeric))%>%
    mutate(gene_ratio = genes/all_genes)%>%
    dplyr::filter(p.adjust < 0.05)%>%
    dplyr::slice(1:10)

#down
    
APF_AF_down <- enrich_APF_AF %>%
  mutate(fold_change = `Mean of M10_1`/`Mean of AAM_1`)%>%
  dplyr::filter(fold_change <= 0.5)%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
APF_AF_down <- APF_AF_down$KEGG_ko
APF_AF_down <- str_replace_all(APF_AF_down, c("ko:" = ""))

enrichment_APF_AF_down <- enrichKEGG(APF_AF_down, organism = "ko")

plot_APF_AF_down <- enrichment_APF_AF_down@result %>% 
    separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
    mutate(across(genes:all_genes, as.numeric))%>%
    mutate(gene_ratio = genes/all_genes)%>%
    dplyr::filter(p.adjust < 0.05)%>%
    dplyr::slice(1:10)    
plot_APF_AF_down$gene_ratio <-plot_APF_AF_down$gene_ratio*(-1)
    
plot_APF_AF <- dplyr::bind_rows(plot_APF_AF_up, plot_APF_AF_down) %>% arrange (desc(gene_ratio))

ggplot(plot_APF_AF, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "APF vs. AF", caption = c("up in APF", "up in AF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment APF vs. AF.pdf", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment APF vs. AF.png", width = 9.86, height = 6.15, unit = "in")


###enrichment and volcano plot of mouse vs. AF
t_mouse_AF <- t_mouse_AF %>%
            rename(locus_tag = X1)%>%
            left_join(annotations, relationship = "many-to-many")%>%
            left_join(data,relationship = "many-to-many")%>%
            select(locus_tag, org, KEGG_ko, koeggNOG, `Discovery?`, `P value`, `Mean of mouse_1`, `Mean of AAM_1`, Difference,
                   `SE of difference`, `q value`)%>%
            dplyr::filter(!is.na(org))

##volcano plot
volcano_mouse_AF <- t_mouse_AF %>%
                  mutate(log_q = -log(`q value`))%>%
                  mutate(log_fold_change = log(`Mean of mouse_1`/`Mean of AAM_1`))
  
volcano_mouse_AF$colorhelp <- ifelse(volcano_mouse_AF$log_q > 1.301029995664e0, volcano_mouse_AF$org, "ns")
volcano_mouse_AF$color <- ifelse(volcano_mouse_AF$log_fold_change < -0.30102999566398119521373889472449|volcano_mouse_AF$log_fold_change > 0.30102999566398119521373889472449, volcano_mouse_AF$colorhelp, "ns")
volcano_mouse_AF$color <- factor(volcano_mouse_AF$color, levels = c("KB1",  "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58", "ns"))
  

colors = c("#0e8b18", "#cc6a00", "#7bc57b","#11a196", "#6b00bc", "#871900","#365934", "#ff9b05", "#0b6760","lightgrey")
labs <- expression(italic("E. faecalis")~"KB1", italic("M. intestinale")~"YL27",
                                  italic("F. plautii")~"YL31", italic("E. clostridioformis")~"YL32", italic("A. muciniphila")~"YL44", italic("T. muris")~"YL45",
                                  italic("C. innocuum")~"I46", italic("B.caecimuris")~"I48", italic("B.coccoides")~"YL58", "ns", "ns")


ggplot(volcano_mouse_AF, aes(x=log_fold_change, y=log_q, color = color)) + 
  geom_point(shape=16) + 
  scale_color_manual(values= colors , labels = labs) +  
  geom_vline(xintercept = 0, color = "black") + 
  geom_vline(xintercept = 0.30102999566398119521373889472449, color = "grey") + 
  geom_vline(xintercept = -0.30102999566398119521373889472449, color = "grey") + 
  geom_hline(yintercept = 1.301029995664e0, color = "grey") +
  coord_cartesian(xlim = c(-5, 5)) + 
  theme_classic(base_size = 14) + 
  theme(axis.text = element_text(size=12)) + theme(axis.title.x = element_text(size = 12, colour = "grey30")) +
  labs(x = "log fold change", y = "-log(q value)", title = "mouse vs AF - all strains", caption = c("up in mouse", "up in AF"),color = "") + 
  theme(plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14))

ggsave(filename = "volcano_all_mouse_AF.pdf", width = 9.86, height = 6.60, unit = "in")
ggsave(filename = "volcano_all_mouse_AF.png", width = 9.86, height = 6.35, unit = "in")


##enrichment mouse vs. AF
enrich_mouse_AF <- dplyr::filter(t_mouse_AF, `q value` <= 0.05)
#up
mouse_AF_up <- enrich_mouse_AF %>%
  mutate(fold_change = `Mean of mouse_1`/`Mean of AAM_1`)%>%
  dplyr::filter(fold_change >= 2)%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_AF_up <- mouse_AF_up$KEGG_ko
mouse_AF_up <- str_replace_all(mouse_AF_up, c("ko:" = ""))

enrichment_mouse_AF_up <- enrichKEGG(mouse_AF_up, organism = "ko")

plot_mouse_AF_up <- enrichment_mouse_AF_up@result %>% 
    separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
    mutate(across(genes:all_genes, as.numeric))%>%
    mutate(gene_ratio = genes/all_genes)%>%
    dplyr::filter(p.adjust < 0.05)%>%
    dplyr::slice(1:10)

#down
    
mouse_AF_down <- enrich_mouse_AF %>%
  mutate(fold_change = `Mean of mouse_1`/`Mean of AAM_1`)%>%
  dplyr::filter(fold_change <= 0.5)%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_AF_down <- mouse_AF_down$KEGG_ko
mouse_AF_down <- str_replace_all(mouse_AF_down, c("ko:" = ""))

enrichment_mouse_AF_down <- enrichKEGG(mouse_AF_down, organism = "ko")

plot_mouse_AF_down <- enrichment_mouse_AF_down@result %>% 
    separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
    mutate(across(genes:all_genes, as.numeric))%>%
    mutate(gene_ratio = genes/all_genes)%>%
    dplyr::filter(p.adjust < 0.05)%>%
    dplyr::slice(1:10)    
plot_mouse_AF_down$gene_ratio <-plot_mouse_AF_down$gene_ratio*(-1)
    
plot_mouse_AF <- dplyr::bind_rows(plot_mouse_AF_up, plot_mouse_AF_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_AF, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. AF", caption = c("up in mouse", "up in AF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. AF.pdf", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. AF.png", width = 9.86, height = 6.15, unit = "in")


###enrichment and volcano plot of mouse vs. APF
t_mouse_APF <- t_mouse_APF %>%
            rename(locus_tag = X1)%>%
            left_join(annotations, relationship = "many-to-many")%>%
            left_join(data,relationship = "many-to-many")%>%
            select(locus_tag, org, KEGG_ko, koeggNOG, `Discovery?`, `P value`, `Mean of mouse_1`, `Mean of M10_1`, Difference,
                   `SE of difference`, `q value`)%>%
            dplyr::filter(!is.na(org))

##volcano plot
volcano_mouse_APF <- t_mouse_APF %>%
                  mutate(log_q = -log(`q value`))%>%
                  mutate(log_fold_change = log(`Mean of mouse_1`/`Mean of M10_1`))
  
volcano_mouse_APF$colorhelp <- ifelse(volcano_mouse_APF$log_q > 1.301029995664e0, volcano_mouse_APF$org, "ns")
volcano_mouse_APF$color <- ifelse(volcano_mouse_APF$log_fold_change < -0.30102999566398119521373889472449|volcano_mouse_APF$log_fold_change > 0.30102999566398119521373889472449, volcano_mouse_APF$colorhelp, "ns")
volcano_mouse_APF$color <- factor(volcano_mouse_APF$color, levels = c("KB1", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58", "ns"))
  

colors = c("#0e8b18", "#97b497", "#cc6a00", "#7bc57b","#11a196", "#6b00bc", "#871900","#365934", "#ff9b05", "#0b6760","lightgrey")
labs <- expression(italic("E. faecalis")~"KB1", italic("A. muris")~"KB18", italic("M. intestinale")~"YL27",
                                  italic("F. plautii")~"YL31", italic("E. clostridioformis")~"YL32", italic("A. muciniphila")~"YL44", italic("T. muris")~"YL45",
                                  italic("C. innocuum")~"I46", italic("B.caecimuris")~"I48", italic("B.coccoides")~"YL58", "ns", "ns")


ggplot(volcano_mouse_APF, aes(x=log_fold_change, y=log_q, color = color)) + 
  geom_point(shape=16) + 
  scale_color_manual(values= colors , labels = labs) +  
  geom_vline(xintercept = 0, color = "black") + 
  geom_vline(xintercept = 0.30102999566398119521373889472449, color = "grey") + 
  geom_vline(xintercept = -0.30102999566398119521373889472449, color = "grey") + 
  geom_hline(yintercept = 1.301029995664e0, color = "grey") +
  coord_cartesian(xlim = c(-5, 5)) + 
  theme_classic(base_size = 14) + 
  theme(axis.text = element_text(size=12)) + theme(axis.title.x = element_text(size = 12, colour = "grey30")) +
  labs(x = "log fold change", y = "-log(q value)", title = "mouse vs APF - all strains", caption = c("up in mouse", "up in APF"),color = "") + 
  theme(plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14))

ggsave(filename = "volcano_all_mouse_APF.pdf", width = 9.86, height = 6.60, unit = "in")
ggsave(filename = "volcano_all_mouse_APF.png", width = 9.86, height = 6.35, unit = "in")


##enrichment mouse vs. APF
enrich_mouse_APF <- dplyr::filter(t_mouse_APF, `q value` <= 0.05)
#up
mouse_APF_up <- enrich_mouse_APF %>%
  mutate(fold_change = `Mean of mouse_1`/`Mean of M10_1`)%>%
  dplyr::filter(fold_change >= 2)%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_up <- mouse_APF_up$KEGG_ko
mouse_APF_up <- str_replace_all(mouse_APF_up, c("ko:" = ""))

enrichment_mouse_APF_up <- enrichKEGG(mouse_APF_up, organism = "ko")

plot_mouse_APF_up <- enrichment_mouse_APF_up@result %>% 
    separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
    mutate(across(genes:all_genes, as.numeric))%>%
    mutate(gene_ratio = genes/all_genes)%>%
    dplyr::filter(p.adjust < 0.05)%>%
    dplyr::slice(1:10)

#down
    
mouse_APF_down <- enrich_mouse_APF %>%
  mutate(fold_change = `Mean of mouse_1`/`Mean of M10_1`)%>%
  dplyr::filter(fold_change <= 0.5)%>%
  select(KEGG_ko) %>%
  dplyr::filter(!is.na(KEGG_ko)) %>%  dplyr::filter(KEGG_ko != "-") %>% dplyr::filter(KEGG_ko != "") %>% 
  separate_longer_delim(KEGG_ko, ",") %>% dplyr::distinct(KEGG_ko, .keep_all = TRUE)
mouse_APF_down <- mouse_APF_down$KEGG_ko
mouse_APF_down <- str_replace_all(mouse_APF_down, c("ko:" = ""))

enrichment_mouse_APF_down <- enrichKEGG(mouse_APF_down, organism = "ko")

plot_mouse_APF_down <- enrichment_mouse_APF_down@result %>% 
    separate_wider_delim(GeneRatio, delim = "/", names = c("genes", "all_genes"), cols_remove = FALSE)%>%
    mutate(across(genes:all_genes, as.numeric))%>%
    mutate(gene_ratio = genes/all_genes)%>%
    dplyr::filter(p.adjust < 0.05)%>%
    dplyr::slice(1:10)    
plot_mouse_APF_down$gene_ratio <-plot_mouse_APF_down$gene_ratio*(-1)
    
plot_mouse_APF <- dplyr::bind_rows(plot_mouse_APF_up, plot_mouse_APF_down) %>% arrange (desc(gene_ratio))

ggplot(plot_mouse_APF, aes(x=gene_ratio, y=reorder(Description, Count), size = Count, color = p.adjust), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) + 
  geom_point(shape=16) + 
  scale_color_gradient(low = "darkorchid4", high = "white", name = "p.adjust", limits = c(0, 0.06), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05), labels = c(0.01, 0.02, 0.03, 0.04, 0.05)) +
  scale_size_continuous(limits = c(0,51)) +
  coord_cartesian(xlim = c(-0.5, 0.5)) + 
  geom_vline(xintercept = 0, color = "grey") + 
  theme_bw(base_size = 12) + 
  labs(x = "gene ratio", y = "", title = "mouse vs. APF", caption = c("up in mouse", "up in APF")) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = 10), axis.text = element_text(size = 14), axis.title.x = element_text(size = 14, colour = "grey30"), 
        plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 135, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"))

ggsave(filename = "enrichment mouse vs. APF.pdf", width = 9.86, height = 6.50, unit = "in")
ggsave(filename = "enrichment mouse vs. APF.png", width = 9.86, height = 6.15, unit = "in")
