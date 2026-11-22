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


###volcano plot of APF vs. AF
t_APF_AF <- t_APF_AF %>%
            rename(locus_tag = X1)%>%
            left_join(annotations, relationship = "many-to-many")%>%
            left_join(data,relationship = "many-to-many")%>%
            select(locus_tag, org, KEGG_ko, koeggNOG, `Discovery?`, `P value`, `Mean of M10_1`, `Mean of AAM_1`, Difference,
                   `SE of difference`, `q value`)%>%
            dplyr::filter(!is.na(org))

volcano_APF_AF <- t_APF_AF %>%
                  mutate(log_q = -log(`q value`))%>%
                  mutate(log_fold_change = log(`Mean of M10_1`/`Mean of AAM_1`))
  
volcano_APF_AF$colorhelp <- ifelse(volcano_APF_AF$log_q > 1.301029995664e0, volcano_APF_AF$org, "ns")
volcano_APF_AF$color <- ifelse(volcano_APF_AF$log_fold_change < -0.30102999566398119521373889472449|volcano_APF_AF$log_fold_change > 0.30102999566398119521373889472449, volcano_APF_AF$colorhelp, "ns")
volcano_APF_AF$color <- factor(volcano_APF_AF$color, levels = c("KB1", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58", "ns"))
  

colors = c("#0e8b18", "#97b497", "#cc6a00", "#7bc57b","#11a196", "#6b00bc", "#871900","#365934", "#ff9b05", "#0b6760","lightgrey")
labs <- expression(italic("E. faecalis")~"KB1", italic("A. muris")~"KB18", italic("M. intestinale")~"YL27",
                                  italic("F. plautii")~"YL31", italic("E. clostridioformis")~"YL32", italic("A. muciniphila")~"YL44", italic("T. muris")~"YL45",
                                  italic("C. innocuum")~"I46", italic("B.caecimuris")~"I48", italic("B.coccoides")~"YL58", "", "")


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


###volcano plot of mouse vs. AF
t_mouse_AF <- t_mouse_AF %>%
            rename(locus_tag = X1)%>%
            left_join(annotations, relationship = "many-to-many")%>%
            left_join(data,relationship = "many-to-many")%>%
            select(locus_tag, org, KEGG_ko, koeggNOG, `Discovery?`, `P value`, `Mean of mouse_1`, `Mean of AAM_1`, Difference,
                   `SE of difference`, `q value`)%>%
            dplyr::filter(!is.na(org))

volcano_mouse_AF <- t_mouse_AF %>%
                  mutate(log_q = -log(`q value`))%>%
                  mutate(log_fold_change = log(`Mean of mouse_1`/`Mean of AAM_1`))
  
volcano_mouse_AF$colorhelp <- ifelse(volcano_mouse_AF$log_q > 1.301029995664e0, volcano_mouse_AF$org, "ns")
volcano_mouse_AF$color <- ifelse(volcano_mouse_AF$log_fold_change < -0.30102999566398119521373889472449|volcano_mouse_AF$log_fold_change > 0.30102999566398119521373889472449, volcano_mouse_AF$colorhelp, "ns")
volcano_mouse_AF$color <- factor(volcano_mouse_AF$color, levels = c("KB1",  "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58", "ns"))
  

colors = c("#0e8b18", "#cc6a00", "#7bc57b","#11a196", "#6b00bc", "#871900","#365934", "#ff9b05", "#0b6760","lightgrey")
labs <- expression(italic("E. faecalis")~"KB1", italic("M. intestinale")~"YL27",
                                  italic("F. plautii")~"YL31", italic("E. clostridioformis")~"YL32", italic("A. muciniphila")~"YL44", italic("T. muris")~"YL45",
                                  italic("C. innocuum")~"I46", italic("B.caecimuris")~"I48", italic("B.coccoides")~"YL58", "", "")


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

###volcano plot of mouse vs. APF
t_mouse_APF <- t_mouse_APF %>%
            rename(locus_tag = X1)%>%
            left_join(annotations, relationship = "many-to-many")%>%
            left_join(data,relationship = "many-to-many")%>%
            select(locus_tag, org, KEGG_ko, koeggNOG, `Discovery?`, `P value`, `Mean of mouse_1`, `Mean of M10_1`, Difference,
                   `SE of difference`, `q value`)%>%
            dplyr::filter(!is.na(org))

volcano_mouse_APF <- t_mouse_APF %>%
                  mutate(log_q = -log(`q value`))%>%
                  mutate(log_fold_change = log(`Mean of mouse_1`/`Mean of M10_1`))
  
volcano_mouse_APF$colorhelp <- ifelse(volcano_mouse_APF$log_q > 1.301029995664e0, volcano_mouse_APF$org, "ns")
volcano_mouse_APF$color <- ifelse(volcano_mouse_APF$log_fold_change < -0.30102999566398119521373889472449|volcano_mouse_APF$log_fold_change > 0.30102999566398119521373889472449, volcano_mouse_APF$colorhelp, "ns")
volcano_mouse_APF$color <- factor(volcano_mouse_APF$color, levels = c("KB1", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58", "ns"))
  

colors = c("#0e8b18", "#97b497", "#cc6a00", "#7bc57b","#11a196", "#6b00bc", "#871900","#365934", "#ff9b05", "#0b6760","lightgrey")
labs <- expression(italic("E. faecalis")~"KB1", italic("A. muris")~"KB18", italic("M. intestinale")~"YL27",
                                  italic("F. plautii")~"YL31", italic("E. clostridioformis")~"YL32", italic("A. muciniphila")~"YL44", italic("T. muris")~"YL45",
                                  italic("C. innocuum")~"I46", italic("B.caecimuris")~"I48", italic("B.coccoides")~"YL58", "", "")


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

rm(colors, labs, t_mouse_APF, t_mouse_AF, t_APF_AF, volcano_APF_AF, volcano_mouse_AF, volcano_mouse_APF)
