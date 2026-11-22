library(openxlsx)
library(tidyverse)
library(readxl)

annotations <- read.delim2("protein_annotation.tsv")
species_short <- c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46","I48", "I49", "YL58")


###for vitro-vivo datasets

mouse_vitro <- read_excel("reanalysed_mouse_vs_invitro.xlsx", 
                                  col_types = c("skip", "text", "text", "skip", "skip", "text", "text", "skip", 
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                                                "numeric", "numeric", "numeric"))%>%
                      separate_wider_delim(2, delim = " ", names = c("genus", "species", "strain"))



vitro_vivo_list_ranks <- list()
correlations <- data.frame(strain = character(), comparison = character(), rho = numeric(),p = numeric())

for (i in 1:12){

table <- mouse_vitro %>% 
  filter(strain == species_short[[i]])
  #if there is only one value for a protein, set this to NA also  
  cols <- c(7, 8, 9)  # AF
  table <- table %>% mutate(across(all_of(cols),
                                   ~ replace(.,rowSums(!is.na(pick(all_of(cols)))) < 2,NA)))
  cols <- c(10, 11, 12)  # APF
  table <- table %>% mutate(across(all_of(cols),
                                   ~ replace(.,rowSums(!is.na(pick(all_of(cols)))) < 2,NA)))
  
  cols <- c(13, 14, 15)  # mouse
  table <- table %>% mutate(across(all_of(cols),
                                   ~ replace(.,rowSums(!is.na(pick(all_of(cols)))) < 2,NA)))
  rm(cols)
  
  table <- table %>%
  mutate(rank_AF_1 = rank(-.[[7]], na.last = "keep")) %>%
  mutate(rank_AF_2 = rank(-.[[8]], na.last = "keep")) %>%
  mutate(rank_AF_3 = rank(-.[[9]], na.last = "keep")) %>%
  mutate(rank_APF_1 = rank(-.[[10]], na.last = "keep")) %>%
  mutate(rank_APF_2 = rank(-.[[11]], na.last = "keep")) %>%
  mutate(rank_APF_3 = rank(-.[[12]], na.last = "keep")) %>%
  mutate(rank_mouse_1 = rank(-.[[13]], na.last = "keep")) %>%
  mutate(rank_mouse_2 = rank(-.[[14]], na.last = "keep")) %>%
  mutate(rank_mouse_3 = rank(-.[[15]], na.last = "keep")) %>%
  mutate(mean_rank_AF = rowMeans(select(., rank_AF_1, rank_AF_2, rank_AF_3), na.rm = TRUE)) %>%
  mutate(mean_rank_APF = rowMeans(select(., rank_APF_1, rank_APF_2, rank_APF_3), na.rm = TRUE)) %>%
  mutate(mean_rank_mouse = rowMeans(select(., rank_mouse_1, rank_mouse_2, rank_mouse_3), na.rm = TRUE))%>%
  mutate(across(starts_with("mean_"), ~replace_na(., NA)))

axis_lim = max(c(table$mean_rank_AF, table$mean_rank_mouse), na.rm = TRUE)*1.1
p_AF_mouse<- ggplot(table, aes(x = mean_rank_AF, y = mean_rank_mouse)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  scale_y_continuous(limits = c(0, axis_lim), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, axis_lim), expand = c(0, 0)) +
  labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
       x = "Mean Rank in AF",
       y = "Mean Rank in Mouse") +
  theme_minimal()

axis_lim = max(c(table$mean_rank_AF, table$mean_rank_APF), na.rm = TRUE)*1.1
p_AF_APF<- ggplot(table, aes(x = mean_rank_AF, y = mean_rank_APF)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  scale_y_continuous(limits = c(0, axis_lim), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, axis_lim), expand = c(0, 0)) +
  labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
       x = "Mean Rank in AF",
       y = "Mean Rank in APF") +
  theme_minimal()

axis_lim = max(c(table$mean_rank_APF, table$mean_rank_mouse), na.rm = TRUE)*1.1
p_APF_mouse<- ggplot(table, aes(x = mean_rank_APF, y = mean_rank_mouse)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  scale_y_continuous(limits = c(0, axis_lim), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, axis_lim), expand = c(0, 0)) +
  labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
       x = "Mean Rank in APF",
       y = "Mean Rank in Mouse") +
  theme_minimal()

ggsave(  filename = paste0(paste(species_short[i], "AF", "mouse", sep = "_"), ".png"),
  plot = p_AF_mouse,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300)
ggsave(  filename = paste0(paste(species_short[i], "AF", "mouse", sep = "_"), ".pdf"),
  plot = p_AF_mouse,  bg = "white",  width = 8,  height = 6,  units = "in")

ggsave(  filename = paste0(paste(species_short[i], "APF", "mouse", sep = "_"), ".png"),
         plot = p_APF_mouse,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300)
ggsave(  filename = paste0(paste(species_short[i], "APF", "mouse", sep = "_"), ".pdf"),
         plot = p_APF_mouse,  bg = "white",  width = 8,  height = 6,  units = "in")

ggsave(  filename = paste0(paste(species_short[i], "AF", "APF", sep = "_"), ".png"),
         plot = p_AF_APF,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300)
ggsave(  filename = paste0(paste(species_short[i], "AF", "APF", sep = "_"), ".pdf"),
         plot = p_AF_APF,  bg = "white",  width = 8,  height = 6,  units = "in")


tmp <- select(table, mean_rank_AF, mean_rank_mouse)
  colnames(tmp) <- c("x", "y")  
tmp <- na.omit(tmp[, c("x", "y")])
  
if (nrow(tmp) < 3) {
    rho <- NA
    pval <- NA} 
else {
    test <- cor.test(tmp$x,tmp$y,method = "spearman",exact = FALSE)
    rho <- unname(test$estimate)
    pval <- test$p.value}
  
  correlations <- rbind(correlations,
    data.frame(strain = species_short[i],comparison = "AF_mouse", rho = rho, p = pval))
rm(tmp, rho, pval, test)

tmp <- select(table, mean_rank_APF, mean_rank_mouse)
colnames(tmp) <- c("x", "y")  
tmp <- na.omit(tmp[, c("x", "y")])

if (nrow(tmp) < 3) {
  rho <- NA
  pval <- NA} 
else {
  test <- cor.test(tmp$x,tmp$y,method = "spearman",exact = FALSE)
  rho <- unname(test$estimate)
  pval <- test$p.value}

correlations <- rbind(correlations,
                      data.frame(strain = species_short[i],comparison = "APF_mouse", rho = rho, p = pval))
rm(tmp, rho, pval, test)

tmp <- select(table, mean_rank_AF, mean_rank_APF)
colnames(tmp) <- c("x", "y")  
tmp <- na.omit(tmp[, c("x", "y")])

if (nrow(tmp) < 3) {
  rho <- NA
  pval <- NA} else {
  test <- cor.test(tmp$x,tmp$y,method = "spearman",exact = FALSE)
  rho <- unname(test$estimate)
  pval <- test$p.value}

correlations <- rbind(correlations,
                      data.frame(strain = species_short[i],comparison = "AF_APF", rho = rho, p = pval))
rm(tmp, rho, pval, test)

vitro_vivo_list_ranks[[species_short[i]]] <- table

filename <- paste(species_short[[i]], "ranks", "csv", sep = ".")
write.csv(table, filename)
rm(table, filename, axis_lim, p_AF_mouse, p_AF_APF, p_APF_mouse)
}








###for single vs. community datasets

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

single_community_list_ranks <- list()

for (i in 1:12){
  
  table <- single_vs_community %>% 
    filter(strain == species_short[[i]]) %>% 
    select(`locus tag_old`,genus, species, strain, `refseq annotation`, koeggNOG,
           `community(batch)_1`, `community(batch)_2`, `community(batch)_3`, starts_with(species_short[[i]]))
  #if there is only one value for a protein, set this to NA also  
  cols <- c(7, 8, 9)  # community
  table <- table %>% mutate(across(all_of(cols),
                                   ~ replace(.,rowSums(!is.na(pick(all_of(cols)))) < 2,NA)))
  cols <- c(10, 11, 12)  # single
  table <- table %>% mutate(across(all_of(cols),
                                   ~ replace(.,rowSums(!is.na(pick(all_of(cols)))) < 2,NA)))
  
  rm(cols)
  
  table <- table %>%
    mutate(rank_community_1 = rank(-.[[7]], na.last = "keep")) %>%
    mutate(rank_community_2 = rank(-.[[8]], na.last = "keep")) %>%
    mutate(rank_community_3 = rank(-.[[9]], na.last = "keep")) %>%
    mutate(rank_single_1 = rank(-.[[10]], na.last = "keep")) %>%
    mutate(rank_single_2 = rank(-.[[11]], na.last = "keep")) %>%
    mutate(rank_single_3 = rank(-.[[12]], na.last = "keep")) %>%
    mutate(mean_rank_community = rowMeans(select(., rank_community_1, rank_community_2, rank_community_3), na.rm = TRUE)) %>%
    mutate(mean_rank_single = rowMeans(select(., rank_single_1, rank_single_2, rank_single_3), na.rm = TRUE))%>%
    mutate(across(starts_with("mean_"), ~replace_na(., NA)))
  
  axis_lim = max(c(table$mean_rank_single, table$mean_rank_community), na.rm = TRUE)*1.1
  p<- ggplot(table, aes(x = mean_rank_single, y = mean_rank_community)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
    scale_y_continuous(limits = c(0, axis_lim), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, axis_lim), expand = c(0, 0)) +
    labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
         x = "Mean Rank in monoculture",
         y = "Mean Rank in community (AF)") +
    theme_minimal()
  
  ggsave(  filename = paste0(paste(species_short[i], "single", "community", sep = "_"), ".png"),
           plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300)
  ggsave(  filename = paste0(paste(species_short[i], "single", "community", sep = "_"), ".pdf"),
           plot = p,  bg = "white",  width = 8,  height = 6,  units = "in")
  
  
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
  
  single_community_list_ranks[[species_short[i]]] <- table
  
  filename <- paste(species_short[[i]], "ranks", "csv", sep = ".")
  write.csv(table, filename)
  rm(table, filename, axis_lim, p)
}


###plot correlations
correlations$strain <- factor(correlations$strain, levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46","I48", "I49", "YL58"))
correlations$comparison <- factor(correlations$comparison, levels = c("single_community", "AF_APF","AF_mouse", "APF_mouse"))

ggplot(filter(correlations, (comparison == "single_community" | comparison == "AF_mouse") & !(strain %in% c("YL2", "KB18", "I49"))), aes(x = comparison, y = rho, fill = strain))+
  geom_col()+
  scale_fill_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc", "#871900", "#365934","#ff9b05","#0b6760"))+
  labs(title = "Spearman Correlation Coefficients for Different Comparisons",
       x = "Comparison",
       y = "Spearman's rho") +
  scale_y_continuous(limits = c(-0.5, 1), expand = c(0, 0)) +
  facet_wrap(~strain)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(correlations, aes(x = comparison, y = rho, fill = strain))+
  geom_col()+
  scale_fill_manual(values = c("#0e8b18","#0000cc","#97b497", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc", "#871900", "#365934","#ff9b05", "#004000","#0b6760"))+
  labs(title = "Spearman Correlation Coefficients for Different Comparisons",
       x = "Comparison",
       y = "Spearman's rho") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(-0.5, 1), expand = c(0, 0)) +
  facet_wrap(~strain)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(filter(correlations, !(strain %in% c("YL2", "KB18", "I49"))), aes(x = comparison, y = rho, fill = strain))+
  geom_col()+
  scale_fill_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc", "#871900", "#365934","#ff9b05","#0b6760"))+
  labs(title = "Spearman Correlation Coefficients for Different Comparisons",
       x = "Comparison",
       y = "Spearman's rho") +
  scale_y_continuous(limits = c(-0.1, 1), expand = c(0, 0)) +
  facet_wrap(~strain)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#plot correlations averaged over species for single-AF and AF-mouse comparisons

library(ggpubr)

ggplot(filter(correlations, !strain %in% c("YL2", "KB18", "I49") & comparison %in% c("single_community", "AF_mouse")), aes(x = comparison, y = rho))+
  geom_boxplot()+ geom_jitter(width = 0.2, aes(color = strain))+
  labs(title = "Spearman Correlation Coefficients across all species for comparisons single_community and AF_mouse",
       x = "",
       y = "Spearman's rho") +
  stat_compare_means(method = "wilcox.test",label = "p.signif") +
  scale_y_continuous(limits = c(-0.3, 1.5)) +
  scale_colour_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc", "#871900", "#365934","#ff9b05","#0b6760"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

