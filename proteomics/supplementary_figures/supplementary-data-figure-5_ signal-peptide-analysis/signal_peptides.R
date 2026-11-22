library(tidyverse)
library(readr)

dbcan_signalp_proteins <- read_csv("dbcan-signalp-proteins.csv")
vitro_vivo_normalised <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/in vitro vs. in vivo/enrichment normalisation by total intensity/all_normalised_annotated.csv")
I46_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/I46.direct.csv")
I48_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/I48_direct.csv")
I49_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/I49.direct.csv")
KB1_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/KB1_.direct.csv")
KB18_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/KB18.direct.csv")
YL2_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/YL2_.direct.csv")
YL27_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/YL27.direct.csv")
YL31_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/YL31.direct.csv")
YL32_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/YL32.direct.csv")
YL44_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/YL44.direct.csv")
YL45_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/YL45.direct.csv")
YL58_ranks <- read_csv("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/rank-based analysis/direct/YL58.direct.csv")

signalp <- left_join(dbcan_signalp_proteins, protein_annotation, by = c('Gene ID' = 'prodigal'))

I46_signalp_single_comm <- inner_join(I46_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_single_1, rank_single_2, rank_single_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

I48_signalp_single_comm <- inner_join(I48_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_I48_1, rank_I48_2, rank_I48_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

I49_signalp_single_comm <- inner_join(I49_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_single_1, rank_single_2, rank_single_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

KB1_signalp_single_comm <- inner_join(KB1_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_single_1, rank_single_2, rank_single_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

KB18_signalp_single_comm <- inner_join(KB18_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_single_1, rank_single_2, rank_single_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

YL2_signalp_single_comm <- inner_join(YL2_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_single_1, rank_single_2, rank_single_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

YL27_signalp_single_comm <- inner_join(YL27_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_single_1, rank_single_2, rank_single_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

YL31_signalp_single_comm <- inner_join(YL31_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_single_1, rank_single_2, rank_single_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

YL32_signalp_single_comm <- inner_join(YL32_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_single_1, rank_single_2, rank_single_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

YL44_signalp_single_comm <- inner_join(YL44_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_single_1, rank_single_2, rank_single_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

YL45_signalp_single_comm <- inner_join(YL45_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_single_1, rank_single_2, rank_single_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

YL58_signalp_single_comm <- inner_join(YL58_ranks, signalp, by = c('locus tag_old' = 'locus_tag'))%>%
  rowwise()%>%mutate(mean_rank_single = mean(c(rank_single_1, rank_single_2, rank_single_3)))%>%
  rowwise()%>%mutate(mean_rank_comm = mean(c(rank_comm_1, rank_comm_2, rank_comm_3)))%>%
  mutate(rank_diff = mean_rank_single - mean_rank_comm)

#single - community all 12
all_signalp_single_comm <- bind_rows(KB1_signalp_single_comm, YL2_signalp_single_comm, KB18_signalp_single_comm, YL27_signalp_single_comm,
                                     YL31_signalp_single_comm, YL32_signalp_single_comm, YL44_signalp_single_comm, YL45_signalp_single_comm, 
                                     I46_signalp_single_comm, I48_signalp_single_comm, I49_signalp_single_comm, YL58_signalp_single_comm)%>%
                                    mutate(mean_rank_comm = ifelse(is.na(mean_rank_comm), 1700, mean_rank_comm))%>%
                                    mutate(mean_rank_single = ifelse(is.na(mean_rank_single), 1700, mean_rank_single))%>%
                                    pivot_longer(cols = c(mean_rank_single, mean_rank_comm),
                                           names_to = "Condition",
                                           values_to = "Rank")

all_signalp_single_comm$org.x <- factor(all_signalp_single_comm$org.x, 
                                           levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", 
                                                      "YL44", "YL45", "I46", "I48", "I49", "YL58"))

all_signalp_single_comm$Condition <- factor(all_signalp_single_comm$Condition, 
                                        levels = c("mean_rank_single", "mean_rank_comm"))

ggplot(all_signalp_single_comm, aes(x = Condition, y = Rank, group = `locus tag_old`, color = org.x)) +
  geom_line(color = "lightgrey", size = 0.3) +  # thin line between points
  geom_point(aes(color = org.x))+
  theme_minimal()+
  theme(axis.title.x = element_blank(), legend.position = "none")+
  scale_y_reverse(breaks = c(0, 500, 1000, 1500, 1700), labels = c("0", "500", "1000", "1500", "not detected"))+
  geom_hline(yintercept = 1700, size =0.5, color ="grey")+
  scale_x_discrete(labels = c("mean_rank_single" = "monoculture", "mean_rank_comm" = "community"))+
  facet_wrap(~org.x)+
  scale_color_manual(values=c("#0e8b18","#0000cc", "#97b497","#cc6a00", "#7bc57b", "#11a196","#6b00bc", "#871900", "#365934",  "#ff9b05","#004000", "#0b6760"))

#single - community 8
all8_signalp_single_comm <- bind_rows(KB1_signalp_single_comm, YL27_signalp_single_comm,
                                     YL31_signalp_single_comm, YL32_signalp_single_comm, YL45_signalp_single_comm, 
                                     I46_signalp_single_comm, I48_signalp_single_comm, YL58_signalp_single_comm)%>%
  mutate(mean_rank_comm = ifelse(is.na(mean_rank_comm), 1700, mean_rank_comm))%>%
  mutate(mean_rank_single = ifelse(is.na(mean_rank_single), 1700, mean_rank_single))%>%
  pivot_longer(cols = c(mean_rank_single, mean_rank_comm),
               names_to = "Condition",
               values_to = "Rank")

all8_signalp_single_comm$org.x <- factor(all8_signalp_single_comm$org.x, 
                                        levels = c("KB1", "YL27", "YL31", "YL32", 
                                                   "YL45", "I46", "I48", "YL58"))

all8_signalp_single_comm$Condition <- factor(all8_signalp_single_comm$Condition, 
                                            levels = c("mean_rank_single", "mean_rank_comm"))

ggplot(all8_signalp_single_comm, aes(x = Condition, y = Rank, group = `locus tag_old`, color = org.x)) +
  geom_line(color = "lightgrey", size = 0.3) +  # thin line between points
  geom_point(aes(color = org.x))+
  theme_minimal()+
  theme(axis.title.x = element_blank(), legend.position = "none")+
  scale_y_reverse(breaks = c(0, 500, 1000, 1500, 1700), labels = c("0", "500", "1000", "1500", "not detected"))+
  geom_hline(yintercept = 1700, size =0.5, color ="grey")+
  scale_x_discrete(labels = c("mean_rank_single" = "monoculture", "mean_rank_comm" = "community"))+
  facet_wrap(~org.x, nrow = 2, ncol = 4)+
  scale_color_manual(values=c("#0e8b18","#cc6a00", "#7bc57b", "#11a196", "#871900", "#365934",  "#ff9b05", "#0b6760"))

ggplot(all8_signalp_single_comm, aes(x = Condition, y = Rank, group = `locus tag_old`, color = org.x)) +
  geom_line(color = "lightgrey", size = 0.3) +  # thin line between points
  geom_point(aes(color = org.x))+
  theme_minimal()+
  theme(axis.title.x = element_blank(), legend.position = "none")+
  scale_y_reverse(breaks = c(0, 500, 1000, 1500, 1700), labels = c("0", "500", "1000", "1500", "not detected"))+
  geom_hline(yintercept = 1700, size =0.5, color ="grey")+
  scale_x_discrete(labels = c("mean_rank_single" = "monoculture", "mean_rank_comm" = "community"))+
  facet_wrap(~org.x, nrow = 4, ncol = 2)+
  scale_color_manual(values=c("#0e8b18","#cc6a00", "#7bc57b", "#11a196", "#871900", "#365934",  "#ff9b05", "#0b6760"))


#vitro - vivo
signalp_vitro_vivo <- inner_join(vitro_vivo_normalised, signalp)%>%
  rowwise()%>%mutate(AF = mean(c(AAM_1, AAM_2, AAM_3), na.rm = TRUE))%>%
  rowwise()%>%mutate(APF = mean(c(M10_1, M10_2, M10_3), na.rm = TRUE))%>%
  rowwise()%>%mutate(mouse = mean(c(mouse_1, mouse_2, mouse_3), na.rm = TRUE))%>%
  pivot_longer(cols = c(AF:mouse), names_to = "Condition", values_to = "Value")%>%
  mutate(Value = ifelse(is.na(Value), -0.004, Value))
  

ggplot(signalp_vitro_vivo, aes(x = Condition, y = Value, group = locus_tag, color = Condition)) +
  geom_line(color = "lightgrey", size = 0.3) +  # thin line between points
  geom_point(aes(color = Condition))+
  theme_minimal()+
  scale_y_continuous(breaks = c(-0.004,0, 0.005, 0.010, 0.015), labels = c("not detected","0.000", "0.005", "0.010", "0.015"))+
  geom_hline(yintercept = -0.004, size =0.5, color ="grey")+
  theme(axis.title.x = element_blank(), legend.position = "none")+
  facet_wrap(~species)+
  scale_color_manual(values=c("#FFA040","#FF6000", "#800080"))


write.csv(all_signalp_single_comm, file = "signalpeptides_single_comm.csv")
write.csv(signalp_vitro_vivo, file = "signalpeptides_vitro_vivo.csv")
save.image("C:/Users/go82jux/OneDrive - TUM/metabolic models AG Kaleta/proteomics data/signal peptide analysis/20250828.RData")