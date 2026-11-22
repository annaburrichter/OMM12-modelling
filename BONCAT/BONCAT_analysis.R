# load packages
library(tidyverse)
library(ggpubr)
library(rstatix)


# comparison of positive BONCAT % -----------------------------------------

BONCAT_fractions <- read_delim("BONCAT_fractions.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
BONCAT_fractions_means <- BONCAT_fractions %>%
  group_by(across(-c(sort_replicate, `BONCAT+`, `BONCAT-`))) %>%
  summarise(
    BONCAT_plus_mean = mean(`BONCAT+`, na.rm = TRUE),
    BONCAT_minus_mean = mean(`BONCAT-`, na.rm = TRUE),
    .groups = "drop")%>%
  group_by(community, replicate) %>%
  mutate(
    control_BONCAT_plus = BONCAT_plus_mean[carbon_source == "-C"][1],
    BONCAT_plus_norm = BONCAT_plus_mean - control_BONCAT_plus) %>%  ungroup()

BONCAT_fractions_pairwise <-BONCAT_fractions_means %>%
  group_by(carbon_source) %>%
  pairwise_t_test(
    BONCAT_plus_norm ~ community,
    p.adjust.method = "BH")

BONCAT_fractions_anova <- BONCAT_fractions_means %>%
  filter(carbon_source != "-C") %>%
  group_by(carbon_source) %>%
  anova_test(BONCAT_plus_norm ~ community)

BONCAT_fractions_tukey <- BONCAT_fractions_means %>%
  filter(carbon_source != "-C") %>%
  group_by(carbon_source) %>%
  tukey_hsd(BONCAT_plus_norm ~ community)

ggplot(filter(BONCAT_fractions_means, carbon_source != "-C"), aes(x=community, y=BONCAT_plus_norm, colour = community))+
  geom_boxplot()+
  scale_color_manual(values = c("#FFA040","#FF6000","#800080"))+
  facet_wrap(~carbon_source)+
  theme_bw(base_size = 12)

# calculate BONCAT+/BONCAT total ------------------------------------------

BONCAT_means_long <- BONCAT_fractions_means %>%
  select(community,replicate, carbon_source, BONCAT_plus_mean, BONCAT_minus_mean) %>%
  pivot_longer(
    cols = c(BONCAT_plus_mean, BONCAT_minus_mean),
    names_to = "mean_type",
    values_to = "BONCAT_mean") %>%
  mutate(
    BONCAT = if_else(
      mean_type == "BONCAT_plus_mean","B+","B-")) %>%  select(-mean_type)%>%
  mutate(BONCAT_mean)

BONCAT_qPCR_fractions <- 
  qPCR_absolute_long%>%
  filter(community != "LO")%>%
  mutate(replicate = as.numeric(replicate))%>%
  left_join(BONCAT_means_long,by = c("community", "replicate", "carbon_source", "BONCAT"))%>%
  filter(!(community == "AF" & replicate == 3 & carbon_source == "Pyr"))%>%
  filter(!(community == "CC" & replicate == 3 & probe == "YL2" & carbon_source == "Pyr"))
#removed samples: AF Pyruvate replicate 3 - very high qPCR abundances and no corresponding DNA concentration 
#-> unrealistically high values because clearly assuming an average DNA amount does not match
#further removed: YL2 value from CC Pyruvate replicate 3 - unrealistically high and only present in one replicate

BONCAT_calc <- BONCAT_qPCR_fractions %>%
  mutate(
    BONCAT_fraction = BONCAT_mean / 100,
    weighted_abundance = abundance * BONCAT_fraction) %>%
  group_by(community, replicate, carbon_source, probe) %>%
  mutate(virtual_total = sum(weighted_abundance)) %>%
  ungroup() %>%
  mutate(
    Bplus_ratio = if_else(
      BONCAT == "B+",
      weighted_abundance / virtual_total,
      NA_real_))%>%
  mutate(log_ratio = log2(Bplus_ratio))%>%
  group_by(community, probe, carbon_source) %>%
  summarise(
    log_ratio = mean(log_ratio, na.rm = TRUE),
    .groups = "drop")



# plots -------------------------------------------------------------------

library(viridis)

ggplot(filter(BONCAT_calc, !probe %in% c("YL2", "KB18", "I49", "total")),
       aes(x = community, y = forcats::fct_rev(probe), fill = log_ratio)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "D")+
  labs(
    title = "enrichment BONCAT+ over BONCAT total", 
    x = "Community",
    y = "Species",
    fill = "BONCAT+/BONCAT_total ratio"
  ) +
  theme_minimal()+  theme(panel.grid = element_blank(), legend.position = "bottom", legend.direction = "horizontal")+
  facet_wrap(~carbon_source, labeller = labeller(carbon_source = c("-C" = "no carbon source","Pyr"= "pyruvate","Tre"= "trehalose")))




# relative plots of B+ fractions ------------------------------------------

BONCAT_relative <- qPCR_absolute_long%>%
  group_by(community, replicate, carbon_source, BONCAT) %>%
  mutate(
    sample_total = abundance[probe == "total"][1],
    relative_value = abundance/sample_total) %>%  ungroup()%>%
  filter(community != "LO")

BONCAT_relative_mean <- BONCAT_relative %>%
  filter(probe != "total") %>%
  group_by(community, carbon_source, BONCAT, probe) %>%
  summarise(relative_value = mean(relative_value, na.rm = TRUE),.groups = "drop")

#load all strain designations and Anna w.'s colour scheme
load("all_probes.rda")
col_annaw <- load("Colorscheme_AnnaW_ordered.rda")

ggplot(filter(BONCAT_relative_mean, probe != "total"), aes(x = interaction(community, carbon_source, BONCAT), y = relative_value, fill=factor(probe, levels=all_probes))) + 
  geom_col() + scale_fill_manual(values = col_annaw_ordered) +
  ylab("relative abundance") + 
  theme_bw()+
  theme(legend.title=element_blank(), axis.title.x =element_blank(), legend.position = "left", axis.text.x = element_text(angle=90))

ggplot(filter(BONCAT_relative_mean, probe != "total" & carbon_source == "Tre" & BONCAT =="B+"), aes(x = interaction(community, carbon_source, BONCAT), y = relative_value, fill=factor(probe, levels=all_probes))) + 
  geom_col() + scale_fill_manual(values = col_annaw_ordered) +
  ylab("relative abundance") + 
  theme_bw()+
  theme(legend.title=element_blank(), axis.title.x =element_blank(), legend.position = "left", axis.text.x = element_text(angle=90))

  

