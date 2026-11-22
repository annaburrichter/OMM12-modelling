library(openxlsx)
library(tidyverse)
library(readxl)

## data processing for mixed monocultures

mix_vs_community <- read_excel("Result_2026-06-09 09_58_01_P200_09_MaxQuant.xlsx", 
                                  sheet = "log10_prot_int_LFQ",
                                  col_types = c("skip", "skip", "skip", "skip", "text", "skip", "skip", "skip", "skip", 
                                                "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))%>%
  separate_wider_delim(1, delim = "::", names = c("locus_tag", NA, "locus_tag_new", "WP", "refseq annotation", NA, "koeggNOG"), too_many = "drop", too_few = "align_start")%>%
  left_join(annotations)
  
species_short <- c("KB1", "YL2", "KB18","YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58")


# measured protein numbers ------------------------------------------------

#original data
common_cols <- c("community(batch)_1", "community(batch)_2", "community(batch)_3")
protein_numbers_original <- map_dfr(unique(single_vs_community$strain), function(strain_) {
  strain_cols <- names(single_vs_community)[startsWith(names(single_vs_community), strain_)]
  counts <- single_vs_community %>%filter(strain == strain_) %>%  
    summarise(across(all_of(common_cols), ~sum(!is.na(.x))),across(all_of(strain_cols), ~sum(!is.na(.x))))
  names(counts) <- c(common_cols,paste0("single_", seq_along(strain_cols)))
  counts %>%mutate(strain = strain_, .before = 1)})
protein_numbers_original <- select(protein_numbers_original, strain:single_3)%>%
  pivot_longer(cols = !strain, names_to = c("condition", "replicate"), names_sep = "_", values_to = "count")
summary_numbers_original <- protein_numbers_original %>%
  group_by(strain, condition) %>%
  summarise(mean_count = mean(count, na.rm = TRUE),sd_count = sd(count, na.rm = TRUE),.groups = "drop")
  


#mixed monocultures

common_cols <- c("OMM12_1", "OMM12_2", "OMM12_3", "Mix_1", "Mix_2", "Mix_3")
protein_numbers_mix <- map_dfr(unique(mix_vs_community$org), function(strain_) {
  counts <- mix_vs_community %>%filter(org == strain_) %>%  
    summarise(across(all_of(common_cols), ~sum(!is.na(.x))))
  names(counts) <- c(common_cols)
  counts %>%mutate(strain = strain_, .before = 1)})
protein_numbers_mix <- filter(protein_numbers_mix, !is.na(strain))%>%
  pivot_longer(cols = !strain, names_to = c("condition", "replicate"), names_sep = "_", values_to = "count")
summary_numbers_mix <- protein_numbers_mix %>%
  group_by(strain, condition) %>%
  summarise(mean_count = mean(count, na.rm = TRUE),sd_count = sd(count, na.rm = TRUE),.groups = "drop")


rm(common_cols)

#statistics for mixed monocultures
ttest_results <- protein_numbers_mix %>%
  group_by(strain) %>%
  summarise(p.value = tryCatch(t.test(count ~ condition)$p.value,error = function(e) NA_real_),.groups = "drop")
label_pos <- summary_numbers_mix %>%
  group_by(strain) %>%
  summarise(xpos = max(mean_count + sd_count) * 1.05,.groups = "drop") %>%
  left_join(ttest_results, by = "strain")

#plot as bar graphs

#original dataset
ggplot(summary_numbers_original, aes(x = factor(strain, levels = c("YL58", "I49", "I48", "I46", "YL45", "YL44","YL32" , "YL31", "YL27", "KB18","YL2" ,"KB1")), 
                                     y = mean_count, 
                                    fill = factor(strain, levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32","YL44" , "YL45", "I46", "I48","I49" ,"YL58")),
                                    alpha = condition)) +
  scale_alpha_manual(values = c(1, 0.5))+
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() + # Flip coordinates for horizontal bars
  geom_errorbar(aes(ymin = mean_count - sd_count,ymax = mean_count + sd_count), width = 0.2,position = position_dodge(width = 0.9))+
  labs(
    title = "original dataset monocultures vs. community",
    x = "species",
    y = "number of detected proteins")+
  theme_bw(base_size =14) + theme(legend.title=element_blank(), 
                                  plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 100, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"),
                                  plot.title = element_text(vjust = 10, size = 14))+
  scale_fill_manual(values = colors)

ggsave(filename = "number of detected proteins original dataset.png", width = 8, height = 6, units = "in")
ggsave(filename = "number of detected proteins original dataset.pdf", width = 8, height = 6, units = "in")


#mixed dataset
ggplot(summary_numbers_mix, aes(x = factor(strain, levels = c("YL58", "I49", "I48", "I46", "YL45", "YL44","YL32" , "YL31", "YL27", "KB18","YL2" ,"KB1")), 
                                     y = mean_count, 
                                     fill = factor(strain, levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32","YL44" , "YL45", "I46", "I48","I49" ,"YL58")),
                                     alpha = condition)) +
  scale_alpha_manual(values = c(1, 0.5))+
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() + # Flip coordinates for horizontal bars
  geom_errorbar(aes(ymin = mean_count - sd_count,ymax = mean_count + sd_count), width = 0.2,position = position_dodge(width = 0.9))+
  geom_text(data = label_pos,aes(x = strain, y = xpos, label = sig),inherit.aes = FALSE, size = 5)+ 
  labs(
    title = "new dataset mixed monocultures vs. community",
    x = "species",
    y = "number of detected proteins")+
  theme_bw(base_size =14) + theme(legend.title=element_blank(), 
                                  plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 100, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"),
                                  plot.title = element_text(vjust = 10, size = 14))+
  scale_fill_manual(values = colors)

ggsave(filename = "number of detected proteins new(mix) dataset_stat.png", width = 8, height = 6, units = "in")
ggsave(filename = "number of detected proteins new(mix) dataset_stat.pdf", width = 8, height = 6, units = "in")

#both datasets

summary_numbers_both <- bind_rows(summary_numbers_mix, summary_numbers_original)

summary_numbers_both$condition <- str_replace_all(summary_numbers_both$condition, 
                            c("Mix" = "mixed monocultures (new dataset)", 
                               "single" = "original monocultures",
                               "OMM12" = "community (new dataset)",
                               "community\\(batch\\)" = "community (original dataset)"))


                                                  


ggplot(summary_numbers_both, aes(x = factor(strain, levels = c("YL58", "I49", "I48", "I46", "YL45", "YL44","YL32" , "YL31", "YL27", "KB18","YL2" ,"KB1")), 
                                y = mean_count, 
                                fill = factor(strain, levels = c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32","YL44" , "YL45", "I46", "I48","I49" ,"YL58")),
                                alpha = condition)) +
  scale_alpha_manual(values = c(1, 0.75, 0.5, 0.25))+
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_count - sd_count,ymax = mean_count + sd_count), width = 0.2,position = position_dodge(width = 0.9))+
  guides(fill = "none")+
labs(
    title = "detected protein numbers for all datasets",
    x = "species",
    y = "number of detected proteins")+
  theme_bw(base_size =14) + theme(legend.title=element_blank(), 
                                  plot.caption = element_text(hjust=c(0.9, 0.1), vjust = 100, size = 14), plot.margin = margin(40, 10, 0, 10, "pt"),
                                  plot.title = element_text(vjust = 10, size = 14))+
  scale_fill_manual(values = colors)

ggsave(filename = "number of detected proteins both datasets.png", width = 11, height = 4, units = "in")
ggsave(filename = "number of detected proteins both datasets.pdf", width = 11, height = 4, units = "in")


# protein overlap -------------------------------------------------

library(eulerr)
library(colorspace)

#separate plots
for (i in 1:12){
  
original_dataset_comm <- single_community_list_ranks[[i]]%>%
    filter(!is.na(mean_rank_community))
original_dataset_comm<- unique(original_dataset_comm$`locus tag_old`)

mix_dataset_comm <- mix_community_list_ranks[[i]]%>%
  filter(!is.na(mean_rank_community))
mix_dataset_comm<- unique(mix_dataset_comm$locus_tag)

original_dataset_mono <- single_community_list_ranks[[i]]%>%
  filter(!is.na(mean_rank_single))
original_dataset_mono<- unique(original_dataset_mono$`locus tag_old`)

mix_dataset_mono <- mix_community_list_ranks[[i]]%>%
  filter(!is.na(mean_rank_single))
mix_dataset_mono<- unique(mix_dataset_mono$locus_tag)


fit_comm <- euler(list(
    original = original_dataset_comm,
    new_dataset = mix_dataset_comm))
fit_mono <- euler(list(
  original = original_dataset_mono,
  new_dataset = mix_dataset_mono))
  
 p <- plot(fit_comm,
            fills = c(colors[i], lighten(colors[i], amount =0.5)),
            alpha = 1, 
            quantities = TRUE,
            main = paste(names(single_community_list_ranks)[i], "_communities"))
  
  ggsave(  filename = paste0(paste(names(single_community_list_ranks)[i], "_overlap_all_community"), ".png"),
           plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300, 
           path =  "./overlap_original_mix/all")
  ggsave(  filename = paste0(paste(names(single_community_list_ranks)[i], "_overlap_all_community"), ".pdf"),
           plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",
           path =  "./overlap_original_mix/all")
  
  
  q <- plot(fit_mono,
            fills = c(colors[i], lighten(colors[i], amount =0.5)),
            alpha = 1, 
            quantities = TRUE,
            main = paste(names(single_community_list_ranks)[i], "_monocultures"))
  
  ggsave(  filename = paste0(paste(names(single_community_list_ranks)[i], "_overlap_all_monoculture"), ".png"),
           plot = q,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300, 
           path =  "./overlap_original_mix/all")
  ggsave(  filename = paste0(paste(names(single_community_list_ranks)[i], "_overlap_all_monoculture"), ".pdf"),
           plot = q,  bg = "white",  width = 8,  height = 6,  units = "in",
           path =  "./overlap_original_mix/all") 
  
rm(original_dataset_comm, original_dataset_mono, mix_dataset_comm, mix_dataset_mono, fit_comm, fit_mono, p, q)
}

#plot all on the same scale

fit_list <- list()
for (i in 1:12){
  
  original_dataset_comm <- single_community_list_ranks[[i]]%>%
    filter(!is.na(mean_rank_community))
  original_dataset_comm<- unique(original_dataset_comm$`locus tag_old`)
  
  mix_dataset_comm <- mix_community_list_ranks[[i]]%>%
    filter(!is.na(mean_rank_community))
  mix_dataset_comm<- unique(mix_dataset_comm$locus_tag)
  
  original_dataset_mono <- single_community_list_ranks[[i]]%>%
    filter(!is.na(mean_rank_single))
  original_dataset_mono<- unique(original_dataset_mono$`locus tag_old`)
  
  mix_dataset_mono <- mix_community_list_ranks[[i]]%>%
    filter(!is.na(mean_rank_single))
  mix_dataset_mono<- unique(mix_dataset_mono$locus_tag)
  
  
  fit_comm <- euler(list(
    original = original_dataset_comm,
    new_dataset = mix_dataset_comm))
  fit_mono <- euler(list(
    original = original_dataset_mono,
    new_dataset = mix_dataset_mono))
  
fit_list[[paste0("fit_mono_", single_community_list_ranks[[i]]$strain[1])]]<- fit_mono
fit_list[[paste0("fit_comm_", single_community_list_ranks[[i]]$strain[1])]]<- fit_comm

rm(original_dataset_comm, original_dataset_mono, mix_dataset_comm, mix_dataset_mono, fit_comm, fit_mono)
}

#extract circle parameters

circle_pars <- map2_df(fit_list, names(fit_list),  ~{
    ell <- as.data.frame(.x$ellipses)
    ell$set <- rownames(ell)
    ell$comparison <- ifelse(grepl("mono", .y), "mono", "comm")
    ell$strain <- sub("fit_(mono|comm)_", "", .y)
    ell })%>%
  filter(is.finite(h),is.finite(k),is.finite(a),is.finite(b))

#calculate scaling factor
global_max_radius <- max(circle_pars$a, na.rm = TRUE)
target_radius <- 25
scale_factor <- target_radius / global_max_radius
#scale all parameters
circle_pars_scaled <- circle_pars %>%
  mutate(h = h * scale_factor,k = k * scale_factor, a = a * scale_factor, b = b * scale_factor)

#create circles and plotting coordinates
make_circle <- function(h, k, r, n = 200){
  theta <- seq(0, 2*pi, length.out = n)
  tibble(x = h + r*cos(theta),y = k + r*sin(theta))}

circle_coords <- pmap_dfr(circle_pars_scaled,
  function(h, k, a, b, phi,set, comparison, strain){
    make_circle(h, k, a) %>%
      mutate(set = set,comparison = comparison,strain = strain)})%>%
  mutate(panel = paste(strain, comparison, sep = "_"))

circle_coords$strain <- factor(circle_coords$strain,
                               levels=c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58"))

max_r <- max(circle_pars$a, na.rm = TRUE)

#extract counts
count_df <- purrr::map2_df(fit_list,names(fit_list),  ~{
    ell <- as.data.frame(.x$ellipses)
    ell$set <- rownames(ell)
    ell$comparison <- ifelse(grepl("mono", .y),"mono","comm")
    ell$strain <- sub("fit_(mono|comm)_", "",.y)
    
    ell$count <- c(.x$original.values["original"], .x$original.values["new_dataset"])
    ell })

count_df <- count_df %>%
  mutate(h = h * scale_factor,k = k * scale_factor) %>%
  filter(is.finite(h),is.finite(k))%>%
  mutate(edge_x = ifelse(set == "original",h - a,h + a),
         label_x = ifelse(set == "original",h - a - 2,h + a + 2),
         label_y = k)
count_df$strain <- factor(count_df$strain,
                               levels=c("KB1", "YL2", "KB18", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "I49", "YL58"))


#plot with faceting

ggplot(circle_coords, aes(x,y,group = interaction(set, strain, comparison),fill = strain, linetype = set)) +
  geom_polygon(aes(fill = strain, alpha = set, linetype = set), colour = "black",linewidth = 0.2)  +
  scale_fill_manual(values = colors) +
  scale_alpha_manual(values = c(original = 0.6,new_dataset = 0.4))+
  scale_linetype_manual(values = c(original = "solid",new_dataset = "dashed"))+
  geom_text(data = count_df,aes(x = label_x, y = label_y, label = count),inherit.aes = FALSE,size = 3)+
  geom_segment(data = count_df,aes(x = edge_x,y = k,xend = label_x,yend = label_y),inherit.aes = FALSE)+
  facet_grid( strain ~ comparison) +
  coord_equal(xlim = c(-max_r * 1.2, max_r * 1.2),ylim = c(-max_r * 1.2, max_r * 1.2))+
  theme_void()+ theme(legend.position = "bottom")

