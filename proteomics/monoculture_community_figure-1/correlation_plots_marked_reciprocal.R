library(openxlsx)
library(tidyverse)
library(readxl)



# interactive plots for marking choice ------------------------------------
#plot correlation, make interactive plot for labelling choice (unclear if opening plotly object works in a loop)
library(plotly)
for (i in 1:12){
axis_lim = max(c(single_community_list_ranks[[i]]$mean_rank_single, single_community_list_ranks[[i]]$mean_rank_community), na.rm = TRUE)*1.1
p<- ggplot(single_community_list_ranks[[i]], aes(x = mean_rank_single, y = mean_rank_community, text = koeggNOG)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  scale_y_continuous(limits = c(-1, axis_lim), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-1, axis_lim), expand = c(0, 0)) +
  labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
       x = "Mean Rank in monoculture",
       y = "Mean Rank in community (AF)") +
  theme_minimal()

ggplotly(p, tooltip = "text")
}


# unmarked coloured plots -------------------------------------------------

#plot and save correlations in respective color for monocultures vs. community 

for (i in 1:12){

p<- ggplot(single_community_list_ranks[[i]], aes(x = 1/mean_rank_single, y = 1/mean_rank_community)) +
  geom_point(colour= colors[i]) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(limits = c(1e-3, 1.1)) +
  scale_x_log10(limits = c(1e-3, 1.1)) +
  labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
       x = "1/Mean Rank in monoculture",
       y = "1/Mean Rank in community (AF)") +
  theme_minimal()+
  theme(legend.position = "none")
  
 
ggsave(  filename = paste0(paste(species_short[i], "single", "community", "reciprocal", sep = "_"), ".png"),
         plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300, 
         path =  "./correlation plots single vs. community marked reciprocal")
ggsave(  filename = paste0(paste(species_short[i], "single", "community", "reciprocal", sep = "_"), ".pdf"),
         plot = p,  width = 8,  height = 6,  units = "in",
         path =  "./correlation plots single vs. community marked reciprocal")

rm(p)
}

# marked plots ------------------------------------------------------------
##create marked plots for YL27, YL31, YL32, I46, I48, YL58
#marking color options: "#d81b60"(magenta), "#0072ff" (blue, for YL27 and I48)
# YL27 --------------------------------------------------------------------
i<- 4

#single vs. community
table <- single_community_list_ranks[[i]]%>%
  mutate(mark = str_detect(koeggNOG, regex("SusC|SuSD|amylase|chitinase|osaminidase", ignore_case = TRUE)))

p<- ggplot(table, aes(x = 1/mean_rank_single, y = 1/mean_rank_community, colour = mark)) +
  geom_point() +
  scale_colour_manual(values = c("FALSE" = colors[i], "TRUE" = "#d81b60"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(limits = c(1e-3, 1.1)) +
  scale_x_log10(limits = c(1e-3, 1.1)) +
  labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
       x = "1/Mean Rank in monoculture",
       y = "1/Mean Rank in community (AF)") +
  theme_minimal()+
  theme(legend.position = "none")

ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".png"),
         plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300, 
         path =  "./correlation plots single vs. community marked reciprocal")
ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".pdf"),
         plot = p,  width = 8,  height = 6,  units = "in",
         path =  "./correlation plots single vs. community marked reciprocal")
export_table <- filter(table, mark == TRUE)
filename <- paste("./correlation plots single vs. community marked reciprocal/", species_short[[i]],"_single vs. comm", "_marked", ".csv", sep = "")
write.csv(export_table, file = filename, )

rm(p, table, filename, export_table)


# YL31 --------------------------------------------------------------------
i<- 5
#filtering put propionate metabolism: 
#I5Q84_01935|I5Q84_01940|I5Q84_01950|I5Q84_01955|I5Q84_01970|I5Q84_01980

#single vs. community
table <- single_community_list_ranks[[i]]%>%
  mutate(mark = str_detect(`locus tag_old`, regex("I5Q84_01935|I5Q84_01940|I5Q84_01950|I5Q84_01955|I5Q84_01970|I5Q84_01980", ignore_case = TRUE)))

p<- ggplot(table, aes(x = 1/mean_rank_single, y = 1/mean_rank_community, colour = mark)) +
  geom_point() +
  scale_colour_manual(values = c("FALSE" = colors[i], "TRUE" = "#d81b60"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(limits = c(1e-3, 1.1)) +
  scale_x_log10(limits = c(1e-3, 1.1)) +
  labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
       x = "1/Mean Rank in monoculture",
       y = "1/Mean Rank in community (AF)") +
  theme_minimal()+
  theme(legend.position = "none")

ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".png"),
         plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300, 
         path =  "./correlation plots single vs. community marked reciprocal")
ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".pdf"),
         plot = p,  width = 8,  height = 6,  units = "in",
         path =  "./correlation plots single vs. community marked reciprocal")
export_table <- filter(table, mark == TRUE)
filename <- paste("./correlation plots single vs. community marked reciprocal/", species_short[[i]],"_single vs. comm", "_marked", ".csv", sep = "")
write.csv(export_table, file = filename, )

rm(p, table, filename, export_table)


# YL32 --------------------------------------------------------------------
i<- 6
#filtering out butyrate metabolism: 
#I5Q83_03225|I5Q83_03230|I5Q83_03240|I5Q83_03295|I5Q83_03300|I5Q83_03305|I5Q83_03310|I5Q83_03315

#single vs. community
table <- single_community_list_ranks[[i]]%>%
  mutate(mark = str_detect(`locus tag_old`, regex("I5Q83_03225|I5Q83_03230|I5Q83_03240|I5Q83_03295|I5Q83_03300|I5Q83_03305|I5Q83_03310|I5Q83_03315", ignore_case = TRUE)))

p<- ggplot(table, aes(x = 1/mean_rank_single, y = 1/mean_rank_community, colour = mark)) +
  geom_point() +
  scale_colour_manual(values = c("FALSE" = colors[i], "TRUE" = "#d81b60"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(limits = c(1e-3, 1.1)) +
  scale_x_log10(limits = c(1e-3, 1.1)) +
  labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
       x = "1/Mean Rank in monoculture",
       y = "1/Mean Rank in community (AF)") +
  theme_minimal()+
  theme(legend.position = "none")

ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".png"),
         plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300, 
         path =  "./correlation plots single vs. community marked reciprocal")
ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".pdf"),
         plot = p,  width = 8,  height = 6,  units = "in",
         path =  "./correlation plots single vs. community marked reciprocal")
export_table <- filter(table, mark == TRUE)
filename <- paste("./correlation plots single vs. community marked reciprocal/", species_short[[i]],"_single vs. comm", "_marked", ".csv", sep = "")
write.csv(export_table, file = filename, )

rm(p, table, filename, export_table)


# I46 --------------------------------------------------------------------
i<- 9

#single vs. community
table <- single_community_list_ranks[[i]]%>%
  mutate(mark = str_detect(koeggNOG, regex("PTS", ignore_case = FALSE)))

p<- ggplot(table, aes(x = 1/mean_rank_single, y = 1/mean_rank_community, colour = mark)) +
  geom_point() +
  scale_colour_manual(values = c("FALSE" = colors[i], "TRUE" = "#d81b60"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(limits = c(1e-3, 1.1)) +
  scale_x_log10(limits = c(1e-3, 1.1)) +
  labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
       x = "1/Mean Rank in monoculture",
       y = "1/Mean Rank in community (AF)") +
  theme_minimal()+
  theme(legend.position = "none")

ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".png"),
         plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300, 
         path =  "./correlation plots single vs. community marked reciprocal")
ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".pdf"),
         plot = p,  width = 8,  height = 6,  units = "in",
         path =  "./correlation plots single vs. community marked reciprocal")
export_table <- filter(table, mark == TRUE)
filename <- paste("./correlation plots single vs. community marked reciprocal/", species_short[[i]],"_single vs. comm", "_marked", ".csv", sep = "")
write.csv(export_table, file = filename, )

rm(p, table, filename, export_table)

# I48 --------------------------------------------------------------------
i<- 10

#single vs. community
table <- single_community_list_ranks[[i]]%>%
  mutate(mark = str_detect(koeggNOG, regex("SusC|SuSD|mannosidase", ignore_case = TRUE)))

p<- ggplot(table, aes(x = 1/mean_rank_single, y = 1/mean_rank_community, colour = mark)) +
  geom_point() +
  scale_colour_manual(values = c("FALSE" = colors[i], "TRUE" = "#d81b60"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(limits = c(1e-3, 1.1)) +
  scale_x_log10(limits = c(1e-3, 1.1)) +
  labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
       x = "1/Mean Rank in monoculture",
       y = "1/Mean Rank in community (AF)") +
  theme_minimal()+
  theme(legend.position = "none")

ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".png"),
         plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300, 
         path =  "./correlation plots single vs. community marked reciprocal")
ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".pdf"),
         plot = p,  width = 8,  height = 6,  units = "in",
         path =  "./correlation plots single vs. community marked reciprocal")
export_table <- filter(table, mark == TRUE)
filename <- paste("./correlation plots single vs. community marked reciprocal/", species_short[[i]],"_single vs. comm", "_marked", ".csv", sep = "")
write.csv(export_table, file = filename, )

rm(p, table, filename, export_table)

# YL58 --------------------------------------------------------------------
i<- 12
#to filter out: WLP: I5Q86_21630|I5Q86_21635|I5Q86_21640|I5Q86_03035|I5Q86_06130|I5Q86_15835|I5Q86_20645


#single vs. community
table <- single_community_list_ranks[[i]]%>%
  mutate(mark = str_detect(`locus tag_old`, regex("I5Q86_21630|I5Q86_21635|I5Q86_21640|I5Q86_03035|I5Q86_06130|I5Q86_15835|I5Q86_20645", ignore_case = TRUE)))

p<- ggplot(table, aes(x = 1/mean_rank_single, y = 1/mean_rank_community, colour = mark)) +
  geom_point() +
  scale_colour_manual(values = c("FALSE" = colors[i], "TRUE" = "#d81b60"))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")+
  scale_y_log10(limits = c(1e-3, 1.1)) +
  scale_x_log10(limits = c(1e-3, 1.1)) +
  labs(title = paste("Mean Rank Comparison for", species_short[[i]]),
       x = "1/Mean Rank in monoculture",
       y = "1/Mean Rank in community (AF)") +
  theme_minimal()+
  theme(legend.position = "none")

ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".png"),
         plot = p,  bg = "white",  width = 8,  height = 6,  units = "in",  dpi = 300, 
         path =  "./correlation plots single vs. community marked reciprocal")
ggsave(  filename = paste0(paste(species_short[i], "single", "community","marked", "reciprocal", sep = "_"), ".pdf"),
         plot = p,  width = 8,  height = 6,  units = "in",
         path =  "./correlation plots single vs. community marked reciprocal")
export_table <- filter(table, mark == TRUE)
filename <- paste("./correlation plots single vs. community marked reciprocal/", species_short[[i]],"_single vs. comm", "_marked", ".csv", sep = "")
write.csv(export_table, file = filename, )

rm(p, table, filename, export_table)

