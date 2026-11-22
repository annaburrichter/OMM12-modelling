library(tidyverse)
library(data.table)

data <- read.delim("proteomics_linear.txt")
data_filtered <- data %>% select(9:17)%>% filter((AAM_1 != 0 & AAM_2 != 0)|(AAM_1 != 0 & AAM_3 != 0)|(AAM_2 != 0 & AAM_3 != 0)|(M10_1 != 0 & M10_2 != 0)|(M10_1 != 0 & M10_3 != 0)|(M10_2 != 0 & M10_3 != 0)|(mouse_1 != 0 & mouse_2 != 0)|(mouse_1 != 0 & mouse_3 != 0)|(mouse_2 != 0 & mouse_3 != 0))
data_hc <- t(data_filtered)
row.names(data_hc) <- c("AF_1", "AF_2", "AF_3", "APF_1", "APF_2", "APF_3", "mouse_1", "mouse_2", "mouse_3")

distance_prot <- dist(data_hc, method = "euclidean")
clustering_prot <- hclust(distance_prot, method ="complete")

plot(clustering_prot, labels = NULL, hang = 0.1, check = TRUE,
     axes = TRUE, frame.plot = FALSE, ann = TRUE,
     main = "hierarchical clustering based on proteome",
     sub = NULL, xlab = "", ylab = "distance")

library(plotly)
library(ggdendro)
plot_cluster <- dendro_data(clustering_prot, type = "rectangle")

#horizontal

cp <- ggplot() +
  geom_segment(data = segment(plot_cluster), 
               aes(x = x, y = y, xend = xend, yend = yend)) +
  
  geom_text(data = label(plot_cluster), 
            aes(x = x, y = y, label = label, hjust = -0.5), 
            size = 5, colour = c("#FFA040", "#FFA040", "#FFA040", "#FF6000","#FF6000","#FF6000","#800080", "#800080", "#800080"),
                                 fontface = "bold") +
  
  coord_flip(clip = "off") +  scale_y_reverse(expand = c(0.2, 0)) + 
  theme(plot.margin = margin(t=1, r=4, b=0, l=0, unit = "in")) +
  theme_minimal(base_size = 14) + 
  labs(x ="", y = "", title ="hierarchical cluster analysis based on proteomics data")

ggsave(cp, file="HC.pdf", width=25, height=15, units="cm")
ggsave(cp, file="HC.png", width=25, height=15, units="cm", bg = "white")


#vertical

cpv <- ggplot() +
  geom_segment(data = segment(plot_cluster), 
               aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_y_continuous(limits = c(-0.5e+10, 4e+10)) + 
  
  geom_text(data = label(plot_cluster), 
            aes(x = x, y = y, label = label, hjust = 1.1, angle = 90), 
            size = 5, colour = c("#FFA040", "#FFA040", "#FFA040", "#FF6000","#FF6000","#FF6000","#800080", "#800080", "#800080"),
            fontface = "bold") +
  geom_hline(yintercept = 0, color="darkgrey") + 
  theme_classic(base_size = 14) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.line.x = element_blank())+
  labs(x ="", y = "", title ="hierarchical cluster analysis based on proteomics data")

ggsave(cpv, file="HC_v.pdf")
ggsave(cpv, file="HC_v.png", bg = "white")
