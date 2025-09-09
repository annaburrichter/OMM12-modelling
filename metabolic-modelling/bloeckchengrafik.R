library(tidyverse)
library(sybil)
library(BacArena)

##AF simulation
AF_1 <- sim_OMM12_aam@simlist[[1]]
AF_2 <- sim_OMM12_aam@simlist[[2]]
AF_3 <- sim_OMM12_aam@simlist[[3]]
AF_4 <- sim_OMM12_aam@simlist[[4]]
AF_1$type <- as.factor(AF_1$type)
AF_2$type <- as.factor(AF_2$type)
AF_3$type <- as.factor(AF_3$type)
AF_4$type <- as.factor(AF_4$type)

#timepoint 0
ggplot(AF_1, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05", "#0b6760"),
                    labels = c("KB1", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("AF timepoint 0")

ggsave(filename = "gridplot_AF_tp0.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_AF_tp0.pdf", width = 6.14, height = 5.75, unit = "in") 

#timepoint 1
ggplot(AF_2, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05", "#0b6760"),
                    labels = c("KB1", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("AF timepoint 1")

ggsave(filename = "gridplot_AF_tp1.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_AF_tp1.pdf", width = 6.14, height = 5.75, unit = "in") 

#timepoint 2
ggplot(AF_3, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05", "#0b6760"),
                    labels = c("KB1", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("AF timepoint 2")

ggsave(filename = "gridplot_AF_tp2.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_AF_tp2.pdf", width = 6.14, height = 5.75, unit = "in") 

#timepoint 3
ggplot(AF_4, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05", "#0b6760"),
                    labels = c("KB1", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("AF timepoint 3")

ggsave(filename = "gridplot_AF_tp3.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_AF_tp3.pdf", width = 6.14, height = 5.75, unit = "in") 


##APF simulation
APF_1 <- sim_OMM12_m10@simlist[[1]]
APF_2 <- sim_OMM12_m10@simlist[[2]]
APF_3 <- sim_OMM12_m10@simlist[[3]]
APF_4 <- sim_OMM12_m10@simlist[[4]]
APF_1$type <- as.factor(APF_1$type)
APF_2$type <- as.factor(APF_2$type)
APF_3$type <- as.factor(APF_3$type)
APF_4$type <- as.factor(APF_4$type)

#timepoint 0
ggplot(APF_1, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18",  "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05", "#0b6760"),
                    labels = c("KB1",  "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("APF timepoint 0")

ggsave(filename = "gridplot_APF_tp0.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_APF_tp0.pdf", width = 6.14, height = 5.75, unit = "in") 

#timepoint 1
ggplot(APF_2, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18",  "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05", "#0b6760"),
                    labels = c("KB1",  "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("APF timepoint 1")

ggsave(filename = "gridplot_APF_tp1.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_APF_tp1.pdf", width = 6.14, height = 5.75, unit = "in") 

#timepoint 2
ggplot(APF_3, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18",  "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05", "#0b6760"),
                    labels = c("KB1", "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("APF timepoint 2")

ggsave(filename = "gridplot_APF_tp2.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_APF_tp2.pdf", width = 6.14, height = 5.75, unit = "in") 

#timepoint 3
ggplot(APF_4, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18",  "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05", "#0b6760"),
                    labels = c("KB1",  "YL31", "YL32", "YL44", "YL45", "I46", "I48", "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("APF timepoint 3")

ggsave(filename = "gridplot_APF_tp3.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_APF_tp3.pdf", width = 6.14, height = 5.75, unit = "in") 

##mouse simulation
mouse_1 <- sim_OMM12_mouse@simlist[[1]]
mouse_2 <- sim_OMM12_mouse@simlist[[2]]
mouse_3 <- sim_OMM12_mouse@simlist[[3]]
mouse_4 <- sim_OMM12_mouse@simlist[[4]]
mouse_1$type <- as.factor(mouse_1$type)
mouse_2$type <- as.factor(mouse_2$type)
mouse_3$type <- as.factor(mouse_3$type)
mouse_4$type <- as.factor(mouse_4$type)

#timepoint 0
ggplot(mouse_1, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05","#004000", "#0b6760"),
                    labels = c("KB1", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48","I49",  "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("mouse timepoint 0")

ggsave(filename = "gridplot_mouse_tp0.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_mouse_tp0.pdf", width = 6.14, height = 5.75, unit = "in") 

#timepoint 1
ggplot(mouse_2, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05","#004000", "#0b6760"),
                    labels = c("KB1", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48","I49",  "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("mouse timepoint 1")

ggsave(filename = "gridplot_mouse_tp1.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_mouse_tp1.pdf", width = 6.14, height = 5.75, unit = "in") 

#timepoint 2
ggplot(mouse_3, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05","#004000", "#0b6760"),
                    labels = c("KB1", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48","I49",  "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("mouse timepoint 2")

ggsave(filename = "gridplot_mouse_tp2.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_mouse_tp2.pdf", width = 6.14, height = 5.75, unit = "in") 

#timepoint 3
ggplot(mouse_4, aes(x = x, y = y, fill = type)) +
  geom_tile() + 
  scale_fill_manual(values = c("#0e8b18", "#cc6a00", "#7bc57b", "#11a196", "#6b00bc","#871900", "#365934", "#ff9b05","#004000", "#0b6760"),
                    labels = c("KB1", "YL27", "YL31", "YL32", "YL44", "YL45", "I46", "I48","I49",  "YL58"))+
  theme_bw() + 
  coord_fixed(ratio=1)+
  scale_x_continuous(breaks= seq(0.5,249.5, by=2))+ 
  scale_y_continuous(breaks= seq(0.5,249.5, by=2))+ 
  theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank())+ 
  labs(fill = "species") + ggtitle("mouse timepoint 3")

ggsave(filename = "gridplot_mouse_tp3.png", width = 6.14, height = 5.75, unit = "in") 
ggsave(filename = "gridplot_mouse_tp3.pdf", width = 6.14, height = 5.75, unit = "in") 
