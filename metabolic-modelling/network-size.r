library(sybil)
library(ggplot2)
library(stringr)
library(ggpubr)
library(rstatix)


models <- readRDS("dat/OMM_munich-csm-aam_20241205.RDS")

rxn_num <- sapply(models, function(m) m@react_num)
act_num <- sapply(models, function(x) {
  sol=optimizeProb(x, retOptSol=T, algorithm="mtf")
  sum(getFluxDist(sol)!=0)
})


rxn_num_dt <- data.table(id=names(models), rxn_num, act_num)
rxn_num_dt[, org := str_extract(id, "[A-Z]+[0-9]+(?=-)")]
rxn_num_dt[, condition:=str_extract(id, "(?<=-).*(?=_)")]
rxn_num_dt[org==condition, condition:="mono"]
rxn_num_dt[, model_type:=str_extract(id, "(all|strict)")]
rxn_num_sel_dt <- rxn_num_dt[model_type=="all" & condition!="batch"]
rxn_num_sel_dt[, condition:=factor(condition, levels = c("mono", "AAM", "M10", "mouse"))]

stat_test <- add_xy_position(wilcox_test(rxn_num_sel_dt, rxn_num ~ condition))
ggplot(rxn_num_sel_dt, aes(x=condition, y=rxn_num)) + geom_boxplot() + theme_classic(base_size=14) + xlab("") + ylab("Network size (# reactions)") +
  stat_pvalue_manual(stat_test, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)
ggsave("img/network-size_rxn-num.pdf", width=4, height=3.5)
ggsave("img/network-size_rxn-num.png", width=4, height=3.5)

stat_test2 <- add_xy_position(wilcox_test(rxn_num_sel_dt, act_num ~ condition))
ggplot(rxn_num_sel_dt, aes(x=condition, y=act_num)) + 
  geom_boxplot() + 
  theme_classic(base_size=14) + xlab("") + ylab("Network size (active reactions)") +
  stat_pvalue_manual(stat_test2, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)
ggsave("img/network-size_active-rxn-num.pdf", width=4, height=3.5)
ggsave("img/network-size_active-rxn-num.png", width=4, height=3.5)
