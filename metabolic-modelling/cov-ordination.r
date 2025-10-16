library(logisticPCA)
library(ggplot2)
library(FactoMineR)
library(ggrepel)
library(ggpubr)

stat.dt <- readRDS("dat/core-rxn.RDS")[[1]]
rxn.dt <- readRDS("dat/core-rxn.RDS")[[2]]

org.colors <- fread("/home/jo/uni/OMM_munich/dat/OMM12-colour-code.txt")
org.colors[, org := str_extract(species, "[A-Z]+[0-9]+$")]

rm.org <- c("YL2", "I49", "KB18") # very low abundant species in proteome
rm.cond <- c("fermenter", "batch") # duplicate conditions

colpal2 <- c("#C48736", "#CE54D1", "#96CED5", "#76D73C", "#403552", "#D4477D", "#5A7E36", "#D19EC4", "#CBC594", "#722A2D", "#D0CD47", "#CF4A31", "#7B6FD0", "#597873", "#6CD3A7", "#484125", "#C17E73", "#688EC1", "#844081", "#7DD06F")


# coverage proteome and model genes
stat.dt[, cond := col]
stat.dt[cond == id, cond := "mono"]
stat.melt <- melt(stat.dt[!is.na(cond) & !(id %in% rm.org & cond != "mono")], id.vars = c("id", "col", "cond"))
stat.melt[, cond2 := ifelse(cond == "AAM", "AF", ifelse(cond == "M10", "APF", cond))]
stat.melt[, cond2 := factor(cond2, levels = c("mono", "batch", "AF", "APF", "mouse"), ordered = T)]
stat.melt[, variable2 := ifelse(variable == "mod.cov", "model genes", "proteome")]
ggplot(stat.melt[cond != "batch"], aes(x = variable2, y = value)) +
    geom_boxplot() +
    facet_wrap(~cond2, nrow = 1) +
    theme_bw(base_size = 14) +
    xlab("") +
    ylab("Coverage between proteome\n and model genes") +
    ylim(0, 1) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.title.y = element_text(size=10))
ggsave("img/mod-proteome-coverage.pdf", height = 2.5, width = 4.5)
ggsave("img/mod-proteome-coverage.png", height = 2.5, width = 4.5)


# jaccard distance between conditions
library(combinat)
rxn_dist_dt <- rxn.dt[!col %in% rm.cond & !id %in% rm.org]
rxn_dist_dt[, condition := str_extract(col, "batch|fermenter|AAM|M10|mouse")]
rxn_dist_dt[is.na(condition), condition := "mono"]
rxn_dist_dt[condition == "AAM", condition := "AF"]
rxn_dist_dt[condition == "M10", condition := "APF"]

comb <- combn(unique(rxn_dist_dt$condition), 2)
jaccard_dt <- data.table()
for (o in unique(rxn_dist_dt$id)) {
    for (i in 1:ncol(comb)) {
        c <- comb[, i]
        rxn1 <- rxn_dist_dt[condition == c[1] & id == o, V3]
        rxn2 <- rxn_dist_dt[condition == c[2] & id == o, V3]
        jaccard <- length(intersect(rxn1, rxn2)) / length(union(rxn1, rxn2))
        jaccard_dt <- rbind(jaccard_dt, data.table(org = o, cond1 = c[1], cond2 = c[2], jaccard))
    }
}
jaccard_dt[, cond := paste0(cond1, " vs. ", cond2)]
ggplot(jaccard_dt, aes(x = cond, y = 1-jaccard)) +
    geom_boxplot() +
    ylab("Dissimilarity (Jaccard)") +
    xlab("") +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1))
ggsave("img/core-rxn_similarity.pdf", width = 4, height = 3)
ggsave("img/core-rxn_similarity.png", width = 4, height = 3)


# jaccard distance within conditions
comb2 <- combn(unique(rxn_dist_dt$id), 2)
jaccard_dt2 <- data.table()
for (c in unique(rxn_dist_dt$condition)) {
    for (i in 1:ncol(comb2)) {
        o <- comb2[, i]
        rxn1 <- rxn_dist_dt[condition == c & id == o[1], V3]
        rxn2 <- rxn_dist_dt[condition == c & id == o[2], V3]
        jaccard2 <- length(intersect(rxn1, rxn2)) / length(union(rxn1, rxn2))
        jaccard_dt2 <- rbind(jaccard_dt2, data.table(org1 = o[1], org2 = o[2], cond = c, jaccard2))
    }
}
jaccard_dt2[, cmp := paste0(org1, " vs. ", org2)]
jaccard_dt2[,cond2:=factor(ifelse(cond=="AF", "vitro (AF)", cond), levels=c("mono", "vitro (AF)", "mouse"), ordered=T)]
ggplot(jaccard_dt2[cond!="APF"], aes(x = cond2, y = 1-jaccard2)) +
    geom_boxplot() +
    ylab("Dissimilarity (Jaccard)") + scale_y_continuous(expand = expansion(mult = c(0, 0.1)), breaks=c(0,0.25,0.5,0.75,1)) + 
    xlab("") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) +
    stat_compare_means(label = "p.signif", method = "wilcox", comparisons = as.list(as.data.frame(combinat::combn(1:3, 2))))
compare_means(jaccard2 ~ cond, jaccard_dt2)
jaccard_dt2[,list(median=median(jaccard2), mean=mean(jaccard2)), by=cond]
ggsave("img/core-rxn_similarity2.pdf", width = 2, height = 3)
ggsave("img/core-rxn_similarity2.png", width = 2, height = 3)


# genomic coverage bvy metabolic model
models <- readRDS("dat/OMM_munich-aam_20241205.RDS")
mod_genes <- sapply(models, function(m) length(unique(m@genes)))
genome_files <- list.files("dat/", pattern = "faa.gz", full.names = T)
genome_genes <- sapply(genome_files, function(f) {
    length(unique(Biostrings::readAAStringSet(f)))
})
summary(mod_genes / genome_genes)
genomic_cov_dt <- data.table(org=names(mod_genes), mod_genes,genome_genes)
ggplot(genomic_cov_dt, aes(y=mod_genes/genome_genes)) + scale_x_discrete( ) +
    geom_boxplot() + theme_bw(base_size = 14) + 
  xlab("") + ylab("Coverage of genomic genes\n in OMM12 in models (%)")
ggsave("img/genomic_coverage.pdf", height = 3, width = 2.5)
ggsave("img/genomic_coverage.png", height = 3, width = 2.5)


# logPCA
rxn.wide <- dcast(rxn.dt[!id %in% rm.org & !col %in% rm.cond], formula = id + col ~ V3, fun.aggregate = length, value.var = "V3")
rxn.mat <- as.matrix(rxn.wide[, -c(1, 2)])
rxn.mat <- rxn.mat[, apply(rxn.mat, 2, function(x) {
  (length(unique(x)) > 1)
})] # remove constant rows
rxn.mat <- apply(rxn.mat, 2, as.logical)
rownames(rxn.mat) <- paste0(rxn.wide$id, "-", rxn.wide$col)
logpca <- logisticPCA(rxn.mat)
logpca$prop_deviance_expl
logpca.dt <- data.table(dim1 = logpca$PCs[, 1], dim2 = logpca$PCs[, 2], id = rownames(logpca$PCs))
logpca.dt[, org := str_extract(id, ".*?(?=-)")]
logpca.dt[, condition := str_extract(id, "batch|fermenter|AAM|M10|mouse")]
logpca.dt[is.na(condition), condition := "mono"]
logpca.dt <- merge(logpca.dt, org.colors, by = "org")
logpca.dt[condition == "AAM", condition := "AF"]
logpca.dt[condition == "M10", condition := "APF"]
ggplot(logpca.dt, aes(x = dim1, y = dim2)) +
    geom_point(size = 4, aes(color = `Hex code`, fill = `Hex code`, shape = condition)) +
    xlab("Dim 1") +
    ylab("Dim 2") +
    theme_bw(base_size = 14) +
    scale_fill_identity() +
    scale_color_identity() +
    facet_wrap(~org) +
    geom_label_repel(aes(label = condition)) +
    scale_shape_manual(values = c(21:24)) +
    labs(caption = paste0("Explained deviance: ", round(100 * logpca$prop_deviance_expl, 1), "%")) +
    theme(plot.caption = element_text(hjust = 0.5, size = rel(1.2)))
ggsave("img/logpca-core_rxn_all2.pdf", height = 8, width = 10)
ggsave("img/logpca-core_rxn_all2.png", height = 8, width = 10)

dist_mat <- as.matrix(dist(logpca.dt[,.(dim1,dim2)]))
rownames(dist_mat) <- logpca.dt$id; colnames(dist_mat) <- rownames(dist_mat)
dist_melt <- data.table(reshape2::melt(dist_mat, varnames=c("id1","id2"), value.name="dist"))
dist_melt[, condition1 := str_extract(id1, "batch|fermenter|AAM|M10|mouse")]
dist_melt[is.na(condition1), condition1 := "mono"]
dist_melt[condition1 == "AAM", condition1 := "AF"]
dist_melt[condition1 == "M10", condition1 := "APF"]
dist_melt[, condition2 := str_extract(id2, "batch|fermenter|AAM|M10|mouse")]
dist_melt[is.na(condition2), condition2 := "mono"]
dist_melt[condition2 == "AAM", condition2 := "AF"]
dist_melt[condition2 == "M10", condition2 := "APF"]
dist_melt[,condition:=paste0(condition1,"-",condition2)]
dist_melt[, org1 := str_extract(id1, ".*?(?=-)")]
dist_melt[, org2 := str_extract(id2, ".*?(?=-)")]
dist_sel <- dist_melt[(condition1=="mono" & condition2=="AF") |
                      (condition1=="mono" & condition2=="APF") |
                      (condition1=="mono" & condition2=="mouse") |
                      (condition1=="mouse" & condition2=="AF") |
                      (condition1=="mouse" & condition2=="APF") |
                      (condition1=="AF" & condition2=="APF")]
dist_sel[condition=="mono-AF", condition_new:="mono vs. vitro (AF)"]
#dist_sel[condition=="mono-APF", condition_new:="mono vs. vitro"]
dist_sel[condition=="mono-mouse", condition_new:="mono vs. vivo"]
#dist_sel[condition=="mouse-AF", condition_new:="vivo vs. vitro"]
dist_sel[condition=="mouse-APF", condition_new:="vivo vs. vitro (APF)"]
ggplot(dist_sel[!is.na(condition_new)], aes(x=condition_new,y=dist)) + geom_boxplot() + theme_bw(base_size=14) + 
  stat_compare_means(label = "p.signif", method = "wilcox", comparisons = as.list(as.data.frame(combinat::combn(1:3,2)))) +
  theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1)) + xlab("") + ylab("Distance (Euclidian)") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
ggsave("img/logpca-distances2.pdf", width=2, height=3)
ggsave("img/logpca-distances2.png", width=2, height=3)
