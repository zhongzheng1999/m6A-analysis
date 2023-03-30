
library(tidyverse)
library(ggpubr)
load("H:/Project/Main/20230105m6A/two/drug/predictedBoxdat.Rdata")
tcga_plotp = plotp
tcga_pred = predictedBoxdat

load("H:/Project/Main/20230105m6A/two/drug/GSE106291_predictedBoxdat.Rdata")
gse106_plotp = plotp
gse106_pred = predictedBoxdat

load("H:/Project/Main/20230105m6A/two/drug/GSE146173_predictedBoxdat.Rdata")
gse146_plotp = plotp
gse146_pred = predictedBoxdat

drug_list = list(
  tcga = tcga_pred,
  gse106 = gse106_pred,
  gse146 = gse146_pred
)

drug_list <- lapply(drug_list, function(x){
  dat_drug <- do.call(rbind, x)
  dat_drug$drug <- str_split(rownames(dat_drug), '\\.T|\\.G|\\.N', simplify = T)[,1]
  return(dat_drug)
})


drug_list <- lapply(drug_list, function(x){
  compare_dat <- compare_means(est.ic50~ImmClust,data = x, group.by = 'drug', method = "kruskal.test")
  compare_dat <- compare_dat[compare_dat$p.signif != 'ns',]
  return(compare_dat$drug)
})

drug_id = Reduce(intersect, drug_list)


drug_plot_list = list(
  tcga = tcga_plotp,
  gse106 = gse106_plotp,
  gse146 = gse146_plotp
)
for (i in names(drug_plot_list)) {
  p2 <- plot_grid(plotlist=plotp[drug_id], ncol=5)
  ggsave(plot = p2,paste0(i, "_drug_multiple.pdf"), width = 12, height = 17)
}


drug_list = list(
  tcga = lapply(tcga_pred, function(x){
    x$id = rownames(x)
    return(x)
  }),
  gse106 = lapply(gse106_pred, function(x){
    x$id = rownames(x)
    return(x)
  }),
  gse146 = lapply(gse146_pred, function(x){
    x$id = rownames(x)
    return(x)
  })
)

drug_list <- lapply(drug_list, function(x){
  dat_drug <- do.call(rbind, x)
  dat_drug$drug <- str_split(rownames(dat_drug), '\\.T|\\.G|\\.N', simplify = T)[,1]
  dat_drug <- dat_drug[dat_drug$drug %in% drug_id,]
  return(dat_drug)
})


drug_list <- lapply(drug_list, function(x){
  x <- x %>% 
    select(-"ImmClust") %>% 
    spread(key = 'id', value = 'est.ic50')
  rownames(x) = x$drug
  x <- x[,-1]
  return(x)
})

meta_list = list(
  tcga = lapply(tcga_pred, function(x){
    x$id = rownames(x)
    return(x)
  }),
  gse106 = lapply(gse106_pred, function(x){
    x$id = rownames(x)
    return(x)
  }),
  gse146 = lapply(gse146_pred, function(x){
    x$id = rownames(x)
    return(x)
  })
)

meta_list <- lapply(meta_list, function(x){
  dat_drug <- do.call(rbind, x)
  dat_drug <- dat_drug[!duplicated(dat_drug$id),]
  meta <- data.frame(row.names = dat_drug$id, Clust = dat_drug$ImmClust)
  return(meta)
})

library(pheatmap)

for (i in names(meta_list)) {
  col_anno <- meta_list[[i]] %>% arrange(Clust)
  n <- t(scale(t(drug_list[[i]])))
  n[n > 2] = 2
  n[n < -2] = -2
  mycol <- ggsci::pal_lancet()(2)
  pheatmap(n[drug_id,rownames(col_anno)],
           color = colorRampPalette(c(mycol[1], 'white', mycol[2]))(100),
           cluster_cols = F, annotation_col = col_anno,
           show_colnames = F,cluster_rows = F,
           filename = paste0(i,'_drug_ic50_1.pdf'), width = 7, height = 5)
}



drug_plot_list = list(
  tcga = tcga_plotp,
  gse106 = gse106_plotp,
  gse146 = gse146_plotp
)
drug_boxid = 'Nutlin.3a'
for (i in names(drug_plot_list)) {
  plotp <- drug_plot_list[[i]]
  p2 <- plotp[[drug_boxid]]
  ggsave(plot = p2,filename = paste0(i, '_', drug_boxid,"_ic50.pdf"), width = 3, height = 3)
}
