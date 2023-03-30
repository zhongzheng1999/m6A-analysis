library(tidyverse)

# sce <- readRDS("H:/Project/Main/20230105m6A/Result/code/data/sce.Rds")
# Idents(sce) <- 'CellType_int'
# sce$CellType_m6A <- paste0(sce$CellType_int,'_',sce$m6A_bi)
# n1 <- levels(sce)
# ident1 <- paste0(n1,'_Active')
# ident2 <- paste0(n1,'_Inactive')
# 
# get_diff <- function(i){
#   cat(ident1[i],'vs',ident2[i],'\n')  
#   FindMarkers(sce, group.by="CellType_m6A",
#               ident.1 = ident1[i], ident.2 = ident2[i], 
#               logfc.threshold = 0.01, min.pct = 0.01, test.use = 'MAST') %>% 
#     rownames_to_column(var = "gene") %>%
#     cbind(cluster_id = n1[i], .)
#   
#   
# }
# diff_gene <- map_dfr(c(1:12), get_diff)

diff_gene <- readRDS("H:/Project/Main/20230105m6A/Result/code/data/diff_gene_gsea.Rds")
sce.markers = diff_gene
table(sce.markers$cluster)
sce.markers <- subset(sce.markers, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
table(sce.markers$cluster)

sce.markers$threshold <- as.factor(ifelse(sce.markers$avg_log2FC > 0 , 'Up', 'Down'))
dim(sce.markers)
df <- as.data.frame(table(sce.markers$cluster_id, sce.markers$threshold))
df <- spread(df, Var2, Freq)
df$Down <- paste0('Down:',df$Down,')')
df$Up <- paste0('(Up:',df$Up)
df$Cluster <- paste0(df$Var1, ':', df$Up, ',', df$Down)


sce.markers <- merge(sce.markers,df,by.x = 'cluster_id', by.y = 'Var1')
cluster_sort = names(table(sce.markers$Cluster)[order(table(sce.markers$Cluster), decreasing = T)])
sce.markers$Cluster <- factor(sce.markers$Cluster,
                                 levels = cluster_sort)
mycol <- ggsci::pal_lancet()(2)
library(ggplot2)
ggplot(sce.markers, aes(avg_log2FC, Cluster))+
  geom_point(aes(size = -log10(p_val_adj),fill = threshold), shape = 21)+
  scale_fill_manual(values = c( 'Up' = mycol[2], 'Down' = mycol[1]))+
  theme_bw()
ggsave('DEG.pdf', width = 8, height = 6)

# 保存到文件
write.csv(sce.markers, "output_sce.markers.csv", quote = F)



#TF
dat <- read.csv('H:/Project/Main/20230105m6A/Result/code/data/auc_mtx.csv', header = T)
sce <- readRDS("H:/Project/Main/20230105m6A/Result/code/data/sce.Rds")
table(dat$Cell == rownames(sce@meta.data))
meta <- paste0(sce$CellType_int, '_', sce$m6A_bi)
dat$meta <- meta

# dat_batf <- dat[c('Cell', 'meta', 'BATF...')]
# dat_batf$m6A <- str_split(dat_batf$meta, '_', simplify = T)[,2]
# dat_batf$Cell <- str_split(dat_batf$meta, '_', simplify = T)[,1]
# colnames(dat_batf) <- c("Cell", "meta", "BATF","m6A")
# library(ggpubr)
# for (i in c('GMP-like','cDC-like','Mono-like','ProMono-like')) {
#   dat_batf_sub <- subset(dat_batf, subset = Cell == i)
#   dat_batf_sub$m6A <- factor(dat_batf_sub$m6A, levels = c('Inactive','Active'))
#   ggboxplot(dat_batf_sub, x = 'm6A', y = 'BATF', fill = 'm6A', palette = 'lancet',
#             ylab = 'BATF Regulon AUC value', title = i)+
#     stat_compare_means()+
#     theme_bw()+
#     theme(legend.position = 'null')+
#     xlab(NULL)
#   
#   ggsave(paste0(i, '_BATF.pdf'), width = 3, height = 3.5)
# }


# dat_mean <- aggregate(.~meta, mean, data = dat[,-1])
# rownames(dat_mean) <- dat_mean[,1]
# dat_mean <- dat_mean[,-1]
# 
# 
# n <- scale(dat_mean)
# # n[n > 2] = 2
# # n[n < -2] = -2
# library(reshape2)
# n1 <- melt(n)
# n2 <- n1 %>% group_by(Var1) %>% top_n(value,n = 3) %>% arrange(Var1)
# colnames(n) <- str_split(colnames(n),"\\.",simplify = T)[,1]
# colnames(n) <- paste0(colnames(n),"(+)")
# pheatmap::pheatmap(n[,n2$Var2], cluster_rows = F, cluster_cols = F,
#                    color = colorRampPalette(c(mycol[1],'white', mycol[2]))(100),
#                    filename = 'TF_pheamap.pdf', width = 15, height = 5)
# 
library(reshape2)
library(ggpubr)
dat_n1 <- melt(dat[,-1],'meta',)
n1_list <- split(dat_n1, dat_n1$variable)
df <- list()
for (i in names(n1_list)) {
  df[[i]] <- as.data.frame(compare_means(formula = value ~ meta, data = n1_list[[i]]))
}

df_bind <- do.call(rbind, df)
df_bind$C1 <- str_split(df_bind$group1,'_',simplify = T)[,1]
df_bind$C2 <- str_split(df_bind$group2,'_',simplify = T)[,1]

df_bind_sub <- df_bind[df_bind$C1 == df_bind$C2,]

mycol <- ggsci::pal_lancet()(2)
mycol1 <- colorRampPalette(c('white',mycol[2]))(4)
tf_cell_tab <- as.data.frame(table(df_bind_sub$p.signif, df_bind_sub$C2))
tf_cell_tab <- tf_cell_tab[tf_cell_tab$Var1 != 'ns',]
tf_cell_tab$Var2 <- factor(tf_cell_tab$Var2,
                           levels = c('T','B','NK','Plasma','CTL','cDC','Monocyte',
                                      'Mono-like','cDC-like','ProMono-like','GMP-like', 'LSPC'))
tf_cell_tab$Significance <- tf_cell_tab$Var1
ggbarplot(tf_cell_tab,'Var2','Freq', palette = pal_material("red")(4),
          position = position_dodge(),
          fill = 'Significance', ylab = 'Number')+rotate_x_text(angle = 30)
ggsave('tf_barplot.pdf', width = 12, height = 5)
saveRDS(tf_cell_tab,'tf_cell_Number.Rds')

# #pie
# tf_cell_tab$Percentage <- tf_cell_tab$Freq/sum(tf_cell_tab$Freq)
# cell_ind <- c('LSPC-Primed','Mono-like','LSPC-Quiescent',
#               'LSPC-Cycle','cDC-like','ProMono-like','GMP-like')
# tf_cell_tab$CellType <- ifelse(tf_cell_tab$Var2 %in% cell_ind, 'leukemic cell','immune cell')
# p1 <- ggpie(tf_cell_tab, x = 'Percentage', fill = 'CellType', color = NULL,label = NULL, palette = mycol)
# p2 <- ggpie(tf_cell_tab, x = 'Percentage', fill = 'Significance', color = NULL,label = NULL, palette = pal_material('red')(4))
# p3 <- ggpie(tf_cell_tab[tf_cell_tab$Significance == '****',], x = 'Percentage', fill = 'CellType', color = NULL,label = NULL, palette = mycol)
# p4 <- ggpie(tf_cell_tab[tf_cell_tab$Significance != '****',], x = 'Percentage', fill = 'CellType', color = NULL,label = NULL, palette = mycol)
# library(patchwork)
# (p1+p2)/(p3+p4)


# library(tidyverse)
# dat <- read.csv('../data/auc_mtx.csv', header = T)
# table(dat$Cell == rownames(sce@meta.data))
# meta <- paste0(sce$CellType_int, '_', sce$m6A_bi)
# dat$meta <- meta
# 
# dat_mean <- aggregate(.~meta, mean, data = dat[,-1])
# rownames(dat_mean) <- dat_mean[,1]
# dat_mean <- dat_mean[,-1]
# 
# dat_mean$type <- str_split(rownames(dat_mean),'_',simplify = T)[,2]
# dat_mean$cell <- str_split(rownames(dat_mean),'_',simplify = T)[,1]
# colnames(dat_mean) <- str_split(colnames(dat_mean), '\\.',simplify = T)[,1]
# dat_mean_list <- split(dat_mean, dat_mean$type) %>% lapply(.,function(x){ x <- melt(x,'cell')})
# for (i in names(dat_mean_list)) {
#   colnames(dat_mean_list[[i]]) <- paste0(i,"_",colnames(dat_mean_list[[i]]))
#   dat_mean_list[[i]]$pair <- paste0(dat_mean_list[[i]][,1],'_',dat_mean_list[[i]][,2])
# }
# 
# load("H:/Project/Smart/20230105m6A/step3/tf_cell_deg.Rdata")
# df_bind_sub$gene <- str_split(rownames(df_bind_sub), '\\.',simplify = T)[,1]
# df_bind_sub$pair <- paste0(df_bind_sub$C1,'_', df_bind_sub$gene)
# 
# dat_deg_merge <- merge(dat_mean_list[[1]],df_bind_sub,by = 'pair')
# dat_deg_merge <- merge(dat_mean_list[[2]],dat_deg_merge,by = 'pair')
# 
# dat_deg_merge$FC <- as.numeric(dat_deg_merge$Active_value)/as.numeric(dat_deg_merge$Inactive_value)
# dat_deg_merge$logFC <- log2(dat_deg_merge$FC)
# 
# dat_deg_merge[which(dat_deg_merge$p < 0.01 & dat_deg_merge$logFC <= -0.15 ),'sig'] <- 'Down'
# dat_deg_merge[which(dat_deg_merge$p < 0.01 & dat_deg_merge$logFC >= 0.15),'sig'] <- 'Up'
# dat_deg_merge[which(dat_deg_merge$p >= 0.01 | abs(dat_deg_merge$logFC) < 0.15 ),'sig'] <- 'None'
# sub = paste0('Up:',length(which((dat_deg_merge[,25]=='Up'))),' ', 'Down:',length(which((dat_deg_merge[,25]=='Down'))))
# cell_ind <- c('LSPC','Mono-like',
#               'cDC-like','ProMono-like','GMP-like')
# dat_deg_merge$CellType <- ifelse(dat_deg_merge$C1 %in% cell_ind, 'leukemic cell','immune cell')

# p <- ggplot(dat_deg_merge, aes(x = logFC, y = -log10(p), color = sig)) +
#   geom_point(alpha = 0.6, size = 1.5, aes(shape = CellType)) + 
#   scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('Up', 'Down', 'None')) +
#   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
#   theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
#   geom_vline(xintercept = c( 0.15), color = 'gray', size = 0.3) +
#   geom_hline(yintercept = -log(0.01, 10), color = 'gray', size = 0.3) +
#   #
#   labs(x = '\nLog2 FoldChange', y = '-Log10(Pvalue)\n', color = '', title = 'Active-Inactive\n',subtitle = sub)
# p
# 
# ggsave('A-I_volcano.pdf', p, width = 4, height = 5)
# 
# write.csv(dat_deg_merge,'tf_deg_merge.csv')
# saveRDS(dat_deg_merge,'dat_deg_merge.Rds')


diff_gene <- readRDS("H:/Project/Main/20230105m6A/Result/code/data/diff_gene_gsea.Rds")
sce.markers = diff_gene
table(sce.markers$cluster)
sce.markers <- subset(sce.markers, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
table(sce.markers$cluster)

tf_tg <- readRDS("H:/Project/Main/20230105m6A/Result/code/data/tf_tg.Rds")
tf_tg_list <- split(tf_tg,tf_tg$TG)
tf_tg_list <- lapply(tf_tg_list, function(x){
  x = x$TF
  return(x)
})

sce.markers$TF = NA
for (i in sce.markers$gene) {
  sce.markers[sce.markers$gene == i, 'TF'] <- paste(tf_tg_list[[i]],collapse = ',')
}

sce.markers_row <- separate_rows(sce.markers, 'TF', sep = ',')
sce.markers_list <- split(sce.markers_row,sce.markers_row$cluster_id)
# openxlsx::write.xlsx(sce.markers_list, file = 'DEG_tf.xlsx')

dat_deg_merge <- readRDS("H:/Project/Main/20230105m6A/Result/code/0327Figure/Figure2/dat_deg_merge.Rds")
sce.markers_row$ind <- paste0(sce.markers_row$cluster_id,'_',sce.markers_row$TF)

sce.markers_row <- merge(sce.markers_row,dat_deg_merge,
                         by.x = 'ind',by.y = 'pair',all.x = T)
sce.markers_row = sce.markers_row[c(1:9,12,15,19,22,23,27:29)]
saveRDS(sce.markers_row,'DEG_tf_sig.Rds')

# writeClipboard(unique(sce.markers_list$`ProMono-like`$gene))

library(tidyverse)
setwd("H:/Project/Main/20230105m6A/Result/code/0327Figure/Figure2")
DEG_tf <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/DEG_tf_sig.Rds")

DEG_tf <- subset(DEG_tf, p.signif != 'ns' & CellType == "leukemic cell")
DEG_tf_cyto <- DEG_tf[c(2,9,3)]
write.csv(DEG_tf_cyto, 'DEG_tf_cyto.csv',quote = F)

DEG_tf_list <- split(DEG_tf,DEG_tf$cluster_id)
DEG_tf_cyto_list <- split(DEG_tf_cyto,DEG_tf_cyto$cluster_id)
for (i in c('cDC-like','LSPC','GMP-like','Mono-like','ProMono-like')) {
  tf_id <- names(sort(table(DEG_tf_cyto_list[[i]]$TF), decreasing = T))[1:10]
  tf_tg <- DEG_tf_cyto_list[[i]][DEG_tf_cyto_list[[i]]$TF %in% tf_id,]
  DEG_tf_cyto_list[[i]] <- tf_tg 
}
write.xlsx(DEG_tf_cyto_list, 'DEG_tf_cyto_list.xlsx')

for (i in c('cDC-like','LSPC','GMP-like','Mono-like','ProMono-like')) {
  pathway_id = paste0("H:/Project/Main/20230105m6A/Result/code/test/",i,".xlsx")
  cDC <- read.xlsx(pathway_id, sheet = 2)
  cDC <- cDC[grep('Summary',cDC$GroupID)[1:5],c(4,9)]
  cDC_row <- separate_rows(cDC,'Symbols', sep = ',')
  colnames(cDC_row) = c("TF","gene.x") 
  cDC_row$cluster_id = 'ALL'
  tf_id <- names(sort(table(DEG_tf_cyto_list[[i]]$TF), decreasing = T))[1:10]
  tf_tg <- DEG_tf_cyto_list[[i]][DEG_tf_cyto_list[[i]]$TF %in% tf_id,]
  DEG_tf_cyto_list[[i]] <- rbind(tf_tg,cDC_row)
}

write.xlsx(DEG_tf_cyto_list, 'DEG_tf_cyto_pathway.xlsx')
writeClipboard(unique(DEG_tf_list$`cDC-like`$gene.x))