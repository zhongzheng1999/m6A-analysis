library(ggpubr)
library(tidyverse)
# dat <- read.csv('../step4/result_type/tcga_fpkm_sweep/p/pearson_pam/3.csv',
#                 header = T, row.names = 1)
# dat$Subtype <- paste0('C',dat$Cluster)

meta <- read.csv("H:/Project/Main/20230105m6A/two/workflow/TCGA_sweep_vst_surv/all/pearson_pam/3.csv",
                 header = T, row.names = 1)
meta$TCGA_Subtype <- paste0('C',meta$Subtype)

df <- read.table('tcga_count.estimate_score.txt', header = T, row.names = 1, sep = '\t')
df <- as.data.frame(t(df))

# rownames(df) <- str_sub(df$NAME, 1,15)
rownames(df) <- df$NAME
colnames(df) <- df[1,]
df <- df[-1,-1]

sample_id <- intersect(rownames(meta),rownames(df))
meta <- meta[sample_id,]
df <- df[sample_id,]
rownames(meta) == rownames(df)

#"StromalScore"  "ImmuneScore"   "ESTIMATEScore"
df <- as.data.frame(apply(df, 2, as.numeric))
df$Subtype <- meta$TCGA_Subtype

compare_list <- as.list(as.data.frame(combn(unique(df$Subtype),2)))
ggviolin(df,'Subtype','ESTIMATEScore', fill = 'Subtype',
         palette = 'lancet',
         add = 'boxplot', add.params = list(fill = 'white'))+
  stat_compare_means(comparisons = compare_list, label = 'p.signif')+
  stat_compare_means(label.y = 6000)
ggsave('TCGA_ESTIMATEScore.pdf', width = 6,height = 5)




#####pathway
library(clusterProfiler)
library(GSEABase)
library(GSVA)
library(ggpubr)
library(tidyverse)

load("H:/Project/Main/20230105m6A/two/subtype/data/TCGA_subt.Rdata")
dat <- expr
genelist <- read.table('GeneList.txt', header = T,sep = '\t')
pathway <- split(genelist,genelist$Category)
pathway <- lapply(pathway, function(x) x <- x$Symbol)

# colnames(dat) <- str_sub(colnames(dat),1,15)
ssgsea_result <- gsva(as.matrix(dat), pathway, method = 'gsva')

sample_id <- intersect(colnames(ssgsea_result),rownames(meta))
ssgsea_result <- ssgsea_result[,sample_id]
meta <- meta[sample_id,]
colnames(ssgsea_result) == rownames(meta)

pathway_dat <- cbind(meta, t(ssgsea_result))
pathway_dat$Subtype <- factor(pathway_dat$TCGA_Subtype, levels = c('C1','C2','C3'))
compare_list <- as.list(as.data.frame(combn(unique(meta$Subtype),2)))
compare_list <- compare_list[c(3,2,1)]
for (i in rownames(ssgsea_result)) {
  tmp <- pathway_dat[c('Subtype',i)]
  colnames(tmp) <- c('Subtype','Score')
  ggviolin(tmp, x = 'Subtype', y = 'Score',fill = 'Subtype',
           add = 'boxplot', add.params = list(fill = 'white'),
           ylab = 'Pathway Enrich Score', palette = 'lancet')+
    stat_compare_means(comparisons = compare_list, label = 'p.signif')+
    stat_compare_means(label.y = 0.7)
  ggsave(paste0(i,'_pathway.pdf'), width = 6, height = 5)
}

meta$TCGA_Subtype <- as.factor(meta$TCGA_Subtype)
n <- t(scale(t(ssgsea_result)))
n[n > 2] = 2
n[n < -2] = -2

mycol <- ggsci::pal_lancet()(2)
pheatmap::pheatmap(n[,order(meta$TCGA_Subtype)],annotation_col = meta[c('TCGA_Subtype')],
                   cluster_cols = F, show_colnames = F,
                   color = colorRampPalette(c(mycol[1], 'white', mycol[2]))(100),
                   filename = 'pathway_imm_gsva.pdf',gaps_col = c(52,98,132),
                   width = 10, height = 3)



#immGene
load("H:/Project/Main/20230105m6A/two/subtype/data/TCGA_subt.Rdata")
dat <- expr
gene <- openxlsx::read.xlsx('../1-s2.0-S1074761318301213-mmc7.xlsx',sheet = 1)
gene <- gene[1:78,]

dat <- as.data.frame(dat)
# colnames(dat) <- str_sub(colnames(dat),1,15)

dat_imm <- dat[gene$HGNC.Symbol,]
dat_imm <- na.omit(dat_imm)

sample_id <- intersect(colnames(dat_imm),rownames(meta))
dat_imm <- dat_imm[,sample_id]
meta <- meta[sample_id,]
colnames(dat_imm) == rownames(meta)
dat_imm <- as.data.frame(t(dat_imm))
rownames(dat_imm) == rownames(meta)
dat_imm$Subtype <- meta$TCGA_Subtype
dat_imm_mean <- aggregate(. ~ Subtype,mean,data = dat_imm)

dat_imm_n <- as.data.frame(t(dat_imm_mean))
colnames(dat_imm_n) <- dat_imm_n[1,]
dat_imm_n <- dat_imm_n[-1,]
dat_imm_n <- apply(dat_imm_n, 2, as.numeric)
rownames(dat_imm_n) <- colnames(dat_imm_mean)[-1]

n <- t(scale(t(dat_imm_n)))
n[n < -2] = -2
n[n > 2] = 2
anno_row <- data.frame(row.names = gene$HGNC.Symbol, gene[5:7])
anno_row <- anno_row[rownames(n),]
# pheatmap::pheatmap(n[order(anno_row$Super.Category),order(meta$Subtype)],
#                    annotation_col = meta[c('Subtype')],
#                    annotation_row = anno_row,
#                    cluster_cols = F, show_colnames = F,
#                    filename = 'immgene_pheatmap_row.pdf',
#                    width = 10, height = 11)

pheatmap::pheatmap(n[order(anno_row$Super.Category),],
                   annotation_row = anno_row,
                   cluster_rows = F,
                   color = colorRampPalette(c(mycol[1], 'white', mycol[2]))(100),
                   cluster_cols = F, show_colnames = T,
                   filename = 'immgene_pheatmap_mean.pdf',
                   width = 5, height = 11)


cluster_dir = 'H:/Project/Main/20230105m6A/two/workflow/TCGA_sweep_vst_surv/all/pearson_pam/3.csv'
tide <- read.csv("H:/Project/Main/20230105m6A/two/imm/data/TCGA_TIDE.csv",
                 row.names = 1, header = T)
# load("H:/Project/Smart/20230105m6A/two/subtype_DEG/gene_list.Rdata")
subt <- read.csv(cluster_dir, check.names = F, stringsAsFactors = F, header = T, row.names = 1)

sample_id <- intersect(rownames(tide),rownames(subt))
tide <- tide[sample_id,]
subt <- subt[sample_id,]
subt$Subtype <- paste0('C',subt$Subtype)
rownames(tide) == rownames(subt)


patient_ICBresponse <- tide
ICBresponse <- patient_ICBresponse
##统计免疫响应患者数
table(ICBresponse$Responder=="False")
# FALSE  TRUE 
# 95   107
ICBresponse$Responder <- ifelse(str_detect(ICBresponse$Responder,"False"),"NR","R")
ICBresponse1 <- dplyr::pull(ICBresponse, 3) ##将data.frame中的一列转换为向量
names(ICBresponse1) <- rownames(patient_ICBresponse)
ICBresponse1[1:10]
# save(ICBresponse,file = "ICBresponse.Rdata")
# 
# 
# load("ICBresponse.Rdata")
##发散条形图
## barplot
dat_plot <- data.frame(id = rownames(ICBresponse),
                       t = ICBresponse$TIDE)

# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$t > 0,'NR','R'),levels=c('R','NR'))
# 排序
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)
library(ggthemes)
# install.packages("ggprism")
library(ggprism)
p5 <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  # scale_x_continuous(breaks = seq(0,50,100))+ # X 轴每隔 50 个单位显示一个刻度
  scale_fill_manual(values = c('R'= '#FF1493','NR'='#40E0D0')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('TIDE score') + 
  guides(fill=guide_legend(key.width = 3, key.height = 5, nrow = 2, ncol = 1, byrow = TRUE))+ #显示图例
  theme_prism(border = T) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(.95, .5),#图例位置
    legend.justification = c("right", "top")#固定右上角
  )
p5
ggsave(plot = p5, filename = 'TCGA_TIDE.pdf', width = 10, height = 4)



#####TEST
load("H:/Project/Main/20230105m6A/two/subtype/data/GSE146173_subt.Rdata")

df <- read.table('GSE146173.gene_expression_RNAseq.Salmon_Estimated_Counts.estimate_score.txt', header = T, row.names = 1, sep = '\t')
df <- as.data.frame(t(df))

# rownames(df) <- str_sub(df$NAME, 1,15)
rownames(df) <- df$NAME
colnames(df) <- df[1,]
df <- df[-1,-1]

sample_id <- intersect(rownames(subt),rownames(df))
subt <- subt[sample_id,]
df <- df[sample_id,]
rownames(subt) == rownames(df)

#"StromalScore"  "ImmuneScore"   "ESTIMATEScore"
df <- as.data.frame(apply(df, 2, as.numeric))
df$Subtype <- subt$res_cluster

compare_list <- as.list(as.data.frame(combn(unique(df$Subtype),2)))
ggviolin(df,'Subtype','ImmuneScore', fill = 'Subtype',
         palette = 'lancet',
         add = 'boxplot', add.params = list(fill = 'white'))+
  stat_compare_means(comparisons = compare_list, label = 'p.signif')+
  stat_compare_means(label.y = 5000)
ggsave('GSE146173_ImmuneScore.pdf', width = 6,height = 5)




#####pathway
library(clusterProfiler)
library(GSEABase)
library(GSVA)
library(ggpubr)
library(tidyverse)

load("H:/Project/Main/20230105m6A/two/subtype/data/GSE146173_vst.Rdata")
dat <- GSE146173_fpkm
genelist <- read.table('GeneList.txt', header = T,sep = '\t')
pathway <- split(genelist,genelist$Category)
pathway <- lapply(pathway, function(x) x <- x$Symbol)

# colnames(dat) <- str_sub(colnames(dat),1,15)
ssgsea_result <- gsva(as.matrix(dat), pathway, method = 'gsva')

sample_id <- intersect(colnames(ssgsea_result),rownames(subt))
ssgsea_result <- ssgsea_result[,sample_id]
subt <- subt[sample_id,]
colnames(ssgsea_result) == rownames(subt)

pathway_dat <- cbind(subt, t(ssgsea_result))
pathway_dat$Subtype <- factor(pathway_dat$res_cluster, label = c('C1','C2','C3'), levels = c('1','2','3'))
compare_list <- as.list(as.data.frame(combn(unique(pathway_dat$Subtype),2)))
compare_list <- compare_list[c(3,2,1)]
# for (i in rownames(ssgsea_result)) {
#   tmp <- pathway_dat[c('Subtype',i)]
#   colnames(tmp) <- c('Subtype','Score')
#   ggviolin(tmp, x = 'Subtype', y = 'Score',fill = 'Subtype',
#            add = 'boxplot', add.params = list(fill = 'white'),
#            ylab = 'Pathway Enrich Score', palette = 'lancet')+
#     stat_compare_means(comparisons = compare_list, label = 'p.signif')+
#     stat_compare_means(label.y = 0.7)
#   ggsave(paste0(i,'_pathway.pdf'), width = 6, height = 5)
# }

subt$TCGA_Subtype <- as.factor(subt$res_cluster)
n <- t(scale(t(ssgsea_result)))
n[n > 2] = 2
n[n < -2] = -2
table(subt$TCGA_Subtype)
mycol <- ggsci::pal_lancet()(2)
pheatmap::pheatmap(n[,order(subt$TCGA_Subtype)],annotation_col = subt[c('TCGA_Subtype')],
                   cluster_cols = F, show_colnames = F,
                   color = colorRampPalette(c(mycol[1], 'white', mycol[2]))(100),
                   filename = 'GSE146173_pathway_imm_gsva.pdf',gaps_col = c(103,212),
                   width = 10, height = 3)



#immGene
load("H:/Project/Main/20230105m6A/two/subtype/data/GSE146173_vst.Rdata")
dat <- GSE146173_fpkm
gene <- openxlsx::read.xlsx('../1-s2.0-S1074761318301213-mmc7.xlsx',sheet = 1)
gene <- gene[1:78,]

dat <- as.data.frame(dat)
# colnames(dat) <- str_sub(colnames(dat),1,15)

dat_imm <- dat[gene$HGNC.Symbol,]
dat_imm <- na.omit(dat_imm)

sample_id <- intersect(colnames(dat_imm),rownames(subt))
dat_imm <- dat_imm[,sample_id]
subt <- subt[sample_id,]
colnames(dat_imm) == rownames(subt)
dat_imm <- as.data.frame(t(dat_imm))
rownames(dat_imm) == rownames(subt)
dat_imm$Subtype <- subt$TCGA_Subtype
dat_imm_mean <- aggregate(. ~ Subtype,mean,data = dat_imm)

dat_imm_n <- as.data.frame(t(dat_imm_mean))
colnames(dat_imm_n) <- dat_imm_n[1,]
dat_imm_n <- dat_imm_n[-1,]
dat_imm_n <- apply(dat_imm_n, 2, as.numeric)
rownames(dat_imm_n) <- colnames(dat_imm_mean)[-1]

n <- t(scale(t(dat_imm_n)))
n[n < -2] = -2
n[n > 2] = 2
anno_row <- data.frame(row.names = gene$HGNC.Symbol, gene[5:7])
anno_row <- anno_row[rownames(n),]
# pheatmap::pheatmap(n[order(anno_row$Super.Category),order(meta$Subtype)],
#                    annotation_col = meta[c('Subtype')],
#                    annotation_row = anno_row,
#                    cluster_cols = F, show_colnames = F,
#                    filename = 'immgene_pheatmap_row.pdf',
#                    width = 10, height = 11)

pheatmap::pheatmap(n[order(anno_row$Super.Category),],
                   annotation_row = anno_row,
                   cluster_rows = F,
                   color = colorRampPalette(c(mycol[1], 'white', mycol[2]))(100),
                   cluster_cols = F, show_colnames = T,
                   filename = 'GSE146173_immgene_pheatmap_mean.pdf',
                   width = 5, height = 11)



tide <- read.csv("H:/Project/Main/20230105m6A/two/imm/data/GSE146173_TIDE.csv",
                 row.names = 1, header = T)
# load("H:/Project/Smart/20230105m6A/two/subtype_DEG/gene_list.Rdata")
# subt <- read.csv(cluster_dir, check.names = F, stringsAsFactors = F, header = T, row.names = 1)

sample_id <- intersect(rownames(tide),rownames(subt))
tide <- tide[sample_id,]
subt <- subt[sample_id,]
subt$Subtype <- paste0('C',subt$res_cluster)
rownames(tide) == rownames(subt)


patient_ICBresponse <- tide
ICBresponse <- patient_ICBresponse
##统计免疫响应患者数
table(ICBresponse$Responder=="False")
# FALSE  TRUE 
# 95   107
ICBresponse$Responder <- ifelse(str_detect(ICBresponse$Responder,"False"),"NR","R")
ICBresponse1 <- dplyr::pull(ICBresponse, 3) ##将data.frame中的一列转换为向量
names(ICBresponse1) <- rownames(patient_ICBresponse)
ICBresponse1[1:10]
# save(ICBresponse,file = "ICBresponse.Rdata")
# 
# 
# load("ICBresponse.Rdata")
##发散条形图
## barplot
dat_plot <- data.frame(id = rownames(ICBresponse),
                       t = ICBresponse$TIDE)

# 新增一列 根据t阈值分类
dat_plot$threshold = factor(ifelse(dat_plot$t > 0,'NR','R'),levels=c('R','NR'))
# 排序
dat_plot <- dat_plot %>% arrange(t)
# 变成因子类型
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
# 绘制
library(ggplot2)
library(ggthemes)
# install.packages("ggprism")
library(ggprism)
p5 <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  # scale_x_continuous(breaks = seq(0,50,100))+ # X 轴每隔 50 个单位显示一个刻度
  scale_fill_manual(values = c('R'= '#FF1493','NR'='#40E0D0')) +
  geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('TIDE score') + 
  guides(fill=guide_legend(key.width = 3, key.height = 5, nrow = 2, ncol = 1, byrow = TRUE))+ #显示图例
  theme_prism(border = T) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(.95, .5),#图例位置
    legend.justification = c("right", "top")#固定右上角
  )
p5
ggsave(plot = p5, filename = 'GSE146173_TIDE.pdf', width = 10, height = 4)
