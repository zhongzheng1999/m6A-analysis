rm(list = ls())
require(maftools)
options(stringsAsFactors = F) 
library(data.table)
tmp=fread('TCGA-LAML.mutect2_snv.tsv')
head(tmp)   
colnames(tmp) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                  "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "HGVSp_Short" , 'effect' ,"Consequence",
                  "vaf" ) 
tmp$Entrez_Gene_Id =1
tmp$Center ='ucsc'
tmp$NCBI_Build ='GRCh38'
tmp$NCBI_Build ='GRCh38'
tmp$Strand ='+'
tmp$Variant_Classification = tmp$effect
tail(sort(table(tmp$Variant_Classification )))
tmp$Tumor_Seq_Allele1 = tmp$Reference_Allele
tmp$Variant_Type = ifelse(
  tmp$Reference_Allele %in% c('A','C','T','G') & tmp$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
  'SNP','INS'
)
table(tmp$Variant_Type)

tmp_1 <- tmp %>% as_tibble() %>% 
  separate_rows(Variant_Classification, sep = ";")
tmp_1$Tumor_Sample_Barcode <- str_sub(tmp_1$Tumor_Sample_Barcode,1,12)
tcga.brca = read.maf(maf = tmp_1,isTCGA = T,
                     vc_nonSyn=names(sort(table(tmp$Variant_Classification ))))

oncoplot(maf = tcga.brca, top = 20) # 高频突变的前10个基因

subt <- read.csv("H:/Project/Smart/20230105m6A/step4/result_type/tcga_fpkm_sweep/p/pearson_pam/3.csv", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
head(subt)
subt$TCGA_Subtype <- paste0('C',subt$Cluster)
subt$Tumor_Sample_Barcode <- subt$X_PATIENT

tcga.brca = read.maf(maf = tmp_1,isTCGA = T,clinicalData = subt,
                     vc_nonSyn=names(sort(table(tmp$Variant_Classification ))))

oncoplot(maf = tcga.brca, top = 20, clinicalFeatures = 'TCGA_Subtype', sortByAnnotation = T,
         removeNonMutated = F, logColBar = T) 

dat <- tmb(tcga.brca)
dat <- merge(dat, subt, by.x = 'Tumor_Sample_Barcode', by.y = 'X_PATIENT')

compare_list <- list(c('C1','C2'),c('C1','C3'),c('C2','C3'))
dat$TCGA_Subtype <- factor(dat$TCGA_Subtype, levels = c('C1','C2','C3'))
ggboxplot(dat, x = 'TCGA_Subtype', y = 'total_perMB_log', fill = 'TCGA_Subtype',
          palette = 'lancet', ylab = 'TMB/MB (log10)', title = 'Mutation Burden')+
  stat_compare_means(comparisons = compare_list)
ggsave('TMB_subtype_boxplot.pdf', width = 4, height = 6)


####fisher
dat_c1_subt <- subt[subt$TCGA_Subtype == 'C1',]
dat_c2_subt <- subt[subt$TCGA_Subtype == 'C2',]
dat_c3_subt <- subt[subt$TCGA_Subtype == 'C3',]
tmp_1$sample_id <- str_sub(tmp_1$Tumor_Sample_Barcode,1,15)
tmp_c1 <- tmp_1[tmp_1$sample_id %in% dat_c1_subt$Tumor_Sample_Barcode,]
tmp_c2 <- tmp_1[tmp_1$sample_id %in% dat_c2_subt$Tumor_Sample_Barcode,]
tmp_c3 <- tmp_1[tmp_1$sample_id %in% dat_c3_subt$Tumor_Sample_Barcode,]
maf_c1 <- read.maf(maf = tmp_c1, clinicalData = dat_c1_subt, isTCGA = T,
                   vc_nonSyn=names(sort(table(tmp_c1$Variant_Classification ))))
maf_c2 <- read.maf(maf = tmp_c2, clinicalData = dat_c2_subt, isTCGA = T,
                   vc_nonSyn=names(sort(table(tmp_c2$Variant_Classification ))))
maf_c3 <- read.maf(maf = tmp_c3, clinicalData = dat_c3_subt, isTCGA = T,
                   vc_nonSyn=names(sort(table(tmp_c3$Variant_Classification ))))

mafCompare(maf_c1, maf_c2)
oncoplot(maf_c1)
oncoplot(maf_c2)
oncoplot(maf_c3)


C1_gene_prop <- apply(maf_c1@gene.summary[,-1],1,sum)/length(maf_c1@variants.per.sample)
C1_gene_prop <- data.frame(gene = maf_c1@gene.summary$Hugo_Symbol, prop_c1 = C1_gene_prop)
C2_gene_prop <- apply(maf_c2@gene.summary[,-1],1,sum)/length(maf_c2@variants.per.sample)
C2_gene_prop <- data.frame(gene = maf_c2@gene.summary$Hugo_Symbol, prop_c2 = C2_gene_prop)
C3_gene_prop <- apply(maf_c3@gene.summary[,-1],1,sum)/length(maf_c3@variants.per.sample)
C3_gene_prop <- data.frame(gene = maf_c3@gene.summary$Hugo_Symbol, prop_c3 = C3_gene_prop)

prop_merge <- merge(C1_gene_prop,C2_gene_prop,all = T)
prop_merge <- merge(prop_merge,C3_gene_prop,all = T)
write.csv(prop_merge,'prop_merge.csv')


# 加载focal和broad data
focalload <- read.table("TCGA-LAML.gistic.tsv", header = T, sep = "\t", check.names = F)
rownames(focalload) <- focalload$`Gene Symbol`
head(focalload)[, 1:8]
broadload <- read.table("TCGA-LAML.gistic.tsv", header = T, sep = "\t", check.names = F)
rownames(broadload) <- broadload$`Gene Symbol`

# 加载分组信息
subt <- read.csv("H:/Project/Smart/20230105m6A/step4/result_type/tcga_fpkm_sweep/p/pearson_pam/3.csv", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
head(subt)
subt$TCGA_Subtype <- paste0('C',subt$Cluster)
subt$Tumor_Sample_Barcode <- gsub('\\.','-',rownames(subt))

## focal
focalload <- focalload[, -1]
focalgainload <- focalload
focallossload <- focalload
# gain 
focalgainload[focalgainload > 0]  <- 1
focalgainload[focalgainload < 0]  <- 0
focalgainload <- data.frame(colSums(focalgainload))
colnames(focalgainload) <- "focal_gain_load"
# loss
focallossload[focallossload > 0]  <- 0
focallossload[focallossload < 0]  <- 1
focallossload <- data.frame(colSums(focallossload))
colnames(focallossload) <- "focal_loss_load"

## broad
broadload <- broadload[, -1]
broadload_gain <- broadload
broadload_loss <- broadload
# gain 
broadload_gain[broadload_gain > 0] <- 1
broadload_gain[broadload_gain < 0] <- 0
broadload_gain <- data.frame(colSums(broadload_gain))
colnames(broadload_gain) <- "broad_gain_load"
# loss 
broadload_loss[broadload_loss > 0] <- 0
broadload_loss[broadload_loss < 0] <- 1
broadload_loss <- data.frame(colSums(broadload_loss))
colnames(broadload_loss) <- "broad_loss_load"

## 把focal和broad的gain load和loss load这四列merge在一起
copyloadlist <- list(focalgainload = focalgainload, focallossload = focallossload,
                     broadload_gain = broadload_gain, broadload_loss = broadload_loss)
for (i in 1:length(copyloadlist)){
  tmpdata <- copyloadlist[[i]]
  copyloadlist[[i]] <- data.frame(barcode = rownames(tmpdata), tmpdata)
}
copyload <- Reduce(function(x, y) merge(x = x, y = y, by = "barcode"), 
                   copyloadlist)
head(copyload)


## 再跟分组信息merge在一起
rownames(copyload) <- str_sub(copyload$barcode,1,15)
copyload <- data.frame(sampleID = rownames(copyload),
                       copyload)
copyload <- merge(copyload, subt[, c("TCGA_Subtype", "Tumor_Sample_Barcode")],
                  by.y = "Tumor_Sample_Barcode", by.x = 'sampleID')
copyload$TCGA_Subtype <- as.character(copyload$TCGA_Subtype)
head(copyload)
table(copyload$TCGA_Subtype)

# 让图中的三组数据按顺序排列
copyload$TCGA_Subtype <- factor(copyload$TCGA_Subtype, levels = c("C1", 
                                                                  "C2",
                                                                  "C3"))
# 对比时的分组
my_comparisons <- list(c("C1", "C2"),
                       c("C1", "C3"),
                       c("C2", "C3"))

# 画图，写在list里
copyplot <- lapply(1:4, function(i){
  loadname <- colnames(copyload)[3:6]
  plotdata <- copyload[, c(loadname[i], "TCGA_Subtype")]
  colnames(plotdata) <- c("loadtype", "TCGA_Subtype")
  ylabel <- c("Burden of Copy Number Gain", "Burden of Copy Number Loss",
              "Burden of Copy Number Gain", "Burden of Copy Number Loss")
  mainlabel <- c("Focal", "Focal", "Arm-level", "Arm-level")
  label <- ggdraw() + draw_label(mainlabel[i])
  ymax <- max(plotdata$loadtype)
  plot_grid(label, 
            ggpubr::ggboxplot(plotdata, x = "TCGA_Subtype", y = "loadtype", 
                              fill = "TCGA_Subtype", palette =c("brown", "yellow", "grey"),
                              ggtheme = theme_bw()) +
              theme(panel.grid =element_blank()) +
              xlab("") +
              ylab(ylabel[i]) +
              theme(legend.position = "none")+
              theme(axis.title.y = element_text(size = rel(1.25)),
                    axis.text.x = element_text(size = rel(1.25)),
                    axis.text.y = element_text(size = rel(1.25))) +
              stat_compare_means(comparisons = my_comparisons, label = "p.signif") + 
              stat_compare_means(label.y = 1.45*ymax,
                                 label.x = 0.8),
            ncol=1, rel_heights=c(.2, 1))
})
cowplot::plot_grid(plotlist = copyplot, ncol = 2)



###CNV
cnv <- read.table('./LAML.cnv.tsv', header = T)
rownames(cnv) <- cnv$symbol
cnv <- cnv[4:194]

cnv[cnv == '-2'] <- 'Homozygous_Deletion'
cnv[cnv == '-1'] <- 'Single_Copy_Deletion'
cnv[cnv == '0'] <- ''
cnv[cnv == '1'] <- 'Low_Amplification'
cnv[cnv == '2'] <- 'High_Amplification'

col <- c( "Single_Copy_Deletion" = "#0000CC","Homozygous_Deletion" = "#0000E5",
          "Low_Amplification" = "#CC0000", "High_Amplification" = "#E50000")
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, 
              gp = gpar(fill = "#e8e7e3", col = NA))
  },
  Single_Copy_Deletion = function(x, y, w, h) {
    grid.rect(x, y, w, h*0.9,
              gp = gpar(fill = col["Single_Copy_Deletion"], col = NA))
  },
  Homozygous_Deletion = function(x, y, w, h) {
    grid.rect(x, y, w, h*0.9,
              gp = gpar(fill = col["Homozygous_Deletion"], col = NA))
  },
  Low_Amplification = function(x, y, w, h) {
    grid.rect(x, y, w, h*0.9,
              gp = gpar(fill = col["Low_Amplification"], col = NA))
  },
  High_Amplification = function(x, y, w, h) {
    grid.rect(x, y, w, h*0.9,
              gp = gpar(fill = col["High_Amplification"], col = NA))
  }
)

column_title = "OncoPrint for TCGA LAML"
heatmap_legend_param = list(title = "Alternations", at = c('Single_Copy_Deletion','Homozygous_Deletion','Low_Amplification','High_Amplification'), 
                            labels = c('Single_Copy_Deletion','Homozygous_Deletion','Low_Amplification','High_Amplification'))

subt <- read.csv("H:/Project/Smart/20230105m6A/step4/result_type/tcga_fpkm_sweep/p/pearson_pam/3.csv", check.names = F, stringsAsFactors = F, header = T, row.names = 1)
head(subt)
subt$TCGA_Subtype <- paste0('C',subt$Cluster)
subt$Tumor_Sample_Barcode <- rownames(subt)

colnames(cnv) <- str_sub(colnames(cnv),1,15)
sample_id <- intersect(subt$Tumor_Sample_Barcode, colnames(cnv))

cnv <- as.data.frame(t(cnv[,sample_id]))
subt <- subt[sample_id,]
rownames(cnv) == rownames(subt)
cnv$Cluster <- subt$TCGA_Subtype
cnv_list <- split(cnv,cnv$Cluster)

cnv_num_list <- lapply(cnv_list, function(x){
  x <- t(x)
  sum <- rowSums(x != '')/ncol(x)
  id <- names(sum[order(sum, decreasing = T)])[1:10]
  return(id)
})

gene <- unique(c(cnv_num_list[[1]],cnv_num_list[[2]],cnv_num_list[[3]]))
gene <- gene[-1]

anno <- subt %>% arrange(TCGA_Subtype)
library(ComplexHeatmap)
dat <- as.data.frame(t(cnv)[gene,])
onco <- oncoPrint(dat,
                  alter_fun = alter_fun, col = col,
                  remove_empty_columns = T, remove_empty_rows = T,
                  column_title = column_title, heatmap_legend_param = heatmap_legend_param,
                  column_order = anno$Tumor_Sample_Barcode,
                  bottom_annotation = clini)

# col_os = colorRamp2(c(0, 4000), c("white", "blue"))
clini <- HeatmapAnnotation(SubType = anno$TCGA_Subtype,
                           OS  = factor(anno$OS))

pdf("tcga2016result1.pdf", width = 10, height = 6) ##导出pdf文件
draw(onco,annotation_legend_side = "bottom")
dev.off()

df <- as.data.frame(t(cnv_list[[3]][,-23131]))
onco3 <- oncoPrint(df[gene,],
                   alter_fun = alter_fun, col = col,
                   remove_empty_columns = T, remove_empty_rows = T,
                   # column_title = column_title,
                   heatmap_legend_param = heatmap_legend_param,
                   # column_order = anno$Tumor_Sample_Barcode,
                   # bottom_annotation = clini,
                   show_row_names = F,
                   row_names_side = "left",
                   pct_side = "right",
                   right_annotation = NULL,
                   show_heatmap_legend = T)
onco3

pdf("C3_result1.pdf", width = 3, height = 5) ##导出pdf文件
draw(onco3)
dev.off()

library(patchwork)
p1 <- onco1+ onco2+onco3 +plot_layout(c(3,2,2))
pdf("CNV_result.pdf", width = 12, height = 4) ##导出pdf文件
draw(p1)
dev.off()