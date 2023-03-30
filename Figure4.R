gsea_results_m6Amarker <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/gsea_results_m6Amarker.Rds")
cellmarker_results_m6Amarker <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/cellmarker_results_m6Amarker.Rds")

gsea_results_m6Amarker <- lapply(gsea_results_m6Amarker, function(x){
  x = subset(x, p.adjust < 0.05)
  return(x)
})

cellmarker_results_m6Amarker <- lapply(cellmarker_results_m6Amarker, function(x){
  x = subset(x, p.adjust < 0.05)
  return(x)
})



merge_list = list()
for (i in names(gsea_results_m6Amarker)) {
  pathway_id = intersect(gsea_results_m6Amarker[[i]]$ID, cellmarker_results_m6Amarker[[i]]$ID)
  gsea_results_m6Amarker[[i]] = gsea_results_m6Amarker[[i]][pathway_id,]
  cellmarker_results_m6Amarker[[i]] = cellmarker_results_m6Amarker[[i]][pathway_id,]
  merge_list[[i]] = merge(gsea_results_m6Amarker[[i]],gsea_results_m6Amarker[[i]],by = 'ID')
}

hall_gsea <- lapply(merge_list, function(x){
  x = x[grep('HALLM', x$ID),]
})


hall_gsea <- lapply(gsea_results_m6Amarker, function(x){
  x = x[grep('HALLM', x$ID),]
})
write.xlsx(hall_gsea, 'm6Amarker_results_gsea.xlsx')


hall_gsea <- lapply(cellmarker_results_m6Amarker, function(x){
  x = x[grep('HALLM', x$ID),]
})
write.xlsx(hall_gsea, 'cellmarker_results_gsea.xlsx')

listInput <- lapply(hall_gsea, function(x) x$ID)
upset(fromList(listInput), order.by = "freq",nsets = 12)

library(VennDiagram)
inter <- get.venn.partitions(listInput)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
inter = inter[inter$..count.. != 0,]
write.table(inter[-c(5, 6)], 'hall_deg_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)


###UpsetR
gsea_results_m6Amarker <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/gsea_results_m6Amarker.Rds")
gsea_results_m6Amarker <- lapply(gsea_results_m6Amarker, function(x){
  x = subset(x, p.adjust < 0.05)
  return(x[c('ID','NES','p.adjust')])
})

hall_gsea <- lapply(gsea_results_m6Amarker, function(x){
  x = x[grep('HALLM', x$ID),]
})

listInput <- lapply(hall_gsea, function(x) x$ID)
upset(fromList(listInput), order.by = "freq",nsets = 12)

pathway_id = unique(Reduce(append, listInput))
gsea_results_m6Amarker <- lapply(gsea_results_m6Amarker, function(x){
  return(x[pathway_id,c('ID','NES','p.adjust')])
})


plot_dat <- do.call(rbind,gsea_results_m6Amarker)
plot_dat$Cluster = str_split(rownames(plot_dat), '\\.', simplify = T)[,1]
plot_dat$ID = str_split(plot_dat$ID, '_', simplify = T, n = 2)[,2]
plot_dat$ID = gsub('_', ' ', plot_dat$ID)
plot_dat = na.omit(plot_dat)
table(plot_dat$Cluster)

mycol = ggsci::pal_lancet()(2)
plot_dat$NES = ifelse(plot_dat$p.adjust < 0.05, plot_dat$NES, 0)
plot_dat$Cluster <- factor(plot_dat$Cluster,
                           levels = c('cDC','CTL','NK','Plasma','Monocyte','B','T',
                                      'LSPC','cDC-like','GMP-like',
                                      'Mono-like','ProMono-like'))

p1 <- ggplot(plot_dat,aes(x=Cluster,y=ID)) #热图绘制
p2 <- p1+scale_fill_gradient2(low = mycol[1], mid = 'white', high = mycol[2],)+
  theme_minimal()+
  geom_point(aes(size=-log10(p.adjust),
                 fill=NES), shape= 21)+
  theme(axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)
p2
ggsave(plot = p2,'gsea_deg_bob.pdf', width = 7, height = 7)



cellmarker_results_m6Amarker <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/cellmarker_results_m6Amarker.Rds")
cellmarker_results_m6Amarker <- lapply(cellmarker_results_m6Amarker, function(x){
  x = subset(x, p.adjust < 0.05)
  return(x[c('ID','NES','p.adjust')])
})

hall_gsea <- lapply(cellmarker_results_m6Amarker, function(x){
  x = x[grep('HALLM', x$ID),]
})

listInput <- lapply(hall_gsea, function(x) x$ID)
upset(fromList(listInput), order.by = "freq",nsets = 12)

pathway_id = unique(Reduce(append, listInput))
cellmarker_results_m6Amarker <- lapply(cellmarker_results_m6Amarker, function(x){
  return(x[pathway_id,c('ID','NES','p.adjust')])
})


plot_dat <- do.call(rbind,cellmarker_results_m6Amarker)
plot_dat$Cluster = str_split(rownames(plot_dat), '\\.', simplify = T)[,1]
plot_dat$ID = str_split(plot_dat$ID, '_', simplify = T, n = 2)[,2]
plot_dat$ID = gsub('_', ' ', plot_dat$ID)
plot_dat = na.omit(plot_dat)
table(plot_dat$Cluster)

mycol = ggsci::pal_lancet()(2)
plot_dat$NES = ifelse(plot_dat$p.adjust < 0.05, plot_dat$NES, 0)
plot_dat$Cluster <- factor(plot_dat$Cluster,
                           levels = c('cDC','CTL','NK','Plasma','Monocyte','B','T',
                                      'LSPC','cDC-like','GMP-like',
                                      'Mono-like','ProMono-like'))

p1 <- ggplot(plot_dat,aes(x=Cluster,y=ID)) #热图绘制
p2 <- p1+scale_fill_gradient2(low = mycol[1], mid = 'white', high = mycol[2],)+
  theme_minimal()+
  geom_point(aes(size=-log10(p.adjust),
                 fill=NES), shape= 21)+
  theme(axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)
p2
ggsave(plot = p2,'gsea_marker_bob.pdf', width = 7, height = 7)

#SSGSEA
library(clusterProfiler)
library(GSVA)
library(org.Hs.eg.db)
library(GSEABase)
library(survival)
library(survminer)

gmtfile ="H:/Project/Resource/msigdb_v7.5.1/msigdb_v7.5.1_GMTs/h.all.v7.5.1.symbols.gmt"
geneset <- getGmt(gmtfile)  

load("H:/Project/Main/20230105m6A/Result/code/data/TCGA_vst.Rdata")
gsva_dat <- gsva(as.matrix(expr), geneset, method = 'gsva')
cox_sig_id <- Coxoutput(mat = as.data.frame(t(gsva_dat)), subt = meta)

# load("H:/Project/Main/20230105m6A/Result/code/data/GSE106291_vst.Rdata")
# gsva_dat <- gsva(as.matrix(expr), geneset, method = 'gsva')
# cox_sig_106 <- Coxoutput(mat = as.data.frame(t(gsva_dat)), subt = meta)
# 
# load("H:/Project/Main/20230105m6A/Result/code/data/GSE146173_vst.Rdata")
# gsva_dat <- gsva(as.matrix(expr), geneset, method = 'gsva')
# cox_sig_146 <- Coxoutput(mat = as.data.frame(t(gsva_dat)), subt = meta)

# Coxoutput <- read.csv('./Coxoutput_TCGA.csv',row.names = 1)
# Coxoutput <- arrange(Coxoutput,pvalue)  %>% #按照p值排序
#   filter(pvalue < 0.05)
# 
# 
# hrtable <- Coxoutput[c(1,2,5,6,4)]
# hrtable$gene <- str_split(hrtable$gene, '_', simplify = T, n = 2)[,2]
# hrtable$gene <- gsub('_', ' ', hrtable$gene)
# tabletext <- cbind(c("Gene",hrtable$gene),
#                    c("HR",format(round(as.numeric(hrtable$HR),3),nsmall = 3)),
#                    c("lower 95%CI",format(round(as.numeric(hrtable$lower),3),nsmall = 3)),
#                    c("upper 95%CI",format(round(as.numeric(hrtable$upper),3),nsmall = 3)),
#                    c("pvalue",formatC(as.numeric(hrtable$pvalue), format = "e", digits = 2)))
# library(forestplot)
# pdf("forestplot of risk table.pdf", width = 12, height = 10)
# forestplot(labeltext=tabletext,
#            mean=c(NA,as.numeric(hrtable$HR)),#log2(HR)
#            lower=c(NA,as.numeric(hrtable$lower)), #log2(95%置信区间下限)
#            upper=c(NA,as.numeric(hrtable$upper)),#log2(95%置信区间上限)
#            graph.pos=6,#图在表中的列位置
#            graphwidth = unit(.25,"npc"),#图在表中的宽度比
#            fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石
#            col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),#box颜色
#            boxsize=0.3,#box大小固定
#            lwd.ci=1,
#            ci.vertices.height = 0.1,ci.vertices=F,#不显示区间
#            zero=1,#zero线横坐标
#            lwd.zero=2,#zero线宽
#            xticks = c(0,1,5,10,15,20,25,30),#横坐标刻度根据需要可随意设置
#            lwd.xaxis=2,
#            xlab=expression("HR"),
#            hrzl_lines=list("1" = gpar(lwd=2, col="black"),#第二行顶部加黑实线
#                            "2" = gpar(lwd=1, col="grey50", lty=2),#第二行顶部加灰色虚线
#                            # "9" = gpar(lwd=1, col="grey50", lty=2),#第九行顶部加灰色虚线
#                            "27" = gpar(lwd=2, col="black")),#最后一行底部加黑线，""中数字为nrow(tabletext) + 1
#            txt_gp=fpTxtGp(label=gpar(cex=1.2),#各种字体大小设置
#                           ticks=gpar(cex=0.85),
#                           xlab=gpar(cex=1),
#                           title=gpar(cex=1.5)),
#            lineheight = unit(.75,"cm"),#固定行高
#            colgap = unit(0.3,"cm"),
#            mar=unit(rep(1.5, times = 4), "cm"),
#            new_page = F
# )
# invisible(dev.off())


####merge
gsea_results_m6Amarker <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/gsea_results_m6Amarker.Rds")
cellmarker_results_m6Amarker <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/cellmarker_results_m6Amarker.Rds")

gsea_results_m6Amarker <- lapply(gsea_results_m6Amarker, function(x){
  x = subset(x, p.adjust < 0.05)
  return(x)
})

cellmarker_results_m6Amarker <- lapply(cellmarker_results_m6Amarker, function(x){
  x = subset(x, p.adjust < 0.05)
  return(x)
})

merge_list = list()
for (i in names(gsea_results_m6Amarker)) {
  pathway_id = intersect(gsea_results_m6Amarker[[i]]$ID, cellmarker_results_m6Amarker[[i]]$ID)
  gsea_results_m6Amarker[[i]] = gsea_results_m6Amarker[[i]][pathway_id,]
  cellmarker_results_m6Amarker[[i]] = cellmarker_results_m6Amarker[[i]][pathway_id,]
  merge_list[[i]] = merge(gsea_results_m6Amarker[[i]],gsea_results_m6Amarker[[i]],by = 'ID')
}


hall_gsea <- lapply(merge_list, function(x){
  x = x[grep('HALLM', x$ID),]
})

listInput <- lapply(hall_gsea, function(x) x$ID)
upset(fromList(listInput), order.by = "freq",nsets = 12)
pathway_id = unique(Reduce(append, listInput))

gsea_results_m6Amarker <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/gsea_results_m6Amarker.Rds")
cellmarker_results_m6Amarker <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/cellmarker_results_m6Amarker.Rds")


gsea_results_m6Amarker <- lapply(gsea_results_m6Amarker, function(x){
  x = subset(x, p.adjust < 0.05)
  x = x[grep('HALLM', x$ID),]
  return(x)
})
listInput <- lapply(gsea_results_m6Amarker, function(x) x$ID)
pathway_id = unique(Reduce(append, listInput))

gsea_results_m6Amarker <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/gsea_results_m6Amarker.Rds")
cellmarker_results_m6Amarker <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/cellmarker_results_m6Amarker.Rds")

gsea_results_m6Amarker <- lapply(gsea_results_m6Amarker, function(x){
  return(x[pathway_id,c('ID','NES','p.adjust')])
})

plot_dat <- do.call(rbind,gsea_results_m6Amarker)
plot_dat$Cluster = str_split(rownames(plot_dat), '\\.', simplify = T)[,1]
plot_dat$ID = str_split(plot_dat$ID, '_', simplify = T, n = 2)[,2]
plot_dat$ID = gsub('_', ' ', plot_dat$ID)
plot_dat = na.omit(plot_dat)
table(plot_dat$Cluster)

mycol = ggsci::pal_lancet()(2)
plot_dat$NES = ifelse(plot_dat$p.adjust < 0.05, plot_dat$NES, 0)
plot_dat$Cluster <- factor(plot_dat$Cluster,
                           levels = c('cDC','CTL','NK','Plasma','Monocyte','B','T',
                                      'LSPC','cDC-like','GMP-like',
                                      'Mono-like','ProMono-like'))

p1 <- ggplot(plot_dat,aes(x=Cluster,y=ID)) #热图绘制
p2 <- p1+scale_fill_gradient2(low = mycol[1], mid = 'white', high = mycol[2],)+
  theme_minimal()+
  geom_point(aes(size=-log10(p.adjust),
                 fill=NES), shape= 21)+
  theme(axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)
p2
ggsave(plot = p2,'merge_gsea_deg_bob.pdf', width = 7, height = 7)


cellmarker_results_m6Amarker <- lapply(cellmarker_results_m6Amarker, function(x){
  return(x[pathway_id,c('ID','NES','p.adjust')])
})


plot_dat <- do.call(rbind,cellmarker_results_m6Amarker)
plot_dat$Cluster = str_split(rownames(plot_dat), '\\.', simplify = T)[,1]
plot_dat$ID = str_split(plot_dat$ID, '_', simplify = T, n = 2)[,2]
plot_dat$ID = gsub('_', ' ', plot_dat$ID)
plot_dat = na.omit(plot_dat)
table(plot_dat$Cluster)

mycol = ggsci::pal_lancet()(2)
plot_dat$NES = ifelse(plot_dat$p.adjust < 0.05 & plot_dat$NES > 0, plot_dat$NES, NA)
plot_dat$p.adjust = ifelse(plot_dat$p.adjust < 0.05 & plot_dat$NES > 0, 1, NA)
plot_dat$Cluster <- factor(plot_dat$Cluster,
                           levels = c('cDC','CTL','NK','Plasma','Monocyte','B','T',
                                      'LSPC','cDC-like','GMP-like',
                                      'Mono-like','ProMono-like'))

p1 <- ggplot(plot_dat,aes(x=Cluster,y=ID)) #热图绘制
p2 <- p1+scale_fill_gradient2(low = mycol[1], mid = 'white', high = mycol[2],)+
  theme_minimal()+
  geom_point(aes(size=-log10(p.adjust),
                 fill=NES), shape= 8)+
  theme(axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)
p2
ggsave(plot = p2,'merge_gsea_marker_bob_label.pdf', width = 7, height = 7)


cellmarker_raw_gsea <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/cellmarker_raw_gsea.Rds")
m6A_raw_gsea <- readRDS("H:/Project/Main/20230105m6A/Result/code/test/m6A_raw_gsea.Rds")

devtools::load_all("H:/Education/myenrichplot-master/myenrichplot-master/")
gseaplot2(cellmarker_raw_gsea$LSPC,'HALLMARK_MYC_TARGETS_V1', pvalue_table = T,
          title = '(LSPC) MYC_TARGETS_V1')
ggsave('LSPC_MYC_TARGETS_V1.pdf', width = 6, height = 5)

gseaplot2(m6A_raw_gsea$LSPC,'HALLMARK_E2F_TARGETS', pvalue_table = T,
          title = '(LSPC) E2F TARGETS', color = 'green')
ggsave('LSPC_E2F_TARGETS_m6A.pdf', width = 6, height = 5)


gseaplot2(cellmarker_raw_gsea$`cDC-like`,'HALLMARK_INTERFERON_GAMMA_RESPONSE', pvalue_table = T,
          title = '(cDC-like) INTERFERON GAMMA RESPONSE')
ggsave('cDC-like_INTERFERON_GAMMA_RESPONSE.pdf', width = 6, height = 5)

gseaplot2(m6A_raw_gsea$`cDC-like`,'HALLMARK_INTERFERON_GAMMA_RESPONSE', pvalue_table = T,
          title = '(cDC-like) INTERFERON GAMMA RESPONSE', color = 'green')
ggsave('cDC-like_INTERFERON_GAMMA_RESPONSE_m6A.pdf', width = 6, height = 5)


sce <- readRDS("H:/Project/Main/20230105m6A/Result/code/data/sce.Rds")
library(ggpubr)
dat <- data.frame(State = sce$CellType_new, m6A_modification_activity = sce$m6A_gene)
dat <- subset(dat, State %in% c('LSPC-Cycle','LSPC-Primed','LSPC-Quiescent'))
ggboxplot(dat, x = 'State', y = 'm6A_modification_activity', fill = 'State',
          palette = 'lancet')+
  stat_compare_means()+
  theme_bw()+
  theme(legend.position = 'none')
ggsave('LSPC_cell_cycle.pdf', width = 7, height = 5)
