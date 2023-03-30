#Figure1
library(Seurat)
library(tidyverse)
#Figure1A
sce <- readRDS("H:/Project/Main/20230105m6A/Result/code/data/sce.Rds")
m6A_gene <- readRDS("H:/Project/Main/20230105m6A/Result/code/data/m6A_gene.Rds")
# VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
DimPlot(sce, group.by = 'CellType_int', reduction = 'tsne', label = T)
ggsave('Dimplot_CellType_new.pdf', height = 4.5,width = 6)


# sce <- NormalizeData(sce)
# sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce,features = m6A_gene)
dat <- sce@assays$RNA@scale.data
# dat <- DoHeatmap(sce, features = m6A_gene, group.by = 'CellType_new', label = T)

# library(pheatmap)
# dat[dat > 2] = 2
# ann_col <- data.frame(row.names = rownames(sce@meta.data),Celltype = sce@meta.data['CellType_new'])
# ann_row <- read.csv('../../../m6A/m6a_genesets1.csv',row.names = 1)
# dat <- dat[rownames(ann_row),]
# dat <- dat[,order(ann_col$CellType_new)]
# pheatmap::pheatmap(dat,
#                    show_colnames = F, cluster_rows = F, cluster_cols = F,
#                    annotation_row = ann_row, annotation_col = ann_col)
library(ComplexHeatmap)
ann_col <- data.frame(row.names = rownames(sce@meta.data),Celltype = sce@meta.data['CellType_int'])
dat[dat > 2] = 2
dat[dat < -2] = -2
p1 <- Heatmap(
  dat, name = "expression", 
  col = c('grey','white','red'),
  column_split = factor(ann_col$CellType_int, levels = c('LSPC','ProMono-like','GMP-like',
                                                         'Mono-like','cDC-like','cDC','NK','CTL','T','Plasma','B','Monocyte')),
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  show_column_names = F,
  cluster_rows = F,
  column_title_gp = gpar(
    fill = rainbow(14)[1:514]
  ),
  clustering_distance_columns = 'euclidean',
  
)
pdf('1BHeatmap.pdf',width = 9, height = 5)
draw(p1)
dev.off()


###AUCell
library(Seurat)
library(tidyverse)
library(AUCell)
sce <- readRDS("H:/Project/Main/20230105m6A/step1/sce.Rds")
m6A_gene <- readRDS("H:/Project/Main/20230105m6A/m6A_gene.Rds")

dat <- sce@assays$RNA@counts

cells_rankings <- AUCell_buildRankings(dat,plotStats = F, splitByBlocks=TRUE)  # 关键一步
cells_AUC <- AUCell_calcAUC(as.character(m6A_gene), cells_rankings, aucMaxRank=ceiling(nrow(cells_rankings)*0.03))

cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE,assign=TRUE,nBreaks = 100)#
# ggsave('m6A_AUC_Thresholds.pdf',height = 5,width = 7)

cells_assignment$geneSet$aucThr$thresholds
cells_assignment$geneSet$aucThr$selected


aucs <- t(getAUC(cells_AUC))
aucs_bi <- ifelse(aucs > cells_assignment$geneSet$aucThr$selected, 1,0)

table(colnames(sce) == rownames(aucs))
sce$m6A_gene <- aucs[,1]
sce$m6A_bi <- aucs_bi[,1]

FeaturePlot(sce,features = 'm6A_gene',reduction = 'tsne', min.cutoff = 0, cols = c('grey','red'))+ggtitle('m6A geneset AUC Score')
ggsave('m6A_AUC_tSNE.pdf',width = 6,height = 5)
FeaturePlot(sce,features = 'm6A_bi',reduction = 'tsne', cols = c('grey','red'),pt.size = 0.8)+ggtitle('Pyroptosis AUC Score')
ggsave('m6A_AUCbi_tSNE.pdf',width = 11,height = 10)


library(ggpubr)
library(ggsci)
dat <- sce@meta.data[,c('m6A_gene','CellType_int')]
cg <- aggregate(dat$m6A_gene, by=list(type=dat$CellType_int),median)
dat$CellType_int <- factor(dat$CellType_int,levels = as.character(cg[order(cg$x,decreasing = F),'type']))

ggviolin(dat, x = "CellType_int", y = "m6A_gene",fill = 'CellType_int',
         palette = colorRampPalette (c ("white", "red"))(14),
         add = 'boxplot',
         ylab = 'Cell AUC Score', xlab = '',title = 'AUC Score', rotate = TRUE, flip = T)+
  theme(legend.position = 'none')+
  geom_hline(aes(yintercept=0.048),linetype=6,col="red",size=1)
ggsave('Celltype_AUC_boxplot.pdf',width = 4,height = 6)


###活跃或非活跃细胞比例图
library(ggstatsplot)
dat <- sce@meta.data[c('CellType_int','m6A_bi','m6A_gene')]
cg <- aggregate(dat$m6A_gene, by=list(type=dat$CellType_int),mean)
dat$CellType_int <- factor(dat$CellType_int,levels = as.character(cg[order(cg$x,decreasing = T),'type']))
# dat$m6A_bi <- factor(dat$m6A_bi,levels = c(0,1),labels = c('Inactive','Active'))

colnames(dat)[2] <- 'ActiveState'
p <- ggbarstats(dat,'ActiveState','CellType_int')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(aes(yintercept=0.25),linetype=6,col="red",size=1)
p[["data"]][[".label"]] <- NA
p

dat$CellType_int <- factor(dat$CellType_int,levels = c('ProMono-like','GMP-like','LSPC',
                                                       'cDC','cDC-like','NK','CTL','Mono-like','T','Plasma','B','Monocyte'))
p <- ggbarstats(dat,'ActiveState','CellType_int')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('Celltype_Activate_prop.pdf',width = 9,height = 7)

#将整体细胞划分为两类：活跃和非活跃然后进行差异分析和通路富集
#GOKEGG Pval < 0.05
sce@meta.data[["m6A_bi"]] <- factor(sce@meta.data[["m6A_bi"]],levels = c(0,1),labels = c('Inactive','Active'))
sce_markers <- FindMarkers(sce,ident.1 = 'Inactive',ident.2 = 'Active',group.by = 'm6A_bi')
sce_markers <- sce_markers[sce_markers$p_val_adj < 0.05,]

#通路富集
#全部通路
if(T){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  sce_markers$gene <- rownames(sce_markers)
  df <- bitr(unique(sce_markers$gene), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  markers=merge(sce_markers,df,by.y='SYMBOL',by.x='gene')
  gene_diff= as.character(markers[ ,'ENTREZID'] ) 
  gene_all= rownames(sce)
  df <- bitr(unique(gene_all), fromType = "SYMBOL",
             toType = c( "ENTREZID"),
             OrgDb = org.Hs.eg.db)
  gene_all= df$ENTREZID
  
  R.utils::setOption( "clusterProfiler.download.method",'auto' )
  kk.result <- enrichKEGG(gene         = gene_diff,
                          organism     = 'hsa',
                          universe     = gene_all,
                          pvalueCutoff = 0.9,
                          qvalueCutoff =0.9)
  kk.result <- setReadable(kk.result, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  kk.result@result$'Enrichment Fold'=apply(kk.result@result,1,function(x){
    GeneRatio=eval(parse(text=x["GeneRatio"]))
    BgRatio=eval(parse(text=x["BgRatio"]))
    enrichment_fold=round(GeneRatio/BgRatio,2)
    return(enrichment_fold)
  })
  
  
  go_enrich_results <- lapply( c('BP','MF','CC') , function(ont) {
    ego <- enrichGO(gene          = gene_diff,
                    universe      = gene_all,
                    OrgDb         = org.Hs.eg.db,
                    ont           = ont ,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.99,
                    qvalueCutoff  = 0.99,
                    readable      = TRUE)
    
    ego@result$'Enrichment Fold'=apply(ego@result,1,function(x){
      GeneRatio=eval(parse(text=x["GeneRatio"]))
      BgRatio=eval(parse(text=x["BgRatio"]))
      enrichment_fold=round(GeneRatio/BgRatio,2)
      return(enrichment_fold)
    })
    #print( head(ego) )
    return(ego)
  })
  
  
  #富集结束，文本结果整理
  dat_enrich <- list(GOBP = go_enrich_results[[1]],
                     GOMF = go_enrich_results[[2]],
                     GOCC = go_enrich_results[[3]],
                     KEGG = kk.result)
  
  #文本结果下载
  library(openxlsx)
  write.xlsx(dat_enrich,paste0('m6A_enrich_results.xlsx'))
  
  
  #统一作图
  enrich_plot_enrich <- list()
  for ( i in names(dat_enrich)) {
    dat <- dat_enrich[[i]]@result %>% slice_min(pvalue,n=8)
    
    enrich_plot_enrich[[i]] <- ggplot(dat,aes(`Enrichment Fold`, fct_reorder(Description, `Enrichment Fold`)))+
      geom_segment(aes(xend=0, yend = Description))+
      geom_point(aes(color=p.adjust, size = Count)) +
      scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
      scale_size_continuous(range=c(2, 10)) +
      theme_classic(base_size = 16) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 40))+
      ylab(NULL)+ggtitle(paste(i))
  }
  
  
  p3 <- enrich_plot_enrich[[1]]+enrich_plot_enrich[[2]]+enrich_plot_enrich[[3]]+enrich_plot_enrich[[4]]
  ggsave(plot = p3,filename = paste0('m6A_enrich_results.pdf'),width = 18,height = 10)
  
  save(go_enrich_results,kk.result,file = paste0('m6A_enrich_results.Rdata'))

}
#使用GO语义简化富集结果
if(T){
  R.utils::setOption( "clusterProfiler.download.method",'auto' )
  kk.result <- enrichKEGG(gene         = gene_diff,
                          organism     = 'hsa',
                          universe     = gene_all,
                          pvalueCutoff = 0.05,
                          qvalueCutoff =0.05,)
  kk.result <- setReadable(kk.result, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  kk.result@result$'Enrichment Fold'=apply(kk.result@result,1,function(x){
    GeneRatio=eval(parse(text=x["GeneRatio"]))
    BgRatio=eval(parse(text=x["BgRatio"]))
    enrichment_fold=round(GeneRatio/BgRatio,2)
    return(enrichment_fold)
  })
  
  
  go_enrich_results <- lapply( c('BP','MF','CC') , function(ont) {
    ego <- enrichGO(gene          = gene_diff,
                    universe      = gene_all,
                    OrgDb         = org.Hs.eg.db,
                    ont           = ont ,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
    
    ego@result$'Enrichment Fold'=apply(ego@result,1,function(x){
      GeneRatio=eval(parse(text=x["GeneRatio"]))
      BgRatio=eval(parse(text=x["BgRatio"]))
      enrichment_fold=round(GeneRatio/BgRatio,2)
      return(enrichment_fold)
    })
    #print( head(ego) )
    return(ego)
  })
  
  
  #富集结束，文本结果整理
  dat_enrich <- list(GOBP = clusterProfiler::simplify(go_enrich_results[[1]],cutoff=0.7,by="p.adjust",select_fun=min,measure="Wang"),
                     GOMF = clusterProfiler::simplify(go_enrich_results[[2]],cutoff=0.7,by="p.adjust",select_fun=min,measure="Wang"),
                     GOCC = clusterProfiler::simplify(go_enrich_results[[3]],cutoff=0.7,by="p.adjust",select_fun=min,measure="Wang"),
                     KEGG = kk.result)
  
  #文本结果下载
  openxlsx::write.xlsx(dat_enrich, file = paste0('m6A_smell_enrich_results.xlsx'))
  
  
  #统一作图
  enrich_plot_enrich <- list()
  for ( i in names(dat_enrich)) {
    dat <- dat_enrich[[i]]@result %>% slice_min(p.adjust,n=10)
    
    enrich_plot_enrich[[i]] <- ggplot(dat,aes(`Enrichment Fold`, fct_reorder(Description, `Enrichment Fold`)))+
      geom_segment(aes(xend=0, yend = Description))+
      geom_point(aes(color=p.adjust, size = Count)) +
      scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
      scale_size_continuous(range=c(2, 10)) +
      theme_classic(base_size = 16) +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 40))+
      ylab(NULL)+ggtitle(paste(i))
  }
  
  #展示上调基因富集图
  p3 <- enrich_plot_enrich[[1]]+enrich_plot_enrich[[2]]+enrich_plot_enrich[[3]]+enrich_plot_enrich[[4]]
  ggsave(plot = p3,filename = paste0('m6A_smell_Enrich_plot.pdf'),width = 18,height = 10)
  saveRDS(dat_enrich, file = paste0('m6A_smell_dat_enrich.Rds'))
}


###对全部的50个癌症通路进行GSVA
#GSVA
library(msigdbr)
library(GSVA)
## 表达矩阵
expr=as.matrix(sce@assays$RNA@data)
##通路基因集
msgdH = msigdbr(species = "Homo sapiens", category = "H")
HSet = msgdH %>% split(x = .$gene_symbol, f = .$gs_name)
names(HSet) <- str_split(names(HSet),'_',simplify = T, n = 2)[,2]
names(HSet) <- gsub('_',' ',names(HSet))

h <- gsva(expr, gset.idx.list = HSet, kcdf="Gaussian",method = "zscore")

library(limma)
## limma gsva通路活性评估
de_gsva <- function(exprSet,meta,compare = NULL){
  
  
  allDiff = list()
  design <- model.matrix(~0+factor(meta))
  colnames(design)=levels(factor(meta))
  rownames(design)=colnames(exprSet)
  
  fit <- lmFit(exprSet,design)
  if(length(unique(meta))==2){
    if(is.null(compare)){
      stop("there are 2 Groups,Please set  compare value")
    }
    contrast.matrix<-makeContrasts(contrasts = compare,levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2,adjust='fdr', coef=1, number=Inf)
    allDiff[[compare]] = na.omit(tempOutput)
    
  }else if(length(unique(meta))>2){
    for(g in colnames(design)){
      fm = ""
      for(gother in colnames(design)[which(!colnames(design) %in% g)]){
        fm = paste0(fm,"+",gother)
      } 
      
      fm = paste0(g,"VsOthers = ",g,"-(",substring(fm,2),")/",ncol(design)-1)
      contrast.matrix <- makeContrasts(contrasts = fm,levels=design)
      fit2 <- contrasts.fit(fit, contrast.matrix) 
      fit2 <- eBayes(fit2)
      allDiff[[g]]=topTable(fit2,adjust='fdr',coef=1,number=Inf)
    }
  }else{
    stop("error only have one group")
  }
  
  return(allDiff)
}
meta <- sce@meta.data[,c("m6A_bi")]
Diff =de_gsva(exprSet = h ,meta = meta,compare = "Active-Inactive")

idiff <-Diff[["Active-Inactive"]]
df <- data.frame(ID = rownames(idiff), score = idiff$t )
Padj_threshold = 0.01
df$group =sapply(1:nrow(idiff),function(x){
  if(idiff[x,"logFC"]>0 & idiff[x,"adj.P.Val"]<Padj_threshold){return("up")}
  else if(idiff[x,"logFC"]<0 & idiff[x,"adj.P.Val"]<Padj_threshold){return("down")}
  else{return("noSig")}
})

# 按照score排序
df$hjust = ifelse(df$score>0,1,0)
df$nudge_y = ifelse(df$score>0,-0.1,0.1)
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
limt = max(abs(df$score))
ggplot(sortdf, aes(ID, score,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("down","noSig","up"),
                    values = c("#008020","grey","#08519C"))+
  geom_text(data = df, aes(label = df$ID, y = df$nudge_y),
            nudge_x =0,nudge_y =0,hjust =df$hjust,
            size = 2)+
  labs(x = " pathways",
       y="t value of GSVA score")+
  scale_y_continuous(limits=c(-limt,limt))+
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6)
        #panel.border = element_blank()
  )+
  theme(plot.title = element_text(hjust = 0.5,size = 12),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 12),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = limt
  )
ggsave('m6A_GSVA.pdf', width = 5, height = 7)
save(h,Diff,file = 'GSVA_result.Rdata')

