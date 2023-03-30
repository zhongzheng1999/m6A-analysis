
devtools::load_all("E:/pRRophetic_0.5/pRRophetic")
library(ggplot2)
library(cowplot)
library(ggpubr)
compare_list <- list(c('C1','C2'),
                     c('C1','C3'),
                     c('C2','C3'))

load("H:/Project/Main/20230105m6A/two/subtype/data/GSE146173_subt.Rdata")
load("H:/Project/Main/20230105m6A/two/subtype/data/GSE146173_vst.Rdata")

expr <- GSE146173_fpkm[,rownames(subt)]
colnames(expr) == rownames(subt)

dat <- expr
dat[1:3, 1:3]

ann <- subt
ann$ImmClust <- paste0("C",ann$res_cluster)
head(ann)

table(ann$ImmClust)

#药物名字
GCP.drug <- read.table("drug.txt") #如果要例文的两种药物，就换成drug_eg.txt
GCP.drug <- GCP.drug$V1
#这里以前2种药物为例
# GCP.drug <- GCP.drug[1:2]

#确定是否可行
drug_dat <- data.frame(drug = GCP.drug, type = NA, row.names = GCP.drug)
for (i in GCP.drug) {
  a <- getCGPinfo(i,tissueType = 'blood')
  drug_dat[i,'type'] <- class(a$trainDataOrd)[1]
}

GCP.drug <- rownames(drug_dat[drug_dat$type == 'matrix',])


# 自定义足够多的box的颜色，颜色数量至少等于分组数量
jco <- ggsci::pal_lancet()(3)

### 药物敏感性预测 ###、
GCPinfo <- GCP.IC50 <- GCP.expr <- cvOut <- predictedPtype <- predictedBoxdat <- list() # 初始化列表
plotp <- list()




for (drug in GCP.drug) {
  set.seed(1248103) # 因为预测过程默认10-fold CV，所以设置种子以便结果可重复
  cat(drug," starts!\n") # 提示当前药物已开始分析
  
  # 预测IC50值，采用默认参数，详细可参考??pRRopheticPredict参数列表
  predictedPtype[[drug]] <- pRRopheticPredict(testMatrix = as.matrix(dat[,rownames(ann)]),
                                              drug = drug,
                                              tissueType = "blood",
                                              selection = 1) # 1表示若有重复基因取均值处理
  
  if(!all(names(predictedPtype[[drug]])==rownames(ann))) {stop("Name mismatched!\n")} # 若名字不匹配则报错退出
  
  predictedBoxdat[[drug]] <- data.frame("est.ic50"=predictedPtype[[drug]],
                                        "ImmClust"=ann$ImmClust, # 这里我修改了C1和C2的名字
                                        row.names = names(predictedPtype[[drug]])) 
  predictedBoxdat[[drug]]$ImmClust <- factor(predictedBoxdat[[drug]]$ImmClust,levels = c("C1","C2",'C3'),ordered = T) # 把类改成因子变量
  
  # 绘图
  p <- ggboxplot(data = predictedBoxdat[[drug]], x='ImmClust', y='est.ic50',
                 fill = 'ImmClust', palette = 'lancet',
                 title = drug, xlab = "", ylab = "Estimated IC50")+
    theme_bw()+
    theme(legend.position="none") + # 倾斜字体
    # theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) +
    stat_compare_means(comparisons = compare_list)+
    stat_compare_means(label.y = max(predictedBoxdat[[drug]]$est.ic50)+1)
  # p <- ggplot(data = predictedBoxdat[[drug]], aes(x=ImmClust, y=est.ic50))
  # p <- p + geom_boxplot(aes(fill = ImmClust)) + 
  #   scale_fill_manual(values = jco[1:length(unique(ann$ImmClust))]) + #自定义box的配色
  #   theme(legend.position="none") + # 倾斜字体
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),plot.title = element_text(size = 12, hjust = 0.5)) + 
  #   xlab("") + ylab("Estimated IC50") + 
  #   ggtitle(drug) # 补上title
  
  plotp[[drug]] <- p # 保存在列表里供合并图片用
  cat(drug," has been finished!\n") # 提示当前药物已分析结束
}

save(plotp,predictedBoxdat,file = 'GSE146173_predictedBoxdat.Rdata')
dat_drug <- do.call(rbind, predictedBoxdat)
dat_drug$drug <- str_split(rownames(dat_drug), '\\.T', simplify = T)[,1]

compare_dat <- compare_means(est.ic50~ImmClust,data = dat_drug, group.by = 'drug')
write.csv(compare_dat, 'GSE146173_compare_dat.csv')

# compare_dat_sub <- compare_dat[compare_dat$p < 0.05, ]
# drug_id <- names(which(table(compare_dat_sub$drug) > 2))
drug_id <- c('ABT.263','BI.D1870','BIBW2992','Bleomycin','Pazopanib')
# 合并图片
#适合展示两种药物
# p1 <- plot_grid(plotp[[1]],plotp[[2]],labels = c("A","B"),nrow = 1) # title可以AI下拉到合适位置，就如例文所示
# p1
# ggsave("GSE146173_boxplot of predicted IC50.pdf", width = 6, height = 5)

# 适合展示多种药物
p2 <- plot_grid(plotlist=plotp[drug_id], ncol=5)
p2
ggsave("GSE146173_boxplot of predicted IC50_multiple_1.pdf", width = 18, height = 5)
