

n <- t(scale(t(dat[all,])))
n[n > 2] = 2
n[n < -2] = -2
col_anno <- data.frame(row.names = rownames(meta), TCGA_Subtype = meta$TCGA_Subtype)
pheatmap::pheatmap(n, cluster_rows = F, show_colnames = T,
                   cluster_cols = T, annotation_col = meta)

meta <- read.table('clipboard', header = T)
rownames(meta) <- meta$ID

table(meta$State)


library(limma)
group_list <- meta$State
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
head(design)
exprSet=as.data.frame(dat)
rownames(design)=colnames(exprSet)
design
contrast.matrix<-makeContrasts("CR-NR",
                               levels = design)
contrast.matrix ##这个矩阵声明，我们要把 Tumor 组跟 Normal 进行差异分析比较

deg = function(exprSet,design,contrast.matrix){
  ##step1
  fit <- lmFit(exprSet,design)
  ##step2
  fit2 <- contrasts.fit(fit, contrast.matrix) 
  ##这一步很重要，大家可以自行看看效果
  
  fit2 <- eBayes(fit2)  ## default no trend !!!
  ##eBayes() with trend=TRUE
  ##step3
  tempOutput = topTable(fit2, coef=1, n=Inf)
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)
}

res_deg = deg(exprSet,design,contrast.matrix)

head(res_deg)

save(res_deg,file = 'imm_deg.Rdata')

table(res_deg$P.Value < 0.01)

logfc <- quantile(abs(gene_diff$logFC), 0.75)
gene_diff <- res_deg
#首先这里自定义根据 |log2FC| >= 1 和 P.Value < 0.05 标记差异类型
#新生sig行标记up down none
gene_diff[which(gene_diff$P.Value < 0.05 & gene_diff$logFC <= -logfc),'sig'] <- 'Down'
gene_diff[which(gene_diff$P.Value < 0.05 & gene_diff$logFC >= logfc),'sig'] <- 'Up'
gene_diff[which(gene_diff$P.Value >= 0.05 | abs(gene_diff$logFC) < logfc),'sig'] <- 'None'
sub = paste0('Up:',length(which((gene_diff[,7]=='Up'))),' ', 'Down:',length(which((gene_diff[,7]=='Down'))))

#横轴 log2FC，纵轴 -log10(P.Value)，颜色表示差异

p <- ggplot(gene_diff, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.6, size = 1) + #alpha调节透明度，size是点的大小
  scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('Up', 'Down', 'None')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
  theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
  geom_vline(xintercept = c(-logfc, logfc), color = 'gray', size = 0.3) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.3) +
  #
  labs(x = 'Log2 Fold Change', y = '-Log10(P.Value)', color = '', title = 'CR vs NR',subtitle = sub)
p

ggsave(plot = p, filename = 'CR_NR_volcano.pdf', width = 4, height = 5)

save(dat,meta,file = 'GSE178926_subt.Rdata')

load("H:/Project/Main/20230105m6A/two/imm/GSE178926_subt.Rdata")
library(tidyverse)
library(glmnet)
library(caret)
gene <- subset(res_deg, P.Value < 0.05 & abs(logFC) > 1) %>% rownames()
train <- as.data.frame(t(dat))
standard <- preProcess(train, method = c("center","scale"))
train_std <- predict(standard, train)
#后面svmRFE函数要求group的类型为factor



train_std[1:4,1:4]
###lasso
# 转为lasso需要的格式

rownames(train_std) == rownames(meta)
train_std$group <- meta$State
train_std$group <- ifelse(train_std$group == 'NR',0,1)

# (y <- ifelse(train_std$group == "0", 0,1)) #把分组信息换成01
# x <- as.matrix(train_std[,gene])
#library(glmnet)
library(car) #package to calculate Variance Inflation Factor
library(corrplot) #correlation plots
library(leaps) #best subsets regression
library(glmnet) #allows ridge regression, LASSO and elastic net
library(caret) #this will help identify the appropriate parameters


grid <- expand.grid(.alpha = seq(0,1, by=.2), 
                    .lambda = seq(0.00, 0.2, by = 0.02))
table(grid)
head(grid)
control <- trainControl(method = "LOOCV") #selectionFunction="best"
set.seed(701) #our random seed
train <- train_std[,c('group',gene)]
rownames(train) == rownames(meta)
train$group<- as.factor(train$group)
enet.train = train(group ~ ., data = train, 
                   method = "glmnet", 
                   trControl = control, 
                   tuneGrid = grid)
enet.train

x <- as.matrix(train[2:48])
y <- as.numeric(train$group)

enet <- glmnet(x, y,family = "binomial", 
               alpha = 0.6, 
               lambda = 0.06)
enet.coef <- coef(enet, s = 0.06, exact = TRUE)
enet.coef
newx <- x
test <- train
enet.y <- predict(enet, newx = newx, type = "response",  s= 0.06)
plot(enet.y, test$group, xlab = "Predicted", 
     ylab = "Actual", main = "Elastic Net")
enet.resid <- enet.y - as.numeric(test$group)
mean(enet.resid^2)

library("pROC")
#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2","#990000")

pdf("GSE178926_ROC.pdf",height=5,width=5)
auc.out <- c()

#先画第一条线，此处是miRNA1
x <- plot.roc(as.numeric(test$group),enet.y[,1],ylim=c(0,1),xlim=c(1,0),
              smooth=F, #绘制平滑曲线
              ci=TRUE, 
              main="GSE178926 Cohort",
              #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
              col=mycol[2],#线的颜色
              lwd=2, #线的粗细
              legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1

ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限

auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
auc.out <- rbind(auc.out,auc.ci)


legend.name <- paste("AUC:",auc.out[1],auc.out[2],sep=" ")
legend("bottomright", 
       legend=legend.name,
       col = mycol[2],
       lwd = 2,
       bty="n")
dev.off()

col_anno <- data.frame(row.names = rownames(meta), Response = meta$State, Response_Score = enet.y[,1])
col_anno$Response_Score <- as.numeric(col_anno$Response_Score)
col_anno <- col_anno %>% dplyr::arrange(Response_Score)
n = t(train[order(col_anno$Response),2:48])
n[n < -2] = -2
n[n > 2] = 2
mycol = ggsci::pal_lancet()(2)
pheatmap::pheatmap(n[,rownames(col_anno)], show_colnames = F,
                   annotation_col = col_anno,
                   cluster_cols = F,
                   color = colorRampPalette(colors = c(mycol[1],'white',mycol[2]))(100),
                   filename = 'GSE178926_pheatmap_1.pdf',
                   width = 8, height = 6)

library(ggpubr)

coef_mat <- as.data.frame(as.matrix(enet.coef))
coef_mat$Gene <- rownames(coef_mat)
coef_mat <- coef_mat[coef_mat$s1 != 0,]
colnames(coef_mat) <- c('Coefficient','Feature')
coef_mat$sig <- ifelse(coef_mat$Coefficient > 0, 'Up', 'Down') %>%  as.factor()
coef_mat <- coef_mat %>% dplyr::arrange(Coefficient)
coef_mat <- coef_mat[-1,]
ggbarplot(coef_mat, x = 'Feature', y = 'Coefficient', fill = 'sig',
          palette = 'lancet',
          )+coord_flip()+
  theme_bw()+
  ggtitle('GSE178926 Elastic Net',subtitle = 'alpha = 0.6, lambda = 0.06')+
  theme(legend.position = 'none')
ggsave('GSE178926_ElasticNet_Coef.pdf', width = 4, height = 5)


load("H:/Project/Main/20230105m6A/two/subtype/data/TCGA_subt.Rdata")
x <- expr[gene,]
x[is.na(x)] = 0
rownames(x) = gene
x <- t(x)
standard <- preProcess(x, method = c("center","scale"))
newx <- predict(standard, x)
enet.y <- predict(enet, newx = newx, type = "response",  s= 0.06)

rownames(enet.y) == rownames(meta)
meta$Response_Score <- enet.y[,1]
ggboxplot(meta, x = 'TCGA_Subtype', y = 'Response_Score', fill = 'TCGA_Subtype',
          palette = 'lancet', title = 'TCGA AZA+Pembro Response')+
  theme_bw()+
  theme(legend.position = 'none')
ggsave('TCGA_Response.pdf', width = 5, height = 5)


load("H:/Project/Main/20230105m6A/two/subtype/data/GSE146173_subt.Rdata")
load("H:/Project/Main/20230105m6A/two/subtype/data/GSE146173_vst.Rdata")

x <- GSE146173_fpkm[gene,rownames(subt)]
x[is.na(x)] = 0
rownames(x) = gene
x <- t(x)
standard <- preProcess(x, method = c("center","scale"))
newx <- predict(standard, x)
enet.y <- predict(enet, newx = newx, type = "response",  s= 0.06)

rownames(enet.y) == rownames(subt)
subt$Response_Score <- enet.y[,1]
subt$TCGA_Subtype <- paste0('C',subt$res_cluster)
ggboxplot(subt, x = 'TCGA_Subtype', y = 'Response_Score', fill = 'TCGA_Subtype',
          palette = 'lancet', title = 'GSE146173 AZA+Pembro Response')+
  theme_bw()+
  theme(legend.position = 'none')
ggsave('GSE146173_Response.pdf', width = 5, height = 5)


load("H:/Project/Main/20230105m6A/two/subtype/data/GSE106291_subt.Rdata")
load("H:/Project/Main/20230105m6A/two/subtype/data/GSE106291_vst.Rdata")

x <- GSE106291_fpkm[gene,rownames(subt)]
x[is.na(x)] = 0
rownames(x) = gene
x <- t(x)
standard <- preProcess(x, method = c("center","scale"))
newx <- predict(standard, x)
enet.y <- predict(enet, newx = newx, type = "response",  s= 0.06)

rownames(enet.y) == rownames(subt)
subt$Response_Score <- enet.y[,1]
subt$TCGA_Subtype <- paste0('C',subt$res_cluster)
ggboxplot(subt, x = 'TCGA_Subtype', y = 'Response_Score', fill = 'TCGA_Subtype',
          palette = 'lancet', title = 'GSE106291 AZA+Pembro Response')+
  theme_bw()+
  theme(legend.position = 'none')
ggsave('GSE106291_Response.pdf', width = 5, height = 5)
