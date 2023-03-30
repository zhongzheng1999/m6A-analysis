

library(tidyverse)
library(patchwork)
load("H:/Project/Main/20230105m6A/two/subtype/data/GSE106291_vst.Rdata")
load("H:/Project/Main/20230105m6A/two/workflow/TCGA_sweep_vst_surv/all/pearson_pam/3.csv_GSE106291_cluster.Rdata")

library(survival)
library(survminer)
for (d in c('rf_pred','nb_pred','svm_pred','knn_pred','bp_pred1')) {
  anno2 <- data.frame(Cluster = get(d))
  colnames(anno2) <- d
  clinical_annot <- cbind(clinical_annot,anno2)
}

clinical_annot$res_cluster <- apply(clinical_annot[5:9], 1, function(x){
  x <- ifelse(sum(x == 1) > 2,1,
              ifelse(sum(x == 2) > 2,2,
                     ifelse(sum(x == 3) > 2,3,NA)
  ))
  return(x)
})
table(clinical_annot$res_cluster)
subt <- na.omit(clinical_annot)
# GSE106291_fpkm <- GSE106291_fpkm[rownames(subt),]
# save(GSE106291_fpkm,subt,file = 'GSE106291_subt.Rdata')
surv_object = Surv(subt$OS_Time, subt$OS_Status)
fit1 <- survfit(surv_object ~ res_cluster, data = subt)
summary(fit1)
p <- ggsurvplot(fit1,palette = 'lancet',
                risk.table =TRUE,pval =TRUE,
                conf.int =F,xlab ="Time in Days",
                ggtheme =theme_light(),
                ncensor.plot = F)
p
p1 <- p$plot/p$table + plot_layout(heights = c(3, 1))
ggsave(plot = p1, filename = 'GSE106291_survive.pdf', height=7, width = 8)

tide <- read.csv("H:/Project/Main/20230105m6A/two/imm/data/all_mean_TIDE.csv",
                 row.names = 1, header = T)
tide <- tide[grep(pattern = 'GSE106291', row.names(tide), value = T),]
rownames(tide) <- str_split(rownames(tide),'\\.', simplify = T, n = 2)[,2]

tide$Responder <- ifelse(str_detect(tide$Responder,"False"),"NR","R")
tide <- tide[rownames(subt),]
tide$Cluster <- subt$res_cluster
library(ggstatsplot)
p2 <- ggbarstats(tide, y = 'Cluster', x = 'Responder', )+
  scale_fill_manual(values = c('R'= '#FF1493','NR'='#40E0D0'))
p2
ggsave(plot = p2, filename = 'GSE106291_barstat.pdf', height=7, width = 6)



load("H:/Project/Main/20230105m6A/two/subtype/data/GSE146173_vst.Rdata")
load("H:/Project/Main/20230105m6A/two/workflow/TCGA_sweep_vst_surv/all/pearson_pam/3.csv_GSE146173_cluster.Rdata")

library(survival)
library(survminer)
for (d in c('rf_pred','nb_pred','svm_pred','knn_pred','bp_pred1')) {
  anno2 <- data.frame(Cluster = get(d))
  colnames(anno2) <- d
  clinical_annot <- cbind(clinical_annot,anno2)
}

clinical_annot$res_cluster <- apply(clinical_annot[5:9], 1, function(x){
  x <- ifelse(sum(x == 1) > 2,1,
              ifelse(sum(x == 2) > 2,2,
                     ifelse(sum(x == 3) > 2,3,NA)
              ))
  return(x)
})
table(clinical_annot$res_cluster)

subt <- na.omit(clinical_annot)

# GSE146173_fpkm <- GSE146173_fpkm[rownames(subt),]
# save(GSE146173_fpkm,subt,file = 'GSE146173_subt.Rdata')

surv_object = Surv(subt$OS_Time, subt$OS_Status)
fit1 <- survfit(surv_object ~ res_cluster, data = subt)
summary(fit1)
p <- ggsurvplot(fit1,palette = 'lancet',
                risk.table =TRUE,pval =TRUE,
                conf.int =F,xlab ="Time in Days",
                ggtheme =theme_light(),
                ncensor.plot = F)
p
p1 <- p$plot/p$table + plot_layout(heights = c(3, 1))
ggsave(plot = p1, filename = 'GSE146173_survive.pdf', height=7, width = 8)


tide <- read.csv("H:/Project/Main/20230105m6A/two/imm/data/all_mean_TIDE.csv",
                 row.names = 1, header = T)
tide <- tide[grep(pattern = 'GSE146173', row.names(tide), value = T),]
rownames(tide) <- str_split(rownames(tide),'\\.', simplify = T, n = 2)[,2]

tide$Responder <- ifelse(str_detect(tide$Responder,"False"),"NR","R")
tide <- tide[rownames(subt),]
tide$Cluster <- subt$res_cluster
library(ggstatsplot)
p2 <- ggbarstats(tide, y = 'Cluster', x = 'Responder', )+
  scale_fill_manual(values = c('R'= '#FF1493','NR'='#40E0D0'))
p2
ggsave(plot = p2, filename = 'GSE146173_barstat.pdf', height=7, width = 6)
