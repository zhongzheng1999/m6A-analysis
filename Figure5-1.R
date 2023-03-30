#TCGA
source("H:/Project/Main/20230105m6A/two/workflow/function.R", echo=TRUE)
#step1
load("H:/Project/Main/20230105m6A/two/subtype/data/TCGA_vst.Rdata")
vsd <- t(expr)
# standard <- preProcess(vsd, method = c("center","scale"))
# vsd_std <- predict(standard, vsd)
# dat <- as.data.frame(vsd_std)
dat = sweep(vsd,2, apply(vsd,2,median,na.rm=T))


#加载用于CCP的基因
marker_list <- readRDS("H:/Project/Main/20230105m6A/two/subtype/marker_list_merge.Rds")
marker_list[['all']] <- unique(unlist(marker_list[1:4]))

# saveRDS(marker_list,'marker_list_merge.Rds')

for (filefold in names(marker_list)[c(1:4,18)]) {
  setwd('H:/Project/Main/20230105m6A/two/workflow/TCGA')
  damp <- marker_list[[filefold]]
  gene <- intersect(colnames(dat), damp)
  matrix <- t(dat[,gene])
  clin_dat <- subt
  dir.create(paste(getwd(), filefold, sep = '/'))
  setwd(paste(getwd(), filefold, sep = '/'))
  geneCox = Coxoutput(mat = dat[,gene], subt = subt)
  matrix <- t(dat[,gsub('\\.','-',geneCox)])

  if(all(rownames(clin_dat) == colnames(matrix))){

    for (distance in c('pearson', 'spearman', 'euclidean')){
      for (clusterAlg in c('pam','km', 'hc')){
        rcc <- CCP(matrix, distance, clusterAlg)
        for (i in c(2, 3, 4)){
          
          anno2 <- data.frame(Subtype = rcc[[i]]$consensusClass)
          anno2 <- cbind(clin_dat, anno2)
          anno2$OS.time <- as.numeric(anno2$OS_Time)
          anno2$OS <- as.numeric(anno2$OS_Status)
          surv_object = Surv(anno2$OS.time, anno2$OS)
          fit1 <- survfit(surv_object ~ Subtype, data = anno2)
          summary(fit1)
          write.csv(anno2, paste(getwd(),paste(distance, clusterAlg,sep = '_'), paste(i, '.csv' ,sep = '') ,sep = '/'))
          
          # pdf('/Users/ranpeng/Desktop/Xcell_score/results(2021-1-25)/1.pdf', height=6, width = 8)
          
          p <- ggsurvplot(fit1,palette = 'lancet',
                          risk.table =TRUE,pval =TRUE,
                          conf.int =F,xlab ="Time in Days",
                          ggtheme =theme_light(),
                          ncensor.plot = F)
          p1 <- p$plot/p$table + plot_layout(heights = c(3, 1))
          ggsave(plot = p1, filename = paste(getwd(), paste(distance, clusterAlg,sep = '_'), paste(i, '.pdf' ,sep = '') ,sep = '/'), height=7, width = 8)
          
          ##TIDE
          tide <- read.csv("H:/Project/Main/20230105m6A/two/imm/data/TCGA_TIDE.csv",
                           row.names = 1, header = T)
          tide$Responder <- ifelse(str_detect(tide$Responder,"False"),"NR","R")
          anno2 <- cbind(tide, anno2)
          library(ggstatsplot)
          p2 <- ggbarstats(anno2, y = 'Subtype', x = 'Responder', )+
            scale_fill_manual(values = c('R'= '#FF1493','NR'='#40E0D0'))
          
          ggsave(plot = p2, filename = paste(getwd(), paste(distance, clusterAlg,sep = '_'), paste(i, 'barstat.pdf' ,sep = '') ,sep = '/'), height=7, width = 6)
          
        }
      }
    }
  }
}
