#step2
source("H:/Project/Main/20230105m6A/two/workflow/function.R", echo=TRUE)
file_dir = list.files('H:/Project/Main/20230105m6A/two/workflow/TCGA_stand_vst_surv/all',
                      full.names = T)
file_dir = file_dir[-1]
file_dir
marker_list <- readRDS("H:/Project/Main/20230105m6A/two/subtype_DEG/marker_list_merge.Rds")

for (i in file_dir) {
  setwd(i)
  count_dat <- readRDS("H:/Project/Main/20230105m6A/two/subtype_DEG/tcga_count.Rds")
  count_dat <- count_dat[which(rowSums(count_dat > 10) > 100) ,]
  for (a in list.files('.','csv$')) {
    subt <- read.csv(a, check.names = F, stringsAsFactors = F, header = T, row.names = 1)
    head(subt)
    expr <- count_dat
    
    sample_index <- intersect(colnames(expr),rownames(subt))
    
    expr <- expr[sample_index]
    subt <- subt[sample_index,]
    subt$TCGA_Subtype <- paste0('C',subt$Subtype)
    
    n.sub.label <- unique(subt$TCGA_Subtype) #亚型名称
    n.sub.label
    
    if(min(table(subt$TCGA_Subtype)) > 10){
      n.sub <- length(table(subt$TCGA_Subtype)) #亚型个数
      n.sub
      
      group <- subt$TCGA_Subtype
      names(group) <- rownames(subt) 
      complist <- createList(group=group)
      
      ### 执行配对DESeq2差异表达 ###
      
      # 差异表达分析过程比较慢请耐心等待，这里为了加速分析过程，只选取方差top25%的基因
      # var <- apply(expr,1,sd)
      # index <- var > quantile(var)[4]
      
      twoclassDESeq2(res.path = ".", #所有配对差异表达结果都会输出在res.path路径下
                     countsTable = expr[,intersect(colnames(expr),rownames(subt))],
                     prefix = paste0("AML_TCGA_",a), #文件名以SKCM开头
                     complist = complist,
                     overwt = F)
      
      DEfiles <- list.files(path = '.',pattern = paste0("AML_TCGA_",a))
      
      degs.list <- list()
      for (b in DEfiles) {
        degs <- read.table(b,sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
        head(degs)
        degs.list[[b]] <- as.data.frame(na.omit(degs))
      }
      
      
      df <- lapply(degs.list, function(x){
        x <- x[marker_list[['all']],]
        x$gene <- rownames(x)
        return(x)
      })
      
      df <- do.call(rbind, df)
      df <- subset(df, subset = FDR < 0.01 )
      
      hub_gene <- unique(df$gene)
      
      load("H:/Project/Main/20230105m6A/two/subtype/data/TCGA_vst.Rdata")
      # subt$TCGA_Subtype <- paste0('C',subt$Cluster)
      subt <- read.csv(a, check.names = F, stringsAsFactors = F, header = T, row.names = 1)
      head(subt)
      sample_index <- intersect(colnames(expr),rownames(subt))
      
      expr <- expr[sample_index]
      subt <- subt[sample_index,]
     
      data <- as.data.frame(t(expr[hub_gene,]))
      data$NSP <- subt$Subtype
      
      data <- data %>%  #先全部变数值型,再将三个分类变量因子化
        mutate_if(.predicate = is.character,.funs =
                    as.numeric) %>%
        mutate_at(.vars =
                    vars(NSP),.funs = as.factor) %>%
        glimpse()
      
      
      train_control <- trainControl(method = 'cv', number = 5)#分层抽样 
      set.seed(42) # 设置随机种子
      colnames(data) <- gsub('-','.',colnames(data))# 根据数据的因变量进行7:3的分层抽样
      index <-createDataPartition(data$NSP, p = 0.7, list = FALSE) 
      index <- as.vector(index)
      traindata <- data[index,] # 提取数据中的index所对应行索引的数据作为训练集
      testdata <- data[-index,] # 其余的作为测试集
      
      # 对数据进行标准化
      standard <- preProcess(traindata, method = c("center","scale"))
      traindata_std <- predict(standard, traindata)
      testdata_std <- predict(standard, testdata)
      
      #========================== 1、回归树 ======================
      #建立模型
      rpart_model <- train(NSP~ ., data = traindata_std, trControl = train_control, method = 'rpart')
      rpart_pred <- predict(rpart_model, testdata[-ncol(testdata_std)])
      rpart_result <- confusionMatrix(rpart_pred, testdata_std$NSP) 
      #========================= 2、随机森林==========================
      ## 建立随机森林模型进行预测，并可视化变量重要性
      rf_model <- randomForest(NSP ~., data = traindata_std, importance =T)
      rf_pred <- predict(rf_model, testdata_std[-ncol(testdata_std)], type = 'class')
      rf_result <- confusionMatrix(rf_pred, testdata_std$NSP)
      # =======================3、朴素贝叶斯 ========================
      nb_model <- train(NSP ~., data = traindata_std,trControl = train_control,
                        method = 'nb')
      nb_pred <- predict(nb_model, testdata_std[-ncol(testdata_std)])
      nb_result <- confusionMatrix(nb_pred, testdata_std$NSP) 
      #=======================4、支持向量机 =====================
      svm_model <- svm(NSP~.,data = traindata_std, kernel = 'radial')
      svm_pred <- predict(svm_model, testdata_std[-ncol(testdata_std)], type = 'response')
      svm_result <- confusionMatrix(svm_pred, testdata_std$NSP) 
      #========================= 5、knn ==============================
      results = c()
      for(c in 3:10) {
        set.seed(1234)
        knn_pred <- knn(traindata_std[-ncol(testdata_std)],
                        testdata_std[-ncol(testdata_std)], traindata_std$NSP, c)
        Table <- table(knn_pred, testdata_std$NSP)
        accuracy <- sum(diag(Table))/sum(Table) #diag()
        results <- c(results, accuracy)
      }
      k = which.max(results)+2
      knn_pred <- knn(train = traindata_std[-ncol(testdata_std)], test = testdata_std[-ncol(testdata_std)], cl = traindata_std$NSP, k = k)
      knn_result <- confusionMatrix(knn_pred , testdata_std$NSP)
      #==================== 6、bp神经网络 =====================
      set.seed(1234)
      traindata_std1 <- traindata_std[,1:ncol(testdata_std)-1]
      traindata_std1$one = traindata_std$NSP == 1
      traindata_std1$two = traindata_std$NSP == 2
      traindata_std1$three = traindata_std$NSP == 3
  
      bp_modal <- neuralnet(one+two+three ~.,data=traindata_std1,threshold=0.01,
                            stepmax=100000,err.fct="ce",linear.output=F,hidden=5)
      testdata_std1 <- testdata_std[,1:ncol(testdata_std)-1]
      testdata_std1$one = testdata_std$NSP == 1
      testdata_std1$two = testdata_std$NSP == 2
      testdata_std1$three = testdata_std$NSP == 3
      bp_result <- compute(bp_modal,testdata_std1)
      bp_pred =c("one","two","three")[apply(bp_result$net.result,1, which.max)] %>% as.data.frame()
      bp_pred_num <- bp_pred %>% 
        mutate(pred=case_when(
          bp_pred[,1]=='one' ~ 1,
          bp_pred[,1]=='two' ~ 2,
          bp_pred[,1]=='three' ~ 3,
        ))
      bp_pred1 <- bp_pred_num[[2]] %>% as.factor()
      bp_cm <- confusionMatrix(bp_pred1,testdata_std$NSP)
      save(traindata_std,testdata_std,file = paste0(i,'/',a,'_TCGA_std.Rdata'))
      save(rpart_model, rpart_pred, rpart_result, rf_model, rf_pred, rf_result,
           svm_model, svm_pred, svm_result, nb_model, nb_result, nb_pred,
           knn_result, knn_pred, bp_modal, bp_pred1, bp_result,
           file = paste0(i,'/',a,'_TCGA.Rdata'))
      rm(rpart_pred, rpart_result, rf_pred, rf_result,
         svm_pred, svm_result, nb_result, nb_pred, knn_result,
         knn_pred, bp_pred1, bp_result)
      rm(expr, subt)
      gc()
      
      load("H:/Project/Main/20230105m6A/two/subtype/data/GSE106291_vst.Rdata")
      data <- as.data.frame(t(GSE106291_fpkm[hub_gene,])) #下载数据为excel文件
      colnames(data) <- hub_gene
      data[is.na(data)] = 0
      
      colnames(data) <- gsub('-','.',colnames(data))
      standard <- preProcess(data, method = c("center","scale"))
      testdata_std <- predict(standard, data)
      # rpart_pred <- predict(rpart_model, testdata)
      rf_pred <- predict(rf_model, testdata_std, type = 'class')
      nb_pred <- predict(nb_model, testdata_std)
      svm_pred <- predict(svm_model, testdata_std, type = 'response')
      knn_pred <- knn(train = traindata_std[-ncol(traindata_std)], test = testdata_std, cl = traindata_std$NSP, k = k)
      bp_result <- compute(bp_modal,testdata_std)
      bp_result$net.result
      bp_pred =c("one","two","three")[apply(bp_result$net.result,1, which.max)] %>% as.data.frame()
      bp_pred_num <- bp_pred %>% 
        mutate(pred=case_when(
          bp_pred[,1]=='one' ~ 1,
          bp_pred[,1]=='two' ~ 2,
          bp_pred[,1]=='three' ~ 3,
        ))
      bp_pred1 <- bp_pred_num[[2]] %>% as.factor()
      
      save(rf_pred, nb_pred, svm_pred, knn_pred, bp_pred1,
           file = paste0(i,'/',a,'_GSE106291_cluster.Rdata'))
      
      library(survival)
      library(survminer)
      surv_list <- list()
      tide_list <- list()
      for (d in c('rf_pred','nb_pred','svm_pred','knn_pred','bp_pred1')) {
        anno2 <- data.frame(Cluster = get(d))
        anno2 <- cbind(clinical_annot, anno2)
        anno2$OS.time <- as.numeric(anno2$OS_Time)
        anno2$OS <- as.numeric(anno2$OS_Status)
        surv_object = Surv(anno2$OS.time, anno2$OS)
        fit1 <- survfit(surv_object ~ Cluster, data = anno2)
        summary(fit1)
        # write.csv(anno2, paste(getwd(),paste(distance, clusterAlg,sep = '_'), paste(i, '.csv' ,sep = '') ,sep = '/'))
        
        # pdf('/Users/ranpeng/Desktop/Xcell_score/results(2021-1-25)/1.pdf', height=6, width = 8)
        
        p <- ggsurvplot(fit1,palette = 'lancet',
                        risk.table =TRUE,pval =TRUE,
                        conf.int =F,xlab ="Time in Days",
                        ggtheme =theme_light(),
                        ncensor.plot = F)
        p
        p1 <- p$plot/p$table + plot_layout(heights = c(3, 1))
        surv_list[[d]] <- p1
        
        ##TIDE
        tide <- read.csv("H:/Project/Main/20230105m6A/two/imm/data/GSE106291_TIDE.csv",
                         row.names = 1, header = T)
        tide$Responder <- ifelse(str_detect(tide$Responder,"False"),"NR","R")
        anno2 <- cbind(tide, anno2)
        library(ggstatsplot)
        p2 <- ggbarstats(anno2, y = 'Cluster', x = 'Responder', )+
          scale_fill_manual(values = c('R'= '#FF1493','NR'='#40E0D0'))
        tide_list[[d]] <- p2
      }
      surv_plot <- cowplot::plot_grid(plotlist = surv_list)
      ggsave(plot = surv_plot, filename =  paste0(i,'/',a,'_GSE106291_surv.pdf'),
             width = 12, height = 9)
      
      tide_plot <- cowplot::plot_grid(plotlist = tide_list)
      ggsave(plot = tide_plot, filename = paste0(i,'/',a,'GSE106291_tide.pdf'),
             width = 12, height = 9)
      
      save(surv_list, tide_list, file =  paste0(i,'/',a,'_GSE106291_sur_tide.Rdata'))
      
      ####GSE146173
      load("H:/Project/Main/20230105m6A/two/subtype/data/GSE146173_vst.Rdata")
      data <- as.data.frame(t(GSE146173_fpkm[hub_gene,])) #下载数据为excel文件
      
      rownames(data) == rownames(clinical_annot)
      colnames(data) <- hub_gene
      data[is.na(data)] = 0
      
      colnames(data) <- gsub('-','.',colnames(data))
      standard <- preProcess(data, method = c("center","scale"))
      testdata_std <- predict(standard, data)
      rpart_pred <- predict(rpart_model, testdata)
      rf_pred <- predict(rf_model, testdata_std, type = 'class')
      nb_pred <- predict(nb_model, testdata_std)
      svm_pred <- predict(svm_model, testdata_std, type = 'response')
      knn_pred <- knn(train = traindata_std[-ncol(traindata_std)], test = testdata_std, cl = traindata_std$NSP, k = k)
      bp_result <- compute(bp_modal,testdata_std)
      bp_result$net.result
      bp_pred =c("one","two","three")[apply(bp_result$net.result,1, which.max)] %>% as.data.frame()
      bp_pred_num <- bp_pred %>% 
        mutate(pred=case_when(
          bp_pred[,1]=='one' ~ 1,
          bp_pred[,1]=='two' ~ 2,
          bp_pred[,1]=='three' ~ 3,
        ))
      bp_pred1 <- bp_pred_num[[2]] %>% as.factor()
      
      save(rf_pred, nb_pred, svm_pred, knn_pred, bp_pred1,
           file = paste0(i,'/',a,'_GSE146173_cluster.Rdata'))
      
      library(survival)
      library(survminer)
      surv_list <- list()
      tide_list <- list()
      for (d in c('rf_pred','nb_pred','svm_pred','knn_pred','bp_pred1')) {
        anno2 <- data.frame(Cluster = get(d))
        anno2 <- cbind(clinical_annot, anno2)
        anno2$OS.time <- as.numeric(anno2$OS_Time)
        anno2$OS <- as.numeric(anno2$OS_Status)
        surv_object = Surv(anno2$OS.time, anno2$OS)
        fit1 <- survfit(surv_object ~ Cluster, data = anno2)
        summary(fit1)
        # write.csv(anno2, paste(getwd(),paste(distance, clusterAlg,sep = '_'), paste(i, '.csv' ,sep = '') ,sep = '/'))
        
        # pdf('/Users/ranpeng/Desktop/Xcell_score/results(2021-1-25)/1.pdf', height=6, width = 8)
        
        p <- ggsurvplot(fit1,palette = 'lancet',
                        risk.table =TRUE,pval =TRUE,
                        conf.int =F,xlab ="Time in Days",
                        ggtheme =theme_light(),
                        ncensor.plot = F)
        p
        p1 <- p$plot/p$table + plot_layout(heights = c(3, 1))
        surv_list[[d]] <- p1
        
        ##TIDE
        tide <- read.csv("H:/Project/Main/20230105m6A/two/imm/data/GSE146173_TIDE.csv",
                         row.names = 1, header = T)
        tide$Responder <- ifelse(str_detect(tide$Responder,"False"),"NR","R")
        anno2 <- cbind(tide, anno2)
        library(ggstatsplot)
        p2 <- ggbarstats(anno2, y = 'Cluster', x = 'Responder', )+
          scale_fill_manual(values = c('R'= '#FF1493','NR'='#40E0D0'))
        tide_list[[d]] <- p2
      }
      surv_plot <- cowplot::plot_grid(plotlist = surv_list)
      ggsave(plot = surv_plot, filename =  paste0(i,'/',a,'_GSE146173_surv.pdf'),
             width = 12, height = 9)
      
      tide_plot <- cowplot::plot_grid(plotlist = tide_list)
      ggsave(plot = tide_plot, filename = paste0(i,'/',a,'_GSE146173_tide.pdf'),
             width = 12, height = 9)
      
      save(surv_list, tide_list, file =  paste0(i,'/',a,'_GSE146173_sur_tide.Rdata'))
    }
  }
}


# 
# ##CV值确定
# cal_cv=function(x){  # 自定义函数 标准差/平均值
#   y=na.omit(x)
#   return(sd(y)/mean(y))
# }
# 
# tcga_fpkm <- readRDS("H:/Project/Main/20230105m6A/two/subtype/data/TCGA-VST.Rds")
# tcga_fpkm <- na.omit(tcga_fpkm[marker_list[['all']],])
# tcga_fpkm <- as.data.frame(t(tcga_fpkm))
# 
# # standard <- preProcess(tcga_fpkm , method = c("center","scale"))
# # vsd_std <- predict(standard, tcga_fpkm )
# # tcga_fpkm <- as.data.frame(vsd_std)
# 
# sample_index <- intersect(rownames(tcga_fpkm),rownames(subt))
# 
# tcga_fpkm <- tcga_fpkm[sample_index,]
# tcga_fpkm$Cluster <- subt$TCGA_Subtype
# 
# df_cv <- aggregate(.~Cluster,cal_cv,data = tcga_fpkm)
# 
# 
# ##MEGER
# df$Cluster <- str_split(rownames(df),'_',simplify = T)[,5]
# rownames(df_cv) <- df_cv$Cluster
# df_cv <- df_cv[,-1]
# df_cv <- as.data.frame(t(df_cv))
# df_cv$gene <- rownames(df_cv)
# 
# df_merge <- merge(df, df_cv, by = 'gene')
# 
# df_merge_sig <- subset(df_merge, subset = FDR < 0.05)
# df_list_sig <- split(df_merge_sig, df_merge_sig$Cluster)
# 
# save(df_merge, df_list_sig, file = 'df_cv_fc.Rdata')
# 
# plot_dat <- df_list_sig[['result.C3']]
# ggplot(plot_dat, aes(x = abs(log2FC), y = C3))+
#   geom_point()+
#   geom_vline(xintercept = c(1, 1), color = 'red', size = 0.3) +
#   geom_hline(yintercept = c(1, 1), color = 'red', size = 0.3) +
#   # scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('Up', 'Down', 'None')) +
#   labs(x = 'Log2 Fold Change', y = 'Coefficient of Variation', color = '', title = 'Cluster1')+
#   theme_bw()
# ggsave('Cluster3_CV_FC.pdf', height = 4, width = 4)

# c1_gene <- df_list_sig[['result.C1']][abs(df_list_sig$result.C1$log2FC) > 1 & df_list_sig$result.C1$C1 < 0.3,'gene']
# c2_gene <- df_list_sig[['result.C2']][abs(df_list_sig$result.C2$log2FC) > 1 & df_list_sig$result.C2$C2 < 0.3,'gene']
# c3_gene <- df_list_sig[['result.C3']][abs(df_list_sig$result.C3$log2FC) > 1 & df_list_sig$result.C3$C3 < 0.3,'gene']

# c1_gene <- df_list_sig[['result.C1']][abs(df_list_sig$result.C1$log2FC) > median(abs(df_list_sig$result.C1$log2FC)) & df_list_sig$result.C1$C1 < median(df_list_sig$result.C1$C1),'gene']
# c2_gene <- df_list_sig[['result.C2']][abs(df_list_sig$result.C2$log2FC) > median(abs(df_list_sig$result.C2$log2FC)) & df_list_sig$result.C2$C2 < median(df_list_sig$result.C2$C2),'gene']
# c3_gene <- df_list_sig[['result.C3']][abs(df_list_sig$result.C3$log2FC) > median(abs(df_list_sig$result.C3$log2FC)) & df_list_sig$result.C3$C3 < median(df_list_sig$result.C3$C3),'gene']

# c1_gene <- df_list_sig[['result.C1']][abs(df_list_sig$result.C1$log2FC) > median(abs(df_list_sig$result.C1$log2FC)),'gene']
# c2_gene <- df_list_sig[['result.C2']][abs(df_list_sig$result.C2$log2FC) > median(abs(df_list_sig$result.C2$log2FC)),'gene']
# c3_gene <- df_list_sig[['result.C3']][abs(df_list_sig$result.C3$log2FC) > median(abs(df_list_sig$result.C3$log2FC)),'gene']

# c1_gene <- df_list_sig[['result.C1']][abs(df_list_sig$result.C1$log2FC) > 1 ,'gene']
# c2_gene <- df_list_sig[['result.C2']][abs(df_list_sig$result.C2$log2FC) > 1 ,'gene']
# c3_gene <- df_list_sig[['result.C3']][abs(df_list_sig$result.C3$log2FC) > 1 ,'gene']

# 
# library(VennDiagram)
# venn <- venn.diagram(list(Cluster1 = c1_gene,
#                           Cluster2 = c2_gene,
#                           Cluster3 = c3_gene), filename = NULL)
# 
# pdf('hub_gene_veen.pdf', height = 4, width = 4)
# grid.draw(venn)
# dev.off()
# 
# openxlsx::write.xlsx(list(Cluster1 = c1_gene,
#                           Cluster2 = c2_gene,
#                           Cluster3 = c3_gene),'hub_gene.xlsx')
# 
# # all <- c(c1_gene,c2_gene,c3_gene)
# # all <- all[!all %in% all[duplicated(all)]]
# # hub_gene <- all
# 
# hub_gene <- unique(c(c1_gene,c2_gene,c3_gene))
# save(c1_gene,c2_gene,c3_gene,hub_gene,file = 'gene_list.Rdata')
