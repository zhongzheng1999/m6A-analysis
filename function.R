library(ConsensusClusterPlus)
library(pheatmap)
library(survival)
library(survminer)
library(patchwork)


CCP <- function(matrix, distance, clusterAlg){
  a = ConsensusClusterPlus(
    matrix,
    maxK=7,
    reps=1000,
    pItem=0.8,
    pFeature=1,
    distance=distance,
    clusterAlg=clusterAlg,
    title=paste(distance, clusterAlg, sep='_'),
    plot='pdf'
  )
  
  return(a)
}

library(DESeq2)
createList <- function(group=NULL) {
  
  tumorsam <- names(group)
  sampleList = list()
  treatsamList =list()
  treatnameList <- c()
  ctrlnameList <- c()
  
  if(length(table(group)) >= 2){
    #A-1: 类1 vs 其他
    sampleList[[1]] = tumorsam
    treatsamList[[1]] = intersect(tumorsam, names(group[group=='C1'])) # 亚型名称需要根据情况修改
    treatnameList[1] <- 'C1' # 该亚型的命名
    ctrlnameList[1] <- "Others" # 其他亚型的命名
    
    #A-2: 类2 vs 其他
    sampleList[[2]] = tumorsam
    treatsamList[[2]] = intersect(tumorsam, names(group[group=='C2']))
    treatnameList[2] <- 'C2'
    ctrlnameList[2] <- "Others"}
  if(length(table(group)) >= 3){
      #A-3: 类3 vs 其他
      sampleList[[3]] = tumorsam
      treatsamList[[3]] = intersect(tumorsam, names(group[group=='C3']))
      treatnameList[3] <- 'C3'
      ctrlnameList[3] <- "Others"
    }
  if(length(table(group)) == 4){
        sampleList[[4]] = tumorsam
        treatsamList[[4]] = intersect(tumorsam, names(group[group=='C4']))
        treatnameList[4] <- 'C4'
        ctrlnameList[4] <- "Others"}

    
  
  #如果有更多类，按以上规律继续写
  
  return(list(sampleList, treatsamList, treatnameList, ctrlnameList))
}

# 配对DESeq2
twoclassDESeq2 <- function(res.path=NULL, countsTable=NULL, prefix=NULL, complist=NULL, overwt=FALSE) {
  
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  allsamples <- colnames(countsTable)
  
  options(warn=1)
  for (k in 1:length(sampleList)) { # 循环读取每一次比较的内容
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]] 
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    
    compname <- paste(treatname, "_vs_", ctrlname, sep="") # 生成最终文件名
    tmp = rep("others", times=length(allsamples))
    names(tmp) <- allsamples
    tmp[samples]="control"
    tmp[treatsam]="treatment"
    outfile <- file.path( res.path, paste(prefix, "_deseq2_test_result.", compname, ".txt", sep="") )
    if (file.exists(outfile) & (overwt==FALSE)) { # 因为差异表达分析较慢，因此如果文件存在，在不覆盖的情况下（overwt=F）不再次计算差异表达
      cat(k, ":", compname, "exists and skipped;\n")
      next
    }
    
    saminfo <- data.frame("Type"=tmp[samples],"SampleID"=samples,stringsAsFactors = F)
    cts <- countsTable[,samples]
    coldata <- saminfo[samples,]
    
    # 差异表达过程，具体参数细节及输出结果解释，请参阅相关document
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = as.formula("~ Type")) # 设计矩阵仅包含亚型信息，若有批次效应请修改
    
    dds$Type <- relevel(dds$Type,ref = "control")
    
    dds <- DESeq(dds)
    res <- results(dds, contrast=c("Type","treatment","control"))
    
    resData <- as.data.frame(res[order(res$padj),])
    resData$id <- rownames(resData)
    resData <- resData[,c("id","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
    colnames(resData) <- c("id","baseMean","log2FC","lfcSE","stat","PValue","FDR")
    #输出到文件
    write.table(resData, file=outfile, row.names=F, col.names=T, sep="\t", quote=F)
    cat(k, ",")
  }
  options(warn=0)
}


subtype_specific_gsea <- function(msigdb=NULL,n.top=10,mode=c("up","down"),degs.list=NULL,subtype.label=NULL,nPerm.gsea=1000,minGSSize.gsea=10,maxGSSize.gsea=500,pvalueCutoff.gsea=1){
  
  MSigDB <- read.gmt(msigdb)
  GSEA.list <- top.gs <- list() #初始化结果列表
  
  if(!is.element(mode, c("up", "dn"))) { stop("mode must be up or dn!\n") }
  
  for (i in 1:n.sub) {
    degs <- degs.list[[n.sub.label[i]]]
    geneList <- degs$log2FC; names(geneList) <- rownames(degs)
    geneList <- sort(geneList,decreasing = T) # ranked gene set
    
    # 由于GSEA不可重复，所以保存GSEA对象入列表，方便下次调用
    cat(paste0("GSEA for ",subtype.label[i]," starts!\n"))
    GSEA.list[[subtype.label[i]]] <- GSEA(geneList = geneList,
                                          TERM2GENE=MSigDB,
                                          nPerm = nPerm.gsea,
                                          minGSSize = minGSSize.gsea,
                                          maxGSSize = maxGSSize.gsea,
                                          seed = T,
                                          verbose = F,
                                          pvalueCutoff = pvalueCutoff.gsea) # 输出全部的GESA结果
    
    GSEA.dat <- as.data.frame(GSEA.list[[subtype.label[i]]])
    
    if(mode == "up") {
      GSEA.dat <- GSEA.dat[order(GSEA.dat$NES,decreasing = T),] # 根据NES降序排列，也就是找特异性上调通路
    } else {
      GSEA.dat <- GSEA.dat[order(GSEA.dat$NES,decreasing = F),] # 根据NES升序排列，也就是找特异性下调通路
    }
    
    # 输出每一次GSEA结果
    write.table(GSEA.dat,paste0(subtype.label[[i]],"_degs_",mode,"_gsea.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
    
    # 亚型特异性top基因集保存入列表
    top.gs[[subtype.label[i]]] <- rownames(GSEA.dat)[1:n.top] 
  }
  
  # 构建GSVA分析需要的gene list
  gs <- list()
  for (i in as.character(unlist(top.gs))) {
    gs[[i]] <- MSigDB[which(MSigDB$ont %in% i),"gene"]
  }
  
  return(list(mode=mode,top.gs=top.gs,gs=gs))
}


Coxoutput <- function(subt = NULL, mat = NULL){
  realdata <- data.frame(row.names = rownames(subt),
                         Days = subt$OS_Time,
                         State = subt$OS_Status,
                         mat)
  Coxoutput=data.frame()
  for(i in colnames(realdata[,3:ncol(realdata)])){
    cox <- coxph(Surv(Days, State) ~ realdata[,i], data = realdata)
    pred.dat <- predict(cox, realdata)
    dat <- data.frame(realdata[1:2], exp = pred.dat)
    coxSummary = summary(cox)
    Coxoutput=rbind(Coxoutput,cbind(gene=i,HR=coxSummary$coefficients[,"exp(coef)"],
                                    z=coxSummary$coefficients[,"z"],
                                    pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                                    lower=coxSummary$conf.int[,3],
                                    upper=coxSummary$conf.int[,4]))
  }
  for(i in c(2:6)){
    Coxoutput[,i] <- as.numeric(as.vector(Coxoutput[,i]))
  }
  
  write.csv(Coxoutput, paste(getwd(), 'Coxoutput.csv'  ,sep = '/'))
  
  Coxoutput <- arrange(Coxoutput,pvalue)  %>% #按照p值排序
    filter(pvalue < 0.05)
  return(Coxoutput$gene)
  
}



if(!require(tidyverse)) install.packages("tidyverse")  
if(!require(readxl)) install.packages("readxl")
if(!require(VIM)) install.packages("VIM")
if(!require(randomForest)) install.packages("randomForest")  
if(!require(magrittr )) install.packages("magrittr ")       
if(!require(caret)) install.packages("caret")  
if(!require(e1071)) install.packages("e1071")   
if(!require(pROC)) install.packages("pROC")   
if(!require(PerformanceAnalytics)) install.packages("PerformanceAnalytics")   
if(!require(DT)) install.packages("DT")    
if(!require(partykit)) install.packages("partykit")   
if(!require(class)) install.packages("class")  
if(!require(neuralnet)) install.packages("neuralnet")