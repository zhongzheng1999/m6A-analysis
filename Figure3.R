
library(bigSCale)
library(SingleCellExperiment)
library(Seurat)
library(igraph)
library(tidyverse)

sce <- readRDS("H:/Project/Main/20230105m6A/step1/sce0131.Rds")
Idents(sce)
sce$m6A_bi <- factor(sce$m6A_bi,levels = c(0,1),labels = c('Inactive','Active'))

sub_leu <- subset(sce, subset = CellType_sum == 'leukemic cell')
sub_imm <- subset(sce, subset = CellType_sum == 'immune cell')


sub_leu.ctl <- subset(sub_leu, subset = m6A_bi == 'Inactive')
sub_leu.t2d <- subset(sub_leu, subset = m6A_bi == 'Active')
sub_imm.ctl <- subset(sub_imm, subset = m6A_bi == 'Inactive')
sub_imm.t2d <- subset(sub_imm, subset = m6A_bi == 'Active')


expr.ctl <- sub_leu.ctl@assays$RNA@counts
expr.t2d <- sub_leu.t2d@assays$RNA@counts

model=compute.network.model(expr.data = cbind(expr.ctl,expr.t2d))

gene.names <- rownames(expr.ctl)
results.ctl = compute.network(expr.data = expr.ctl,gene.names = gene.names,model = model, quantile.p = 0.7)
results.t2d = compute.network(expr.data = expr.t2d,gene.names = gene.names,model = model, quantile.p = 0.7)

save(results.ctl,results.t2d,file = 'result_list_leu_0.7.Rdata')


expr.ctl <- sub_imm.ctl@assays$RNA@counts
expr.t2d <- sub_imm.t2d@assays$RNA@counts

model=compute.network.model(expr.data = cbind(expr.ctl,expr.t2d))

gene.names <- rownames(expr.ctl)
results.ctl = compute.network(expr.data = expr.ctl,gene.names = gene.names,model = model, quantile.p = 0.9)
results.t2d = compute.network(expr.data = expr.t2d,gene.names = gene.names,model = model, quantile.p = 0.9)

save(results.ctl,results.t2d,file = 'result_list_imm_0.9.Rdata')




output=homogenize.networks(list(results.ctl,results.t2d))
results.ctl=output[[1]]
results.t2d=output[[2]]
saveRDS(output,'output_leu_0.9.Rds')

comparison=compare.centrality(list(results.ctl$centrality,results.t2d$centrality),c('Control','IPF'))
DT::datatable(comparison$PAGErank)
