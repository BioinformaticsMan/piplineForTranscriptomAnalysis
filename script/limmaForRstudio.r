#!/usrg/bin/bash

## define the parameter
## args[1] exprSet
## args[2] group
## args[3] logFC value, e.g. 1 for two fold change; minimum value is 0.263(1.2 fold change) 
## args[4] prefix for sample name or project name
#args <- commandArgs(T)

##define workstation
getwd()
setwd("E://AD//data//GSE118553//GSE118553_series_matrix.txt//anova//sorted//DEGs/FC")


## data loading 
exprSet <- read.table("GSE118553_series_matrix_FC_merge.txt",sep="\t",header = T,row.names = 1)

# head(exprSet)
# GSM3332689 GSM3332721 GSM3332729 GSM3332733 GSM3332735 GSM3332737 GSM3332769 GSM3332773 GSM3332789 GSM3332790 GSM3332794
# LOC651959   5.833565   5.987739   5.623442   5.982772   6.083557   5.694394   5.790087   5.612731   5.779734   5.952935   5.656121
# BAPX1       2.166212   2.637456   2.679379   2.934317   2.778142   2.676598   2.765600   2.076990   2.617697   2.379449   2.281616
# LOC643637   4.582602   4.356155   4.061952   3.761383   3.700054   3.926259   4.752269   4.457427   4.424250   4.995072   3.967456
# LOC92497    3.858629   3.816531   4.218465   4.007267   4.578450   3.954026   4.260308   4.049230   3.712960   3.706672   3.974329
# LOC651957   2.315253   2.527378   2.388252   2.686880   3.087317   3.431687   2.567497   2.495661   2.601784   2.738926   2.204749
# LOC643634   2.524963   2.688440   2.663570   3.346058   2.544646   2.864712   2.804768   2.932951   2.592211   2.365190   3.374501

#########################################################################################
group <- read.table("GSE118553_series_matrix_FC_merge_group.txt",header=T,row.names = 1)
# head(group)
#             state
# GSM3332689  case
# GSM3332721  case
# GSM3332729  case
# GSM3332733  case
# GSM3332735  case
# GSM3332737  case

group_list = as.character(group[,1])
## set the work domain
#getwd()
#setwd("E://AD//data//GSE118553//GSE118553_series_matrix.txt//anova//sorted//FC")


## correct the value
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))

LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex) }
print("======correct the value========")

## visualize the result before and after data correction
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(exprSet, col = cols,main="expression value",las=2)
suppressMessages(library(limma))


## turn group infomation into numbers
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)

## difference analysis
##step1
fit <- lmFit(exprSet,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
write.csv(nrDEG,paste("FC","all_limma.csv",sep = "_"),quote = F)

diff_gene_limma <-subset(nrDEG, adj.P.Val < 0.05 & abs(logFC) > as.numeric(0.263))## fold change = 1.6

#dim(diff_gene_limma)

write.csv(diff_gene_limma,paste("FC","DEGs_limma.csv",sep = "_"),quote=F)


## volcano file
library(ggpubr)
library(ggthemes)
deg.data = nrDEG
head(deg.data)
deg.data$logP <- -log10(deg.data$adj.P.Val)
#ggscatter(deg.data, x="logFC", y="logP")+theme_base()

deg.data$Group = "not-significant"
head(deg.data)

deg.data$Group[which( (deg.data$adj.P.Val<0.05) &  (deg.data$logFC > as.numeric(0.263)) )] = "up-regulated"
deg.data$Group[which( (deg.data$adj.P.Val<0.05) &  (deg.data$logFC< -as.numeric(0.263)) )] = "down-regulated"
table(deg.data$Group)
deg.data$GeneID = rownames(deg.data)

deg.data$Label = ""
deg.data <- deg.data[order(deg.data$adj.P.Val),]
up.genes <- head(deg.data$GeneID[which(deg.data$Group == "up-regulated")], 10)
down.genes <- head(deg.data$GeneID[which(deg.data$Group == "down-regulated")], 10)
deg.top10.genes <- c(as.character(up.genes), as.character(down.genes))
deg.top10.genes
deg.data$Label[match(deg.top10.genes, deg.data$GeneID)] <- deg.top10.genes


tiff(filename = "FC_volcano_label_top10_genes.tiff",res=300,width = 2000,height = 2000)
ggscatter(deg.data, x="logFC", y="logP", color = "Group", 
          palette = c("#2f5688", "#BBBBBB", "#CC0000"), 
          size = 1,label = deg.data$Label,font.label = 6,repel = T,
          xlab="log2FoldChage",ylab="-log10(adj.P.Val)",alpha=1)+
  theme_base()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  geom_vline(xintercept = c(-as.numeric(0.263),as.numeric(0.263)),linetype="dashed")

dev.off()


## heatmap
## nrDEG: limma results
## exprSet: expression matrix
## %>%  means use the left value to  right expression
library(dplyr)   ## use tibble  (a data structure similar with data.frame)
#head(nrDEG)
up_50 <- nrDEG %>% as_tibble() %>%
  mutate(genename=rownames(nrDEG)) %>%
  dplyr::arrange((desc(logFC))) %>%
  .$genename %>% .[1:50]

down_50 <- nrDEG %>% as_tibble() %>%
  mutate(genename=rownames(nrDEG)) %>%
  dplyr::arrange(logFC) %>%
  .$genename %>% .[1:50]
index<-c(up_50,down_50) 

library(pheatmap)
#pheatmap(exprSet[index,],show_colnames =F,show_rownames = F)
index_matrix <- -t(scale(t(exprSet[index,])))
index_matrix[index_matrix>1]=1 
##

## add annotion
anno=data.frame(group=group_list)
rownames(anno)=colnames(index_matrix)
#head(anno)
tiff(filename = "FC_heatmap_top50_genes.tiff",res=300,width = 4000,height = 4500)
pheatmap(index_matrix,
         show_colnames = F,
         show_rownames = T,
         cluster_cols = T,
         color = colorRampPalette(colors = c("#2f5688","white","#CC0000"))(100),
         annotation_col = anno)
dev.off()
