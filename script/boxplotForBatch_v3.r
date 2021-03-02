rm(list=ls())
# args <- commandArgs(T)


setwd("/home/tonyleao/wkd/AD_manuscript_revise/boxplot/Baak-stage/TC/")

## 01 load data
datExpr0   <- read.table("GSE131617-GPL5175_expr_TC_Baak_stage.txt",header = T,sep = "\t",row.names = 1)
# clinical   <- read.csv("TCGA-ESCA_clinical_clear.csv",header = T,row.names = 1)
targetGene <- read.csv("TC-hubGene_all.csv",header = T)
DEGs       <- read.csv("TC_genes_byANOVA_FDR.csv",header = T,row.names = 1)

## 02 extract data based on gene_list
# datExpr0 <- t(log2(1+datExpr0))
gene_list  <- targetGene[,1]
datExpr    <- datExpr0[rownames(datExpr0)%in%gene_list,] 
# datExpr    <- log2(datExpr+1)
clinical <- as.data.frame(t(datExpr0[nrow(datExpr0),]))

identical(rownames(clinical),colnames(datExpr))

# datExpr[,c(1:23)] = as.numeric(unlist(datExpr[,c(1:23)]))

datExpr    <- rbind(datExpr,as.character(t(clinical$group)))
rownames(datExpr)[nrow(datExpr)] = "group"
datExpr    <- as.data.frame(t(datExpr))

## 03 remove samples with missing value
# datExpr    <- datExpr[datExpr$tumor_stage %in% c("I","II","III","IV"),]
# datExpr$group = rep("NA",nrow(datExpr))
# datExpr$group[datExpr$tumor_stage=='I']=1
# datExpr$group[datExpr$tumor_stage=='II']=2
# datExpr$group[datExpr$tumor_stage=='III']=3
# datExpr$group[datExpr$tumor_stage=='IV']=4
dim(datExpr)

## convert the expression data to numeric
geneNum <- ncol(datExpr)-1
datExpr[,c(1:geneNum)] = as.numeric(unlist(datExpr[,c(1:geneNum)]))
# boxplot(MFAP2~tumor_stage,data=datExpr,col=c("red","yellow","blue","green"))

## 04 order the gene with P.Val
genePval   <- DEGs[rownames(DEGs)%in%gene_list,]
genePval   <- genePval[order(genePval$P_value,decreasing = F),]

prefix = "TC"

myData <- as.data.frame(apply(datExpr[,1:(ncol(datExpr)-1)],2,as.numeric))
rownames(myData) = rownames(datExpr)
group <- data.frame(row.names =rownames(datExpr),group= datExpr[,ncol(datExpr)])
identical(rownames(myData),rownames(group))
# myData = as.data.frame(matrix(unlist(datExpr)))


# ## boxplot
# tiff(file=paste(prefix,"Boxplot of Hub Gene Expression.tiff", sep="_"),res=300,height=3200,width=3300)
#   par(mfrow=c(4,5))
#   for (i in rownames(genePval)){
#     y <- myData[,i]
#     # x <- myData$group
#     x <- group$group
#     #boxplot(y~x,main=paste(as.character(i),"Expression Status",sep=" "),
#     boxplot(y~x,main=paste("Expression of",as.character(i),sep=" "),
#             data=myData,
#             ylab ="Gene Expression",
#             xlab = paste("Baak stage",as.character(i),"\n",
#                          "P.Val = ",round(genePval[as.character(i),"P_value"],8),
#                          "FDR = ",round(genePval[as.character(i),"adj.P.value"],8),
#                          sep=" "),
#             las=1,
#             col=c("#FE4365","#FC9D9A","#F9CDAD","#83AF9B"))
#             # col=c("orange", "grey"))
#   }
# dev.off()


library(beeswarm)
tiff(file=paste(prefix,"Boxplot of Hub Gene Expression_beeswarm.tiff", sep="_"),res=300,height=2400,width=3300)
par(mfrow=c(4,5))
for (i in rownames(genePval)){
  y <- myData[,i]
  # x <- myData$group
  x <- group$group
  #boxplot(y~x,main=paste(as.character(i),"Expression Status",sep=" "),
  beeswarm(y~x,main=paste("Expression of",as.character(i),sep=" "),
          data=myData,
          ylab ="Gene Expression",
          xlab = paste("Baak stage of",as.character(i),"\n",
                       "P.Val = ",round(genePval[as.character(i),"P_value"],8),
                       "FDR = ",round(genePval[as.character(i),"adj.P.value"],8),
                       sep=" "),
          las=1,
          col=c("#FE4365","#FC9D9A","#F9CDAD","#83AF9B"))
  # col=c("orange", "grey"))
}
dev.off()



# 
# 
# tiff(file=paste(prefix,"Boxplot of Hub Gene Expression_beeswarm.tiff", sep="_"),res=300,height=1600,width=3300)
# par(mfrow=c(2,5))
# for (i in rownames(genePval)){
#   y <- myData[,i]
#   x <- myData$group
#   #boxplot(y~x,main=paste(as.character(i),"Expression Status",sep=" "),
#   beeswarm(y~x,main=paste("Expression of",as.character(i),sep=" "),
#           data=myData,
#           ylab ="Gene Expression",
#           xlab = paste("Tumor status of",as.character(i),"\n","adj.P.Val = ",round(genePval[as.character(i),"adj.P.Val"],8),sep=" "),
#           las=1,
#           col=c("#83AF9B","#F9CDAD"))
#   # col=c("orange", "grey"))
# }
# dev.off()
# 
# 
# 
# library(vioplot)
# tiff(file=paste(prefix,"Boxplot of Hub Gene Expression_vioplot.tiff", sep="_"),res=300,height=1600,width=3300)
# par(mfrow=c(2,5))
# for (i in rownames(genePval)){
#   y <- myData[,i]
#   x <- myData$group
#   #boxplot(y~x,main=paste(as.character(i),"Expression Status",sep=" "),
#   vioplot(y~x,main=paste("Expression of",as.character(i),sep=" "),
#            data=myData,
#            ylab ="Gene Expression",
#            xlab = paste("Tumor status of",as.character(i),"\n","adj.P.Val = ",round(genePval[as.character(i),"adj.P.Val"],8),sep=" "),
#            las=1,
#            col=c("#83AF9B","#F9CDAD"))
#   # col=c("orange", "grey"))
# }
# dev.off()
# 
# 
# ##ref https://zhuanlan.zhihu.com/p/142161067
# library(beanplot)
# tiff(file=paste(prefix,"Boxplot of Hub Gene Expression_beanplot.tiff", sep="_"),res=300,height=1600,width=3300)
# par(mfrow=c(2,5))
# for (i in rownames(genePval)){
#   y <- myData[,i]
#   x <- myData$group
#   #boxplot(y~x,main=paste(as.character(i),"Expression Status",sep=" "),
#   beanplot(y~x,main=paste("Expression of",as.character(i),sep=" "),
#           data=myData,
#           ylab ="Gene Expression",
#           xlab = paste("Tumor status of",as.character(i),"\n","adj.P.Val = ",round(genePval[as.character(i),"adj.P.Val"],8),sep=" "),
#           las=1,
#           col=c("#83AF9B","#F9CDAD"))
#   # col=c("orange", "grey"))
# }
# dev.off()
# 
# 
# # 
# # ##beanplot modify
# # library(beanplot)
# # tiff(file=paste(prefix,"Boxplot of Hub Gene Expression_beanplot_modify.tiff", sep="_"),res=300,height=1600,width=3300)
# # par(mfrow=c(2,5))
# # for (i in rownames(genePval)){
# #   y <- myData[,i]
# #   x <- myData$group
# #   #boxplot(y~x,main=paste(as.character(i),"Expression Status",sep =" "),
# #   beanplot(y~ reorder(supp,y,mean)*x,data=myData,
# #            side = "b", col = list("#83AF9B","#F9CDAD", 
# #            border = c("#83AF9B","#F9CDAD"),
# #            main=paste("Expression of",as.character(i),sep=" "),
# #            ylab ="Gene Expression",
# #            xlab = paste("Tumor stage of",as.character(i),"\n","adj.P.Val = ",
# #                         round(genePval[as.character(i),"adj.P.Val"],8),sep=" ")))
# # 
# # }
# # dev.off()
# 
# 
# 
# 
# 
# beanplot(len ~ reorder(supp, len, mean) * dose, ToothGrowth, 
#          side = "b", col = list("yellow", "orange"), border = c("yellow2", 
#                                                                 "darkorange"), main = "Guinea Pigs' Tooth Growth", 
#          xlab = "Vitamin C dose mg", ylab = "tooth length", ylim = c(-1, 
#                                                                      40), yaxs = "i")
# 

