args <- commandArgs(T)

#getwd()
#setwd("E://AD//data//GSE118553//GSE118553_series_matrix.txt//anova//sorted//WGCNA//TC//TC_6//boxplot")

datExpr0 <- read.table(args[1],sep="\t",header = T,row.names = 1);
datExpr0 <- t(log2(1+datExpr0))
head(datExpr0[1:5,1:5])
datExpr <- datExpr0

datTraits <- read.table(args[2],header = T,row.names = 1)
index = c(rownames(datExpr))
datTraits <- as.data.frame(t(datTraits[,index]))

hub <- read.csv(args[3],header = T)
Diffs <- read.csv(args[4],header = T,row.names = 1) ##here just load the p value of gene
#DEGs <- read.csv("TC_DEGs_limma_logFC1.csv",head =T,row.names = 1)  ### the results were bad
prefix=args[5]

diffGene <- rownames(Diffs)
length(diffGene)
probes <- colnames(datExpr)
hubGene1 <- as.character(hub[,1])
length(hubGene1)

hubProbes = probes %in% hubGene1
hubGeneExpr <- as.data.frame(datExpr[,hubProbes])
hubGeneExpr$diagnosis = datTraits[,1]
write.csv(hubGeneExpr,paste(prefix,"hubGeneExpr.csv",sep="_"),quote=F)

hubGeneDEGsProbes <- rownames(Diffs) %in% hubGene1

hubGeneDEGs <- Diffs[hubGeneDEGsProbes,]
hubGeneDEGs <- hubGeneDEGs[order(hubGeneDEGs$adj.P.Val,decreasing=F),]
hubGeneDEGs.rownames <- rownames(hubGeneDEGs)

## boxplot
myData = hubGeneExpr

if (length(hubGeneDEGs.rownames)<=5){
tiff(file=paste(prefix,"Boxplot of Hub Gene Expression.tiff", sep=" "),res=300,height=600,width=3300)
par(mfrow=c(1,5))
for (i in hubGeneDEGs.rownames){
  y <- myData[,i]
  x <- myData[,length(myData)]
  #boxplot(y~x,main=paste(as.character(i),"Expression Status",sep=" "),
  boxplot(y~x,main=paste("Expression of",as.character(i),sep=" "),
          data=myData,
          ylab ="Gene Expression",
          xlab = paste("Diagnosis of",as.character(i),"\n","adj.P.vlaue = ",round(hubGeneDEGs[as.character(i),"adj.P.Val"],8),sep=" "),
          las=1,col=c("orange", "grey"))
}
dev.off()

}else if (length(hubGeneDEGs.rownames)<=10) {
tiff(file=paste(prefix,"Boxplot of Hub Gene Expression.tiff", sep=" "),res=300,height=1200,width=3300)
par(mfrow=c(2,5))
for (i in hubGeneDEGs.rownames){
  y <- myData[,i]
  x <- myData[,length(myData)]
  #boxplot(y~x,main=paste(as.character(i),"Expression Status",sep=" "),
  boxplot(y~x,main=paste("Expression of",as.character(i),sep=" "),
          data=myData,
          ylab ="Gene Expression",
          xlab = paste("Diagnosis of",as.character(i),"\n","adj.P.vlaue = ",round(hubGeneDEGs[as.character(i),"adj.P.Val"],8),sep=" "),
          las=1,col=c("orange", "grey"))
}
dev.off()

}else if (length(hubGeneDEGs.rownames)<=15) {
tiff(file=paste(prefix,"Boxplot of Hub Gene Expression.tiff", sep=" "),res=300,height=1800,width=3300)
par(mfrow=c(3,5))
for (i in hubGeneDEGs.rownames){
  y <- myData[,i]
  x <- myData[,length(myData)]
  #boxplot(y~x,main=paste(as.character(i),"Expression Status",sep=" "),
  boxplot(y~x,main=paste("Expression of",as.character(i),sep=" "),
          data=myData,
          ylab ="Gene Expression",
          xlab = paste("Diagnosis of",as.character(i),"\n","adj.P.vlaue = ",round(hubGeneDEGs[as.character(i),"adj.P.Val"],8),sep=" "),
          las=1,col=c("orange", "grey"))
}
dev.off()


}else if (length(hubGeneDEGs.rownames)<=20) {
tiff(file=paste(prefix,"Boxplot of Hub Gene Expression.tiff", sep=" "),res=300,height=2400,width=3300)
par(mfrow=c(4,5))
for (i in hubGeneDEGs.rownames){
  y <- myData[,i]
  x <- myData[,length(myData)]
  #boxplot(y~x,main=paste(as.character(i),"Expression Status",sep=" "),
  boxplot(y~x,main=paste("Expression of",as.character(i),sep=" "),
          data=myData,
          ylab ="Gene Expression",
          xlab = paste("Diagnosis of",as.character(i),"\n","adj.P.vlaue = ",round(hubGeneDEGs[as.character(i),"adj.P.Val"],8),sep=" "),
          las=1,col=c("orange", "grey"))
}


}else if (length(hubGeneDEGs.rownames)<=25) {
tiff(file=paste(prefix,"Boxplot of Hub Gene Expression.tiff", sep=" "),res=300,height=3000,width=3300)
par(mfrow=c(5,5))
for (i in hubGeneDEGs.rownames){
  y <- myData[,i]
  x <- myData[,length(myData)]
  #boxplot(y~x,main=paste(as.character(i),"Expression Status",sep=" "),
  boxplot(y~x,main=paste("Expression of",as.character(i),sep=" "),
          data=myData,
          ylab ="Gene Expression",
          xlab = paste("Diagnosis of",as.character(i),"\n","adj.P.vlaue = ",round(hubGeneDEGs[as.character(i),"adj.P.Val"],8),sep=" "),
          las=1,col=c("orange", "grey"))
}
dev.off()

} else {
print ("Too many hub gene!")
}


