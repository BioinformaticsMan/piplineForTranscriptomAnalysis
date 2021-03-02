args<-commandArgs(T)

library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DOSE)
library(ggplot2)
library(topGO)

#getwd()
#setwd("E://AD//data//GSE118553//GSE118553_series_matrix.txt//anova//sorted//GO_KEGG//TC//module")
##

#myData <- read.csv("MEgreenyellow_Diagnosis_state_MM.csv",header = T)
myData <- read.csv(args[1],header = T)

head(myData)
DEGs <- read.csv(args[2],header = T)
head(DEGs)


geneModule = as.character(myData[,1])
geneModule
geneDEGs = as.character(DEGs[,1])
geneDEGs
geneList = Reduce(intersect,list(geneModule,geneDEGs))
geneList

prefix = args[3] 

##########################################################################
geneID = geneList
head(geneID)

eg <- bitr(geneID, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db"); head(eg)
#new_degenes$EntrezID <- eg$ENTREZID
genelist <- eg$ENTREZID

#write.csv(genelist,"genelist.csv")
#genelist <- read.csv("genelist.csv")
genelist  <- genelist[!duplicated(genelist)]
head(genelist)
go <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 1 ,keyType = 'ENTREZID')
dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])
print("========================ALL barplot and dotplot=====================")
#ALL barplot and dotplot
tiff(file=paste(prefix,"GOALL_barplot.tiff",sep="-"),res = 300, height = 2000, width = 2000)
barplot(go,showCategory=15,drop=T, title =paste( prefix,"GO-ALL",sep="-"))
dev.off()
tiff(file=paste(prefix,"GOALL_dot.tiff",sep="-"),res = 300, height = 2500, width = 2400)
dotplot(go,showCategory=15, title = paste(prefix, "GO-ALL", sep="-"))
dev.off()

# print("=========================network BP CC MF==================================")
# #network BP CC MF
# #BiocManager::install("topGO")
# go.BP <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2, keyType = 'ENTREZID')
# #BP
# tiff(file=paste(prefix,"GOBP_network.tiff",sep="-"),res = 500, height = 3500, width = 3000)
# plotGOgraph(go.BP)
# dev.off()
# 
# #CC
# go.CC <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='CC',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
# tiff(file=paste(prefix,"GOCC_network.tiff",sep="-"),res = 500, height = 3500, width = 3000)
# plotGOgraph(go.CC)
# dev.off()
# 
# #MF
# go.MF <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID')
# tiff(file=paste(prefix,"GOMF_network.tiff",sep="-"),res = 500, height = 3500, width = 3000)
# plotGOgraph(go.MF)
# dev.off()
#print("=============================barplot&dotplot BP CC MF============================")
##barplot&dotplot BP CC MF
#tiff(file=paste(prefix,"GOBP_barplot.tiff",sep="-"),res = 300, height = 2000, width = 2000)
#barplot(go.BP,showCategory=10,drop=T, title = paste(prefix,"GO-BP",sep="-"))
#dev.off()
#tiff(file=paste(prefix,"GOBP_dot.tiff",sep="-"),res = 300, height = 2500, width = 2400)
#dotplot(go.BP,showCategory=10, title = paste(prefix,"GO-BP",sep="-"))
#dev.off()
#
##CC
#tiff(file=paste(prefix,"GOCC_barplot.tiff",sep="-"),res = 300, height = 2000, width = 2000)
#barplot(go.CC,showCategory=10,drop=T, title = paste(prefix,"GO-CC",sep="-"))
#dev.off()
#tiff(file=paste(prefix,"GOCC_dot.tiff",sep="-"),res = 300, height = 2500, width = 2400)
#dotplot(go.CC,showCategory=10, title = paste(prefix,"GO-CC",sep="-"))
#dev.off()
#
##MF
#tiff(file=paste(prefix,"GOMF_barplot.tiff",sep="-"),res = 300, height = 2000, width = 4000)
#barplot(go.MF,showCategory=10101010101010101010,drop=T, title = paste(prefix,"GO-MF",sep="-"))
#dev.off()
#tiff(file=paste(prefix,"GOMF_dot.tiff",sep="-"),res = 300, height = 2500, width = 4000)
#dotplot(go.MF,showCategory=10, title = paste(prefix,"GO-MF",sep="-"))
#dev.off()
#
#print("============================emapplot BP CC MF=====================================")

##emapplot BP CC MF
##BP
#tiff(file=paste(prefix,"GOBP_emap.tiff",sep="-"),res = 300, height = 2000, width = 3000)
#emapplot(go.BP)
#dev.off()
##CC
#tiff(file=paste(prefix,"GOCC_emap.tiff",sep="-"),res = 300, height = 2000, width = 3000)
#emapplot(go.CC)
#dev.off()
##MF
#tiff(file=paste(prefix,"GOMF_emap.tiff",sep="-"),res = 300, height = 2000, width = 3000)
#emapplot(go.MF)
#dev.off()

#print("==================================goplot BP CC MF===================================")
###goplot BP CC MF
##BP
#tiff(file=paste(prefix,"GOBP_goplot.tiff",sep="-"),res = 300, height = 2000, width = 3000)
#goplot(go.BP,showCategory = 5,layout = "sugiyama",geom = "text", title = paste(prefix,"GO-BP",sep="-"))
#dev.off()
##CC
#tiff(file=paste(prefix,"GOCC_goplot.tiff",sep="-"),res = 300, height = 2000, width = 3000)
#goplot(go.CC,showCategory = 5,layout = "sugiyama",geom = "text", title = paste(prefix,"GO-CC",sep="-"))
#dev.off()
##MF/
#tiff(file=paste(prefix,"GOMF_goplot.tiff",sep="-"),res = 300, height = 2000, width = 3000)
#goplot(go.MF,showCategory = 5,layout = "sugiyama",geom = "text",title = paste(prefix,"GO-MF",sep="-"))
#dev.off()



#KEGG dotplot barplot emapplot

#dotplot
#ek<- enrichKEGG(genelist, pvalueCutoff=0.05)
print("===================================KEGG===============================")
ek <- enrichKEGG(genelist, organism = 'hsa', keyType = 'kegg', pvalueCutoff = 1,pAdjustMethod = 'BH', minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2)

#ek<- enrichKEGG(genelist, pvalueCutoff=0.05, OrgDb = org.Hs.eg.db)
tiff(file=paste(prefix,"KEGG_dotplot.tiff",sep="-"),res = 300, height = 2500, width = 2400)
dotplot(ek,
        showCategory = 10,
        split = NULL,  title = "KEGG")
dev.off()
#barplot
tiff(file=paste(prefix,"KEGG_barplot.tiff",sep="-"),res = 300, height = 2500, width = 2400)
barplot(ek,
        showCategory = 10,
        split = NULL,  title = "KEGG")
dev.off()
##emapplot
#tiff(file=paste(prefix,"KEGG_emapplot.tiff",sep="-"),res = 300, height = 4000, width = 4000)
#emapplot(ek,
#         showCategory = 10,
#         split = NULL, font.size = 12,title = "KEGG")
#dev.off()
