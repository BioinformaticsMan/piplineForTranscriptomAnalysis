
args <- commandArgs(T)

#getwd()
# setwd("/home/tonyleao/wkd/AD_manuscript_revise/WGCNA_PC1/")

## Use all samples and output the network of module gene and filter gene.


library(WGCNA)
library(data.table)
library(stringr)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

#############################################################
## 01 load data

# datExpr0 <- read.table("GSE118553_series_matrix_TC_merge_small_test.txt",sep="\t",header = T,row.names = 1);
datExpr0 <- read.table(args[1],sep="\t",header = T,row.names = 1);
datExpr0 <- t(log2(1+datExpr0))
head(datExpr0[1:5,1:5])
datExpr <- datExpr0

# datTraits <- read.table("GSE118553_series_matrix_trait.txt",header = T,row.names = 1)
datTraits <- read.table(args[2],header = T,row.names = 1)
index = c(rownames(datExpr))
datTraits <- as.data.frame(t(datTraits[,index]))

#save(datExpr,datTraits,file = "FPKM-01-dataInput.RData")

#############################################################
## 02 check data
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK



##############################################################
### 03 detect outliers of datasets
#sampleTree = hclust(dist(datExpr),method = "average")
#tiff(filename = "Sample clustering to detect outliers.tiff",res=300,width = 3000,height = 2000)
#plot(sampleTree,main="Sample clustering to detect outliers",sub="",xlab="")
#abline(h=29,col="red")
#dev.off()
#
#clust = cutreeStatic(sampleTree, cutHeight = 29, minSize=1)
#table(clust)
#keepSamples = (clust==1)
#table(keepSamples)
#datExpr = datExpr[keepSamples,]
#nGenes = ncol(datExpr)
#nSamples = nrow(datExpr)
#dim(datExpr)
##save(datExpr,file = "FPKM-03-dataInput.RData")
#datTraits = datTraits[keepSamples,]
#dim(datTraits)
#
#sampleTree = hclust(dist(datExpr),method = "average")
#tiff(filename = "Sample clustering to after removing outliers.tiff",res=300,width = 3000,height = 2000)
#plot(sampleTree,main="Sample clustering to after removing outliers",sub="",xlab="")
#dev.off()




#############################################################
## 03 Choose a set of soft-thresholding powers
# Call the network topology analysis function
#type="signed"
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,networkType="signed")

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
tiff(filename = "power.tiff",res=300,width = 2000,height = 2000)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

power <- sft$powerEstimate
softPower = power
print("power is :   ")
print (power)


#############################################################
## 05 check weather the network based on soft-thresholding powers close to scale free

# ADJ1_cor <- abs(WGCNA::cor( multiExpr[[1]]$data,use = "p" ))^softPower
## for small number of genes(<5000)
## k <- as.vector(apply(ADJ1_cor,2,sum,na.rm=T))

## for large number of genes??
k <- softConnectivity(datE=datExpr,power=softPower) 
sizeGrWindow(10, 5)
tiff(filename = "scale_free.tiff",res=300,width = 2000,height = 2000)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()




#############################################################
## 04 contruct adjacency matrix and toplogy matrix based on soft-thresholding powers

##1) Co-expression similarity and adjacency
adjacency = adjacency(datExpr, power = softPower)

##2) Turn adjacency into topological overlap
type = "signed"
TOM = TOMsimilarity(adjacency,TOMType = type)
dissTOM = 1-TOM
##3) call the hierachical clustering function
hierTOM = hclust(as.dist(dissTOM),method="average")
geneTree = hclust(as.dist(dissTOM),method = "average")





#############################################################
## 06 cluster gene of topology matrix based on dissimilarity between genes
## Then use dynamic tree cut to identify module

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
tiff(filename = "geneTree.tiff",res=300,height = 3000,width = 3000)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()



# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
tiff(filename = "plotDendroAndColors1.tiff",res=300,width = 2000,height = 2000)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()




#############################################################
# 08 Calculate eigengenes of every module and use the gigengenes as first
## components of genes, which represents the entire level of gene expression
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
#MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
#sizeGrWindow(7, 6)
tiff(filename = "METree.tiff",res=300,height = 3000,width = 3000)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
## chose 75% relationship to merge
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

## plot eigengene adjacency heatmap
tiff(filename = "Eigengene adjacency heatmap.tiff",res=300,height = 3000,width = 3000)
plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90) 
dev.off()


#############################################################
# 07 Network heatmap plot, selected genes 
## randomly choose 400 genes and draw topology overlap heatmap.
## deep yellow and red represent high level of topology overlap

# draw topology
nSelect = 400 
# For reproducibility, we set the random seed 
set.seed(10)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster. 
selectTree = hclust(as.dist(selectTOM), method = "average") 
selectColors = dynamicColors[select] 
# Open a graphical window 
sizeGrWindow(9,9) 
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot 
tiff(filename = "Network heatmap plot, selected genes.tiff",res=300,width = 2000,height = 2000)
plotDiss = selectTOM^softPower
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes") 
dev.off()



#######################################################################################
## 09 merge the module with high relationship (cor>0.8) or small dissimilarity(<0.2)

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

## plot DendroAndColors.tiff
sizeGrWindow(12, 9)
tiff(filename = "plotDendroAndColors2.tiff",res=300,height = 3000,width = 3000)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

##############################################################


#############################################################
## 10 draw sample dendrogram and trait heatmap
#sampleTree = hclust(dist(datExpr0), method = "average")
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
tiff(filename = "Sample dendrogram and trait heatmap.tiff",res=300,width =4000,height = 2000)
plotDendroAndColors(sampleTree2, 
                    traitColors,
                    groupLabels = names(datTraits),
                    rowTextAlignment = "right-justified",
                    addTextGuide = TRUE ,
                    hang = 0.03,
                    dendroLabels = NULL, 
                    addGuide = FALSE,  
                    guideHang = 0.05,
                    main = "Sample dendrogram and trait heatmap")
dev.off()



#############################################################
## 11 draw module-trait relationships heatmap

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


################################################################################
################################################################################
## here we output the eigengenes of the datExpr for further analysis (AD project)
# write.csv(MEs,paste(prefix,"mergedMEs.csv",sep = "_"))
MEs0_out = moduleEigengenes(datExpr, moduleColors)$eigengenes
write.csv(MEs0_out,"mergedMEs.csv")
## output the varExplained A dataframe in which each column corresponds to a module, 
## with the component varExplained[PC, module] giving the variance of module module 
varExplain0 <- moduleEigengenes(datExpr, moduleColors)$varExplained	
colnames(varExplain0) = colnames(MEs0_out)
write.csv(varExplain0,"varExplain0.csv")



###############################
sizeGrWindow(10,6)
# Will display correlations and their p-values
tiff(filename = "Module-trait relationships heatmap.tiff",res=300,height = 3000,width = 3000)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

#dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               cex.lab.x = 1,
               cex.lab.y = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()



#############################################################
## 12 select key module related to traits based on relationship between traits and eigengene
cor_ADR <- signif(WGCNA::cor(datTraits,mergedMEs,use="p",method="pearson"),5)
p.values <- corPvalueStudent(cor_ADR,nSamples=nrow(datTraits))
#diagnosis_MS_max_cor <- which.max(abs(cor_ADR["diagnosis_state",]))
#diagnosis_MS_max_p <- which.min(p.values["diagnosis_state",])
# Freq_MS_max_cor <- which.max(abs(cor_ADR["Freq",-which(colnames(cor_ADR) == "MEgrey")]))
# Freq_MS_max_p <- which.min(p.values["Freq",-which(colnames(p.values) == "MEgrey")])
postiveModule <-  which.max(cor_ADR["diagnosis_state",])
negativeModule <-  which.min(cor_ADR["diagnosis_state",])
pModuleName <- colnames(MEs[postiveModule])
nModuleName <- colnames(MEs[negativeModule])
#extract the color
pModuleColor <- substr(pModuleName,3,nchar(pModuleName))
nModuleColor <- substr(nModuleName,3,nchar(nModuleName))

print (paste("Positive correlation: ",pModuleName,"cor: ",cor_ADR[1,pModuleName],sep=" "))
print (paste("Negative correlation: ",nModuleName,"cor: ",cor_ADR[1,nModuleName],sep=" "))


#############################################################
## 13 identify hub genes by gene network significance
## use diagnosis
GS1 <- as.numeric(WGCNA::cor(datTraits[,1],datExpr,use="p",method="pearson"))
GeneSignificance <- abs(GS1)
## obtained the significance of traits in every module
ModuleSignificance <- tapply(GeneSignificance,mergedColors,mean,na.rm=T)




#############################################################
## 14 identify hub genes

## calculate the correlation matrix between trait and gene
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$diagnosis_state)
#weight
colnames(weight) = "Diagnosis_state"

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
colnames(geneTraitSignificance) = paste("GS.", colnames(weight), sep="")
colnames(GSPvalue) = paste("p.GS.", colnames(weight), sep="")


## calculate correlation matrix between module and gene
# names (colors) of the modules
modNames = substring(colnames(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
colnames(geneModuleMembership) = paste("MM", modNames, sep="")
colnames(MMPvalue) = paste("p.MM", modNames, sep="")


####output GS and MS file for positive module and negative module
## for positive module
module = pModuleColor
pheno = "Diagnosis_state"
column = match(module, modNames)
moduleGenes = moduleColors==module

MM <- abs(geneModuleMembership[moduleGenes,column])
GS <- abs(geneTraitSignificance[moduleGenes,1])
save(MM,GS,file="FC_MM_GC.RData")

sizeGrWindow(7, 7)
tiff(filename = paste(pModuleName,pheno,"Module membership vs gene significance.tiff",sep="_"),res=300,width = 2000,height = 2000)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()


## select hub gene  abs(MM) >median GS >median
#FilterGenes = (geneTraitSignificance > 0.5) & (abs(geneModuleMembership[,paste("MM",module,sep="")]) >0.5)
#FilterGenes = (abs(geneTraitSignificance)) > mean(abs(geneTraitSignificance[moduleGenes,1])) & (abs(geneModuleMembership[,paste("MM",module,sep="")])) > mean(abs(geneModuleMembership[,paste("MM",module,sep="")]))
FilterGenes = (abs(geneTraitSignificance) >= 0.2) & (abs(geneModuleMembership[,paste("MM",module,sep="")]) >=0.8)
table(FilterGenes)
FilterGenes = moduleGenes &FilterGenes
table(FilterGenes)


print("================write====================")
write.csv(abs(geneModuleMembership[moduleGenes,]),paste(pModuleName,pheno,"MM.csv",sep="_"))
write.csv(abs(geneTraitSignificance[moduleGenes,]),paste(pModuleName,pheno,"GS.csv",sep="_"))
write.csv(FilterGenes,paste(pModuleName,pheno,"Filters.csv",sep="_"))

## filter gene network of positive module
trait_filterGenes <- colnames(datExpr)[FilterGenes]
filterGene_TOM <- TOM[FilterGenes,FilterGenes]
dimnames(filterGene_TOM) = list(colnames(datExpr)[FilterGenes], colnames(datExpr)[FilterGenes])
cyt = exportNetworkToCytoscape(filterGene_TOM,                        ## toplogy overlap of gene
                              edgeFile = paste(pModuleName,pheno,"filterGenes","CytoscapeInput-edges.txt",sep="_"),
                              nodeFile = paste(pModuleName,pheno,"filterGenes","CytoscapeInput-nodes.txt",sep="_"),
                              weighted = TRUE,
                              threshold = 0.05,
                              nodeNames = trait_filterGenes,          ### gene name(based on what you want to extract)
                              nodeAttr = moduleColors[FilterGenes])### annoted the color of filterGenes


# output positive module gene network
module = pModuleColor
moduleGenes = moduleColors==module

trait_moduleGenes <- colnames(datExpr)[moduleGenes]
moduleGene_TOM <- TOM[moduleGenes,moduleGenes]
dimnames(moduleGene_TOM) = list(colnames(datExpr)[moduleGenes],colnames(datExpr)[moduleGenes])
cyt = exportNetworkToCytoscape(moduleGene_TOM,                        ## toplogy overlap of gene
                               edgeFile = paste(pModuleName,pheno,"moduleGenes","CytoscapeInput-edges.txt",sep="_"),
                               nodeFile = paste(pModuleName,pheno,"moduleGenes","CytoscapeInput-nodes.txt",sep="_"),
                               weighted = TRUE,
                               threshold = 0.05,
                               nodeNames = trait_moduleGenes,          ### gene name(based on what you want to extract)
                               nodeAttr = moduleColors[moduleGenes])### annoted the color of module gene
print("================writeOver================")




## for negative module
module = nModuleColor
pheno = "Diagnosis_state"
column = match(module, modNames)
moduleGenes = moduleColors==module

sizeGrWindow(7, 7)
tiff(filename = paste(nModuleName,pheno,"Module membership vs gene significance.tiff",sep="_"),res=300,width = 2000,height = 2000)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

## select hub gene  abs(MM) >0.5 GS >0.5
# FilterGenes_n = (abs(geneTraitSignificance)) > mean(abs(geneTraitSignificance[,1])) & (abs(geneModuleMembership[,paste("MM",module,sep="")])) > mean(abs(geneModuleMembership[,paste("MM",module,sep="")]))
# table(FilterGenes_n)

FilterGenes = (abs(geneTraitSignificance) >= 0.2) & (abs(geneModuleMembership[,paste("MM",module,sep="")]) >=0.8)
table(FilterGenes)
FilterGenes = moduleGenes & FilterGenes
table(FilterGenes)


print("================write====================")
write.csv(abs(geneModuleMembership[moduleGenes,]),paste(nModuleName,pheno,"MM.csv",sep="_"))
write.csv(abs(geneTraitSignificance[moduleGenes,]),paste(nModuleName,pheno,"GS.csv",sep="_"))
write.csv(FilterGenes,paste(nModuleName,pheno,"Filters.csv",sep="_"))


## output filter gene network of negative module
trait_filterGenes <- colnames(datExpr)[FilterGenes]
filterGene_TOM <- TOM[FilterGenes,FilterGenes]
dimnames(filterGene_TOM) = list(colnames(datExpr)[FilterGenes], colnames(datExpr)[FilterGenes])
cyt = exportNetworkToCytoscape(filterGene_TOM,                        ## toplogy overlap of gene
                              edgeFile = paste(nModuleName,pheno,"filterGenes","CytoscapeInput-edges.txt",sep="_"),
                              nodeFile = paste(nModuleName,pheno,"filterGenes","CytoscapeInput-nodes.txt",sep="_"),
                              weighted = TRUE,
                              threshold = 0.05,
                              nodeNames = trait_filterGenes,          ### gene name(based on what you want to extract)
                              nodeAttr = moduleColors[FilterGenes])### annoted the color of filterGenes

# output negative module gene network
module = nModuleColor
moduleGenes = moduleColors==module

trait_moduleGenes <- colnames(datExpr)[moduleGenes]
moduleGene_TOM <- TOM[moduleGenes,moduleGenes]
dimnames(moduleGene_TOM) = list(colnames(datExpr)[moduleGenes],colnames(datExpr)[moduleGenes])
cyt = exportNetworkToCytoscape(moduleGene_TOM,                        ## toplogy overlap of gene
                               edgeFile = paste(nModuleName,pheno,"moduleGenes","CytoscapeInput-edges.txt",sep="_"),
                               nodeFile = paste(nModuleName,pheno,"moduleGenes","CytoscapeInput-nodes.txt",sep="_"),
                               weighted = TRUE,
                               threshold = 0.05,
                               nodeNames = trait_moduleGenes,          ### gene name(based on what you want to extract)
                               nodeAttr = moduleColors[moduleGenes])### annoted the color of module gene




print("================writeOver================")




#############################################################
## 15 output hub genes for Cytoscape
# Export the network into edge and node list files Cytoscape can read
##all network
# probes = colnames(datExpr)
# head(probes)
# dimnames(TOM) <- list(probes, probes)
# head(moduleColors)
# cyt = exportNetworkToCytoscape(TOM,                        ## topology overlap for all gene
#                                edgeFile = "CytoscapeInput-edges.txt",
#                                nodeFile = "CytoscapeInput-nodes.txt",
#                                weighted = TRUE,
#                                threshold = 0.05,           ##
#                                nodeNames = probes,         ## the name of the gene
#                                nodeAttr = moduleColors)    ## module color of the gene

