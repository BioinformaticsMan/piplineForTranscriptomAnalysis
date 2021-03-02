getwd()
setwd("E://AD//data//GSE118553//GSE118553_series_matrix.txt//anova//sorted//vennPlot")

library(VennDiagram)
#hubGene <- read.csv("hubGene.csv",header = T)
tchubGene <- read.csv("TC-hubGene.csv",header = T)
fchubGene <- read.csv("FC-hubGene.csv",header = T)
echubGene <- read.csv("EC-hubGene.csv",header = T)
cehubGene <- read.csv("CE-hubGene.csv",header = T)
tcHubGene <- as.character(tchubGene[,1])
fcHubGene <- as.character(fchubGene[,1])
ecHubGene <- as.character(echubGene[,1])
ceHubGene <- as.character(cehubGene[,1])


# venn_list <- list(tcHubGene,fcHubGene,ecHubGene,ceHubGene)
# names(venn_list) <- c('tcHubGene','fcHubGene','ecHubGene','ceHubGene')
# 
# venn.diagram(venn_list, filename = 'hubGene_vennPlot.png', imagetype = 'png', margin = 0.2, 
#              alpha = 0.30, col = 'black', 
#              fill = c('red', 'blue','yellow','green'),
#              cex = 1, fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif')
# 
# 



tcDiff <- read.csv("TC_DEGs_limma.csv",header = T)
fcDiff <- read.csv("FC_DEGs_limma.csv",header = T)
ecDiff <- read.csv("EC_DEGs_limma.csv",header = T)
ceDiff <- read.csv("CE_DEGs_limma.csv",header = T)

tcdiff <- as.character(tcDiff[,1])
fcdiff <- as.character(fcDiff[,1])
ecdiff <- as.character(ecDiff[,1])
cediff <- as.character(ceDiff[,1])

# venn_list <- list(tcdiff,fcdiff,ecdiff,cediff)
# names(venn_list) <- c('tcDiff','fcDiff','ecDiff','ceDiff')
# venn.diagram(venn_list, filename = 'DEGs_vennPlot.png', imagetype = 'png', margin = 0.2, 
#              alpha = 0.30, col = 'black', 
#              fill = c('red', 'blue','yellow','green'),
#              cex = 1, fontfamily = 'serif', cat.cex = 1, cat.fontfamily = 'serif')
# 


tcModuleGene <- read.csv("TC_moduleGene.csv",header = T)
fcModuleGene <- read.csv("FC_moduleGene.csv",header = T)
ecModuleGene <- read.csv("EC_moduleGene.csv",header = T)
ceModuleGene <- read.csv("CE_moduleGene.csv",header = T)

tc_moduleGene <- as.character(tcModuleGene[,1])
fc_moduleGene <- as.character(fcModuleGene[,1])
ec_moduleGene <- as.character(ecModuleGene[,1])
ce_moduleGene <- as.character(ceModuleGene[,1])


# overlap_moduleGene = Reduce(intersect,list(tc_moduleGene,fc_moduleGene,ec_moduleGene,ce_moduleGene))
# overlap_diffGene = Reduce(intersect,list(tcdiff,fcdiff,ecdiff,cediff))
# 
# overlap_diffGene %in%overlap_moduleGene
# 
# write.csv(overlap_moduleGene,file="overlap_moduleGene.csv")

library(UpSetR)
library(dplyr)
library(tidyr)


listinput <- list(tc_moduleGene = tc_moduleGene,
                  fc_moduleGene = fc_moduleGene,
                  ec_moduleGene = ec_moduleGene,
                  ce_moduleGene = ce_moduleGene)

tiff(filename = "moduleGene_venn.tiff",res=300,width = 2000,height = 2000)
p <- upset(fromList(listinput),nsets = 4, order.by = "freq")
p
dev.off()


listinput <- list(tcdiff = tcdiff,
                  fcdiff = fcdiff,
                  ecdiff = ecdiff,
                  cediff = cediff)

tiff(filename = "DEGS_venn.tiff",res=300,width = 2000,height = 2000)
p <- upset(fromList(listinput),nsets = 4, order.by = "freq")
p
dev.off()


listinput <- list(tcHubGene = tcHubGene,
                  fcHubGene = fcHubGene,
                  ecHubGene = ecHubGene,
                  ceHubGene = ceHubGene)

tiff(filename = "hubGene_venn.tiff",res=300,width = 2000,height = 2000)
p <- upset(fromList(listinput),nsets = 4, order.by = "freq")
p
dev.off()




listinput <- list(tcHubGene = tcHubGene,
                  fcHubGene = fcHubGene,
                  ecHubGene = ecHubGene,
                  ceHubGene = ceHubGene,
                  tcdiff = tcdiff,
                  fcdiff = fcdiff,
                  ecdiff = ecdiff,
                  cediff = cediff)

tiff(filename = "DEGs_hubGene_venn.tiff",res=300,width = 2000,height = 2000)
p <- upset(fromList(listinput),nsets = 8, order.by = "freq")
p
dev.off()



listinput <- list(tc_moduleGene = tc_moduleGene,
                  fc_moduleGene = fc_moduleGene,
                  ec_moduleGene = ec_moduleGene,
                  ce_moduleGene = ce_moduleGene,
                  tcdiff = tcdiff,
                  fcdiff = fcdiff,
                  ecdiff = ecdiff,
                  cediff = cediff)

tiff(filename = "DEGs_moduleGene_venn.tiff",res=300,width = 2500,height = 2000)
p <- upset(fromList(listinput),nsets = 8, order.by = "freq")
p
dev.off()