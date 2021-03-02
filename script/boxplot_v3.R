rm(list=ls())
#library(pacman)
#library(ggsignif)
#library(ggpubr)
library(ggplot2)

getwd()
setwd("~/wkd/AD_project/data/boxplot/boxplot/TC/")

d <- read.csv("TC_classfiers_metrics.csv",header = T,sep = ",")
prefix="TC"




classifier = d$classifie
F_1_score = as.numeric(d$F_1_score)
PPV = as.numeric(d$PPV)
NPV = as.numeric(d$NPV)
Sensitivity= as.numeric(d$Sensitivity)
Specificity = as.numeric(d$Specificity)

length(d$F_1_score)

value = c(d$F_1_score,d$PPV,d$NPV,d$Sensitivity,d$Specificity)
metrics <- c(rep("F_1_score",60),rep("PPV",60),rep("NPV",60),rep("Sensitivity",60),rep("Specificity",60))
color <- c(rep("#933848",60),rep("#A17271",60),rep("#A79181",60),rep("Sensitivity",60),rep("#63786E",60))
myData <- data.frame(Classifier = rep(d$classifier,5),Metrics = metrics, Value = value, Color = color)



tiff(filename = paste(prefix,"AD_classifiers_metrics.tiff",sep="_"),res=300,width = 3000,height = 2000)
ggboxplot(myData, 
          x="Classifier", 
          y="Value",
          fill="Metrics",
          xlab = "Classifiers",
          ylab = "Value",
          bxp.errorbar = TRUE,
          error.plot = "linerange",
          palette = "ucscgb")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1),
        legend.position = "top")
dev.off()


