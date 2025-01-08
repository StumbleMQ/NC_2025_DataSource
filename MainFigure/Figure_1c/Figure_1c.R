###
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
###
library(dplyr)
library(tidyverse)
###
##setwd("~/Downloads/NC_2025_RNAseq_DataSource/MainFigure/Figure_1c")
### load the differential expressed genes wild type BMDM cells subjected to LPS 
### treated with and without ML-SA5. 
BMDMLPS_deg <- read.csv("BMDMLPS_DESeq2_deg.csv")
BMDMLPS_deg <- BMDMLPS_deg[, -1]
colnames(BMDMLPS_deg)
colnames(BMDMLPS_deg) <- c("GeneName",  "baseMean","log2FoldChange", "lfcSE",        
                           "stat","pvalue", "padj")
dim(BMDMLPS_deg) #[1] 51732     7
### clean the data set 
BMDMLPS_deg <- na.omit(BMDMLPS_deg)
dim(BMDMLPS_deg) 
View(BMDMLPS_deg) # [1] 19426     7
### save the cleaned data set
write.csv(BMDMLPS_deg, file= "Mcoln1_wt_lps_DESeq2_deg_gene_cleaned.csv",
          row.names = FALSE)
###########################################################################################
data <- read.csv("Mcoln1_wt_lps_DESeq2_deg_gene_cleaned.csv")
row.names(data) <- data$GeneName
View(data)
## Volcano plot for Figure 1c
library(dplyr)
library(ggplot2)                         
library(ggrepel)                               
data$significance="stable"
data$significance[data$log2FoldChange>=1 & data$padj <0.01]="up"
data$significance[data$log2FoldChange<= -1 & data$padj <0.01]="down"

sum(data$significance == "up") # 542
sum(data$significance == "down") # 622
sum(data$significance == "stable") # 18262

### demonstrate several selected genes 
data_selected.c <- data[(data$GeneName == "Il1b"|data$GeneName =="Nfkb1"|data$GeneName =="Il6"|data$GeneName =="Tnf"|data$GeneName =="Slc7a11"),]
View(data_selected.c)

data_selected.stable <- data[(data$GeneName =="Il6"|data$GeneName =="Tnf"|data$GeneName =="Nfkb1"),]
View(data_selected.stable)
data_selected.down <- data[(data$GeneName == "Il1b"),]
View(data_selected.down)
data_selected.up <- data[(data$GeneName =="Slc7a11"),]
View(data_selected.up)
### Sign in a log value into "0" in data_selected.up in order to show the data point in graph:
data_selected.up$padj <- c(6.0e-298)

p1 = ggplot(data,aes(log2FoldChange,-1*log10(padj)))+
     xlim(-6,6)+ylim(0,300)+ 
   geom_point(aes(color=significance), shape=16, size=1)+  #, alpha=0.7
   #geom_point(data=data_selected.up, fill = "firebrick",color = "firebrick", shape=21, size =2)+
   geom_point(data=data_selected.stable, color = "black",fill="gray50", shape=21, size=3)+
   geom_point(data=data_selected.up, color = "firebrick3",fill="firebrick3", shape=21, size=3)+
   geom_point(data=data_selected.down, color = "blue",fill="blue", shape=21, size=3)+
   theme_classic()+
   scale_color_manual(values = c("blue","gray50","firebrick3"))+
   geom_label_repel(data=data_selected.down,
                   family = "Times New Roman", 
                   fontface =4, 
                   color = "dark blue", 
                 aes(label=rownames(data_selected.down))) + 
 
   geom_label_repel(data=data_selected.stable,
                   family = "Times New Roman", 
                   fontface =4, 
                   color ="black",size = 3, 
                   aes(label=rownames(data_selected.stable)),
                   nudge_x = .15,
                   box.padding= 0.5,
                   nudge_y = 1,
                   segment.curvature = -0.1,
                   segment.ncp =3, 
                   segment.angle =20) + 
   geom_label_repel(data=data_selected.up,
                   family = "Times New Roman", 
                   fontface =4, 
                   color ="firebrick3",
                   aes(label=rownames(data_selected.up))) + 
   geom_hline(yintercept = 2,linetype=4,size=0.3)+
   geom_vline(xintercept = c(-1,1),linetype=4,size=0.3)+
   theme(title=element_text(size = 14),text = element_text(size=14))+
   labs(x="log2(FoldChange)",y="-log10(padj)") 
 p1
################## End ####################
