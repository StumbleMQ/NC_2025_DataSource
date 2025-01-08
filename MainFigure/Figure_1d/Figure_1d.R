Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())

### 
library(dplyr)
library(tidyverse)
library(data.table)
###
getwd()
# setwd("~/Downloads/NC_figure_20250102/Figure_1d")
###########################################################################
### The file named "BMDM2024_gseaFigures.csv" contains annotated genes of all samples 
### BMDM cells collected from wild type Mcoln1 mice subjected to LPS with ML-SA5 treatment group: "LPSA5.1" "LPSA5.2" "LPSA5.3" "LPSA5.4"
### BMDM cells collected from wild type Mcoln1 mice subjected to LPS without ML-SA5 treatment group::"LPS.1", "LPS.2", "LPS.3", "LPS.4"  
rawcounts <- read.csv("BMDM_RawCounts_AllSamples.csv")
rawcounts<- rawcounts[, c(1:9)]
dim(rawcounts) # [1] 51732     9
colnames(rawcounts)
#[1] "GeneName"       "LPSA5.1" "LPSA5.2" "LPSA5.3" "LPSA5.4" "LPS.1" "LPS.2" "LPS.3"   "LPS.4"  
colnames(rawcounts) <- c("GeneName","LPSA5.1", "LPSA5.2", "LPSA5.3", "LPSA5.4", "LPS.1","LPS.2","LPS.3","LPS.4")
View(rawcounts)

## select distinct Gene Name row.
rawcounts.1 <- rawcounts %>%
    distinct(GeneName, .keep_all = T)
dim(rawcounts.1) # [1] 51732     9

### 
write.csv(rawcounts.1, file = "BMDMLPS_rawcounts.csv", row.names = FALSE)
### genes names error : march1 should be Marchf1.
###########################################################################
BMDM <- read.csv("BMDMLPS_rawcounts.csv")
dim(BMDM) # [1] 51732    9 ## fix the first column problem
### meta data file for analysis
metadata <- read.csv("BMDMLPS_metadata.csv", stringsAsFactors = T)
str(metadata)
# 'data.frame':	8 obs. of  4 variables:
# $ id     : Factor w/ 8 levels "LPS.1","LPS.2",..: 5 6 7 8 1 2 3 4
# $ dex    : Factor w/ 1 level "LPS": 1 1 1 1 1 1 1 1
# $ treated: Factor w/ 2 levels "Control","MLSA5": 2 2 2 2 1 1 1 1
# $ cell   : Factor w/ 1 level "BMDM": 1 1 1 1 1 1 1 1
class(metadata)
### Differential Analysis 
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = BMDM, 
                              colData = metadata, 
                              design=~treated, 
                               tidy=TRUE)
 
###############################################
### samples distances
rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
###
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(8, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
## PCA Analysis
plotPCA(rld,intgroup=c("treated"))
############ Differential analysis ######################################
dds <- DESeq(dds) 
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

################## ###################
res1 <- results(dds,tidy=TRUE)
View(res1)
#### BMDMLPS_DESeq2_deg.csv indicates the differential genes expression: 
write.csv(res1, file="BMDMLPS_DESeq2_deg.csv")

#################################################################
dds.n <- rlogTransformation(dds)  
# After DESeq2 Analysis, perform data set normalization. 
exprSet_new=assay(dds.n)
View(exprSet_new)
dim(exprSet_new) # 51732 8
###############################################################
par(cex = 0.7)
n.sample=ncol(BMDM)
cols <- rainbow(n.sample*1.2)
colnames(BMDM)
c1 = column_to_rownames(BMDM,"GeneName")
#### counts without normalization ######
boxplot(c1, col = cols,main="expression value",las=2)
#### counts after normalization #########
boxplot(exprSet_new, col = cols,main="expression value",las=2)
write.csv(exprSet_new, file = "exprSetrlog_BMDM_DESeq2.csv")
### this file was used for building the heat map graph Figure 4a
##############################################################
# after DESeq differential gene expression analysis,
# remove rowSums !=0 and use this for heat-map visualization. 
dds2 <- dds[rowSums(counts(dds)) != 0, ]
dds.n.c <- rlogTransformation(dds2) 
exprSet_new.c=assay(dds.n.c)
View(exprSet_new.c)
dim(exprSet_new.c) # 26895  8 
boxplot(exprSet_new.c, col = cols,main="expression value",las=2)
write.csv(exprSet_new.c, file = "exprSetrlog_BMDM_DESeq2.clean.csv")
### this file also can be used for building the heat map graph of Figure 4a.
################################################################################
### compete for Figure_1d:
deg_BMDMLPS <- read.csv("BMDMLPS_DESeq2_deg.csv")
View(deg_BMDMLPS)
## fix problem of the first column name
deg_BMDMLPS <- deg_BMDMLPS[ , -1]
colnames(deg_BMDMLPS)
# [1] row"            "baseMean"       "log2FoldChange"
# [5] "lfcSE"          "stat"           "pvalue"         "padj" 
colnames(deg_BMDMLPS) <- c("GeneName", "baseMean","log2FoldChange", "lfcSE","stat","pvalue", "padj")
View(deg_BMDMLPS)
## Transfer mouse gene names to human gene names
deg_BMDMLPS$GeneName <- toupper(deg_BMDMLPS$GeneName)

### build up genes expression lists 
select.baseMean <- (deg_BMDMLPS$baseMean > 10)
table(select.baseMean)
# select.baseMean
# FALSE  TRUE 
# 36536 15196 
degs.list=as.character(deg_BMDMLPS$GeneName)[select.baseMean]
res3=distinct(deg_BMDMLPS,GeneName,.keep_all = T)
View(res3)
rownames(res3)=res3$GeneName
### 
f1=res3[degs.list,]
View(f1)
## select log2FoldChange and pvalue as parameters to build the differential expressed genes list.
exp.fc=f1[,c(1,3,6)]
View(exp.fc)
##############  id转化##################################
library(msigdbr)
library(GSEABase)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
expm.id <- bitr(exp.fc$GeneName, 
                fromType = "SYMBOL", toType = "ENTREZID", 
                OrgDb = "org.Hs.eg.db")

head(expm.id)
View(expm.id)

exp.fc.id <- merge(exp.fc, expm.id,by.x="GeneName", by.y="SYMBOL", all=F)
head(exp.fc.id, 3)
# GeneName log2FoldChange       pvalue ENTREZID
# 1   A4GALT    -0.61720818 0.0006521589    53947
# 2     AAAS     0.06253631 0.5775995916     8086

dim(exp.fc.id) #[1] 10805     4

exp.fc.id=na.omit(exp.fc.id) 
View(exp.fc.id) 

## Arrange the genes list by log2FoldChange###################################
################ Arrange the gene list by the value of log2FoldChange #########
exp.fc.id.sorted <- exp.fc.id[order(exp.fc.id$log2FoldChange, decreasing = T),]
head(exp.fc.id.sorted, 2)
# GeneName log2FoldChange       pvalue ENTREZID
# 2141    CXCR1       5.133285 2.596725e-21     3577
# 5643   MYO18B       5.009726 3.387916e-32    84700
id.fc <- exp.fc.id.sorted$log2FoldChange  
names(id.fc) <- exp.fc.id.sorted$ENTREZID
head(id.fc,2)
#    3577    84700   
# 5.133285 5.009726 
View(id.fc)
gsea <- gseKEGG(id.fc, organism = "hsa",pvalueCutoff = 0.5)
?gseKEGG

#Reading KEGG annotation online: "https://rest.kegg.jp/link/hsa/pathway"...
#Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/hsa"...
#preparing geneSet collections...
#GSEA analysis...
#leading edge analysis...
# done...
# There were 15 warnings (use warnings() to see them)

#### Check Result
gseKEGG <- gsea@result
write.table(gseKEGG,"deg_BMDMLPS_enrich_gseKEGG.txt",sep = "\t",col.names = NA)
View(gseKEGG)

gseaplot2(gsea, geneSetID = "hsa04216",
          title = "Ferroptosis",
          color ="black",
          pvalue_table = TRUE)

gseaplot2(gsea, geneSetID = "hsa05321",
          title = "Inflammatory bowel disease",
          color = "black",
          base_size = 11,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = TRUE,
          ES_geom = "line")

?gseaplot2

#######################################################################################

## Arrange the genes list by score 
# score =(-log10(pValue)*SIGN(log2Foldchange))
View(exp.fc.id)
###############
exp.fc.id.score <- exp.fc.id
View(exp.fc.id.score)
exp.fc.id.score$score = (-log10(exp.fc.id.score$pValue)*sign(exp.fc.id.score$log2Foldchange))
# Error in log10(exp.fc.id.score$pValue) :  non-numeric argument to mathematical function
## manually construct score column using excel 

exp.fc.id.score <- exp.fc.id.score[order(exp.fc.id.score$pvalue, decreasing = F),]
head(exp.fc.id.score, 5)

### sign in a log vale into "0" in the first four gene name
exp.fc.id.score$pvalue[c(1:4)] <- c(6.0e-298,6.0e-298,6.0e-298,6.0e-298)

head(exp.fc.id.score, 10)

library(dbplyr)
library(data.table)

exp.fc.id.score$score = (-log10(exp.fc.id.score$pvalue)*sign(exp.fc.id.score$log2FoldChange))
################ Arrange the gene list by the value of log2FoldChange #########

exp.fc.id.score.sorted <- exp.fc.id.score[order(exp.fc.id.score$score, decreasing = T),]
head(exp.fc.id.score.sorted, 8)
#     GeneName log2FoldChange        pvalue ENTREZID    score
# 1756     CLIC5       3.006438 6.000000e-298    53405 297.2218
# 8588   SLC7A11       2.282536 6.000000e-298    23657 297.2218
# 6795     PLPP3       1.637385 3.692126e-297     8613 296.4327

id.fc.score <- exp.fc.id.score.sorted$score 
names(id.fc.score) <- exp.fc.id.score.sorted$ENTREZID
head(id.fc.score,2)
#      53405    23657 
#  297.2218 297.2218  

View(id.fc.score)
gsea <- gseKEGG(id.fc.score, organism = "hsa",pvalueCutoff = 0.5)
?gseKEGG

# preparing geneSet collections...
# GSEA analysis...
# leading edge analysis...
# done...
# There were 14 warnings (use warnings() to see them)

#### check result
gseKEGG <- gsea@result
write.table(gseKEGG,"deg_BMDMLPS_enrich_score_gseKEGG.txt",sep = "\t",col.names = NA)
View(gseKEGG)

gseaplot2(gsea, geneSetID = "hsa04216",
          title = "Ferroptosis",
          color ="black",
          pvalue_table = TRUE)
gseaplot2(gsea, geneSetID = "hsa05321",
          title = "Inflammatory bowel disease",
          color = "black",
          base_size = 11,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = TRUE,
          ES_geom = "line")

### Alternatively using the "result" file to build GSEAPreranked files 
### using GSEA software ( Broad Institue GSEA v4.1.3 and GSEA v4.3.3 gene set enrichment analysis)for GSEA analysis. 
################################################################################

