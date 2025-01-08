Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
### 
library(dplyr)
library(tidyverse)
library(data.table)
##
getwd()
setwd("~/Downloads/Supple.Figure/Supple.Figure_2i")
###########################################################################
rawcounts <- read.csv("../Figure_1d/BMDM_RawCounts_AllSamples.csv")
View(rawcounts)
rawcounts<- rawcounts[, c(1,6:12)]
dim(rawcounts) # [1] 51732     8
colnames(rawcounts)
#[1] "GeneName" "LPS.1"    "LPS.2"    "LPS.3"    "LPS.4"    "ko_LPS.1" "ko_LPS.2" "ko_LPS.3"  
colnames(rawcounts) <- c("GeneName","LPS.1","LPS.2","LPS.3","LPS.4","ko_LPS.1","ko_LPS.2","ko_LPS.3")
View(rawcounts)

## select distinct Gene Name row.
rawcounts.1 <- rawcounts %>%
    distinct(GeneName, .keep_all = T)
dim(rawcounts.1) # [1] 51732     8
write.csv(rawcounts.1, file = "BMDM_wtlps_kolps_rawcounts.csv", row.names = FALSE)
View(rawcounts.1)
### genes names error : march1 should be Marchf1.
###########################################################################
BMDM <- read.csv("BMDM_wtlps_kolps_rawcounts.csv")
dim(BMDM) # [1] 51732    8 
View(BMDM)
metadata <- read.csv("BMDM_wtlps_kolps_metadata.csv", stringsAsFactors = T)
str(metadata)
View(metadata)
# 'data.frame':	7 obs. of  4 variables:
# $ id     : Factor w/ 7 levels "ko_LPS.1","ko_LPS.2",..: 4 5 6 7 1 2 3
# $ dex    : Factor w/ 1 level "LPS": 1 1 1 1 1 1 1 1
# $ treated: Factor w/ 2 levels "Mcoln1_KO_LPS",..: 2 2 2 2 1 1 1
# $ cell   : Factor w/ 1 level "BMDM": 1 1 1 1 1 1 1 1
class(metadata)
### Compete differentially expressed genes between Mcoln1 knock out and Wild type BMDM cells subjected to LPS
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = BMDM, 
                              colData = metadata, 
                              design=~treated, 
                               tidy=TRUE)
 
################################################
rld <- rlog(dds, blind = FALSE)
tail(assay(rld), 3)
########### Samples distance ###########################
sampleDists <- dist(t(assay(rld)))
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
############ differentially expressed genes between wild type and knockout Mcoln1 BMDM cells subjected to LPS###################
dds <- DESeq(dds) 
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

################## Results #####################################################
res1 <- results(dds,contrast = c("treated","Mcoln1_KO_LPS","WT_LPS"),
                tidy=TRUE)
View(res1)

#### save the result as the following file, which indicates the differential gene expression: 
write.csv(res1, file="Mcoln1_wtlps_kolps_DESeq2_deg.csv", row.names = FALSE)
################################################################################
## fix problem of the first column name
result <- read.csv("Mcoln1_wtlps_kolps_DESeq2_deg.csv")
View(result)
colnames(result)
# [1] row"            "baseMean"       "log2FoldChange"
# [5] "lfcSE"          "stat"           "pvalue"         "padj" 
colnames(result) <- c("GeneName", "baseMean","log2FoldChange", "lfcSE","stat","pvalue", "padj")
View(result)

result <- na.omit(result)
write.csv(result, file = "Mcoln1_wtlps_kolps_DESeq2_deg.clean.csv")
#################################################################################
### Rename the result 
deg_BMDMLPS <- result
######### Transfer mouse gene names to human gene names #############
deg_BMDMLPS$GeneName <- toupper(deg_BMDMLPS$GeneName)

### build up genes expression lists 
select.baseMean <- (deg_BMDMLPS$baseMean > 10)
table(select.baseMean)
# select.baseMean
# FALSE  TRUE 
# 5248 14874  
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
##############  id transform ##################################
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
# 1   A3GALT2     2.50582897 2.388406e-08   127550
# 2   A4GALT    -0.02581124 9.058730e-01    53947
# 3     AAAS    -0.73408150 1.826555e-05     8086

dim(exp.fc.id) #[1] 10809     4
exp.fc.id=na.omit(exp.fc.id) 
View(exp.fc.id) 

################ Arrange the gene list by the value of log2FoldChange #########
exp.fc.id.sorted <- exp.fc.id[order(exp.fc.id$log2FoldChange, decreasing = T),]
head(exp.fc.id.sorted, 2)

#      GeneName log2FoldChange       pvalue ENTREZID
# 7771   RNF125       11.10436 1.161329e-26    54941
# 4231     IL33       10.15596 8.979809e-09    90865

id.fc <- exp.fc.id.sorted$log2FoldChange  
names(id.fc) <- exp.fc.id.sorted$ENTREZID

head(id.fc,2)
#    54941    90865 
# 11.10436 10.15596 

View(id.fc)
gsea <- gseKEGG(id.fc, organism = "hsa",pvalueCutoff = 0.5)

### Check gseKEGG result ###
gseKEGG <- gsea@result
write.table(gseKEGG,"deg_wtlps_kolps_enrich_gseKEGG.txt",sep = "\t",col.names = NA)
View(gseKEGG)

gseaplot2(gsea, geneSetID = "hsa05321",
          title = "Inflammatory bowel disease",
          color = "black",
          base_size = 11,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = TRUE,
          ES_geom = "line")

#######################################################################################
### Arrange the genes list by score 
### Reference: https://www.baderlab.org/CancerStemCellProject/VeroniqueVoisin/AdditionalResources/GSEA
### score =(-log10(pValue)*SIGN(log2Foldchange))
View(exp.fc.id)
exp.fc.id.score <- exp.fc.id
View(exp.fc.id.score)
exp.fc.id.score$score = (-log10(exp.fc.id.score$pValue)*sign(exp.fc.id.score$log2Foldchange))

### Error in log10(exp.fc.id.score$pValue) :  non-numeric argument to mathematical function
### Either manually construct score column using excel 
### Or sign in a relevant pvalue in order to construct the score value for ranking in this dataset

exp.fc.id.score <- exp.fc.id.score[order(exp.fc.id.score$pvalue, decreasing = F),]
head(exp.fc.id.score, 25)

### sign in a relevant log value into "0" in the first 22 genes name, due to the log10(0) gives a N/A value.  
exp.fc.id.score$pvalue[c(1:22)] <- rep(c(6.0e-306))
head(exp.fc.id.score, 25)
####
library(dbplyr)
library(data.table)
exp.fc.id.score$score = (-log10(exp.fc.id.score$pvalue)*sign(exp.fc.id.score$log2FoldChange))

################ Arrange the gene list by the value of score  #########
exp.fc.id.score.sorted <- exp.fc.id.score[order(exp.fc.id.score$score, decreasing = T),]
head(exp.fc.id.score.sorted, 3)
#     GeneName log2FoldChange        pvalue ENTREZID    score
# 1424     CD83       5.821350 1.460907e-306     9308 305.8354
# 137     ACTN1       3.247347 6.000000e-306       87 305.2218
# 208    ADGRL2       3.956688 6.000000e-306    23266 305.2218

id.fc.score <- exp.fc.id.score.sorted$score 
names(id.fc.score) <- exp.fc.id.score.sorted$ENTREZID
head(id.fc.score,2)
#     9308       87 
# 305.8354 305.2218 

View(id.fc.score)
gsea <- gseKEGG(id.fc.score, organism = "hsa",pvalueCutoff = 0.5)
?gseKEGG

### check the gseKEGG result ###
gseKEGG <- gsea@result

write.table(gseKEGG,"deg_wtlps_kolps_enrich_score_gseKEGG.txt",sep = "\t",col.names = NA)
View(gseKEGG)

gseaplot2(gsea, geneSetID = "hsa05321",
          title = "Inflammatory bowel disease",
          color = "black",
          base_size = 11,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = TRUE,
          ES_geom = "line")

### Alternatively using the "result" file to build GSEAPreranked files, and  
### using GSEA software (Broad Institue GSEA v4.1.3 and GSEA v4.3.3 gene set enrichment analysis)for GSEA analysis. 
################################################################################
sessionInfo()

# R version 4.3.3 (2024-02-29 ucrt)
# R version 4.1.3 (2022-03-10)
# other attached packages:
# [1] dbplyr_2.5.0           org.Hs.eg.db_3.18.0    enrichplot_1.22.0      clusterProfiler_4.10.1
# [5] GSEABase_1.64.0        graph_1.80.0           annotate_1.80.0        XML_3.99-0.16.1       
# [9] AnnotationDbi_1.64.1   IRanges_2.36.0         S4Vectors_0.40.2       Biobase_2.62.0        
# [13] BiocGenerics_0.48.1    msigdbr_7.5.1          data.table_1.15.2      lubridate_1.9.3       
# [17] forcats_1.0.0          stringr_1.5.1          purrr_1.0.2            readr_2.1.5           
# [21] tidyr_1.3.1            tibble_3.2.1           ggplot2_3.5.0          tidyverse_2.0.0       
# [25] dplyr_1.1.4           

