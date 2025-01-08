Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
### 
library(dplyr)
library(tidyverse)
library(data.table)
## set a working directory
getwd()
setwd("~/Downloads/NC_2025_RNAseq_DataSource/Supple.Figure_1e")
## The following data set was downloaded from
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148505
## and the file name is GSE148505_raw_read_count_matrix.txt.gz

rawcounts <- read.csv("GSE148505_raw_count.csv")
dim(rawcounts) # [1] 23420    16
View(rawcounts)
## select distinct Gene.id row.
rawcounts.1 <- rawcounts %>%
    distinct(Geneid, .keep_all = T)
dim(rawcounts.1) # [1] 23418    16

## which raw names are duplicated
test <- rawcounts[duplicated(rawcounts$Geneid),]
View(test)

## Correct the duplicated raws from original data set,and renamed it.
GSE148505 <- read.csv("GSE148505_raw_count_corrected.csv")
dim(GSE148505) # [1] 23418    16
View(GSE148505)

### prepare a meta data file according the samples information for further analysis
metadata <- read.csv("GSE148505_metadata.csv", stringsAsFactors = T)
View(metadata)

##### select the wild type groups with and without DSS treated from the GSE148505 data set.
GSE148505wt <- GSE148505[ ,c(1:5, 11:13)]
View(GSE148505wt)

GSE148505wt <- GSE148505wt %>% mutate_if(is.numeric, round)

GSE148505$Geneid <- as.factor(GSE148505$Geneid)

## build a meta data contains the wild type groups with and without DSS treated 
metadata_wt <- metadata[c(1:4, 10:12),]
View(metadata_wt)

metadata_wt$dex=c("WTc","WTc","WTc","WTc","WTt","WTt","WTt")
metadata_wt$dex <- as.factor(metadata_wt$dex)

GSE148505wt$Geneid <- as.factor(GSE148505wt$Geneid)

### set up a reference group
metadata_wt$dex <- relevel(metadata_wt$dex, ref= "WTc")
library(DESeq2)
dds.wt <- DESeqDataSetFromMatrix(countData=GSE148505wt, 
                                colData=metadata_wt, 
                                design=~dex, 
                                tidy=TRUE)

dds.wt <- DESeq(dds.wt)
resultsNames(dds.wt) # [1] "Intercept"      "dex_WTt_vs_WTc"
#######################################################
res1.wt <- results(dds.wt, name = "dex_WTt_vs_WTc") 
res1.wt=as.data.frame(res1.wt)
View(res1.wt)
res1.wt =arrange(res1.wt,padj)
### save the differential expressed genes as the following file name
write.csv(res1.wt, file = "GSE148505deg_DESeq2_wt.csv")
####################################
library(dplyr)
library(tidyverse)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(data.table)
####################################
### to build a bar graph of differntial expressed selected kegg-pathway
GSE148505_deg <- read.csv(file = "GSE148505deg_DESeq2_wt.csv")
colnames(GSE148505_deg)
### fix the first colomn name 
colnames(GSE148505_deg) <- c("GeneName", "baseMean","log2FoldChange", "lfcSE","stat","pvalue", "padj")
View(GSE148505_deg) 
dim(GSE148505_deg) # [1] 23418     7

### clean the data-set 
GSE148505_deg <- na.omit(GSE148505_deg)
dim(GSE148505_deg) # [1] 18064     7
View(GSE148505_deg)
### save for further analysis 
write.csv(GSE148505_deg, file = "GSE148505_DESeq2_deg_gene_cleaned.csv")
##############  use this data to retrieve data to build the kegg-pathway bar graph #################
data <- read.csv("GSE148505_DESeq2_deg_gene_cleaned.csv")
View(data)
### clean the first column
data <- data[,-1]
row.names(data) <- data$GeneName
dim(data) #[1] 18064     7
data$GeneName <- toupper(data$GeneName)
View(data)
###########################################################
### select significantly expressed gene by the following parameters:
### baseMean > 10 & abs(log2FoldChange) > 1 & pvalue < 0.01 and filter data set
select.baseMean <- (data$baseMean >10)
select.log2FC <- abs(data$log2FoldChange) > 1
select.qval <- (data$pvalue < 0.01)

select.vec <- select.baseMean & select.log2FC & select.qval
table(select.vec)
# select.vec
# FALSE  TRUE 
# 16855  1209  
### 
degs.list = as.character(data$GeneName)[select.vec]
View(degs.list) # character 1209 gene names
#################################################################
library(org.Hs.eg.db) 
keytypes(org.Hs.eg.db)
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys =  degs.list,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

?enrichKEGG
erich.kegg.res <- enrichKEGG(gene = DEG.entrez_id,
                             organism = "hsa",
                             keyType = "kegg")


#Reading KEGG annotation online: "https://rest.kegg.jp/link/hsa/pathway"...
# Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/hsa"...

dotplot(erich.kegg.res)
zzh=erich.kegg.res@result
dir.create("enrich")

write.table(zzh,"./enrich/deg_GSE148505.enrich.kegg.res.txt",sep = "\t",col.names = NA)
kegg=read.table("./enrich/deg_GSE148505.enrich.kegg.res.txt",header = T,sep = "\t")
View(kegg)
################################
k = data.frame(kegg)
### build a graph 
library(ggplot2)
library(dplyr)
before <- as.numeric(sub("/\\d+$", "", k$GeneRatio))
after <- as.numeric(sub("^\\d+/", "", k$GeneRatio))
k$GeneRatio = before /after
font.size =8
### save the result
write.table(k,"./enrich/deg_GSE148505.enrich.kegg_generatio.txt",sep = "\t",
            col.names = NA)
###
k %>% 
  arrange(p.adjust) %>% 
  slice(1:30) %>% 
  ggplot(aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(range=c(3, 8))+
  labs(y=NULL) +
  ggtitle("")+
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))
####################################################################################
### Build up the Supplementary Figure 1e by selecting related Kegg signaling pathways.
### To demonstrate selected Kegg pathways using a bubble graph 
View(k)
k %>% 
  arrange(p.adjust) %>% 
  slice(1,3,6,7,10,16,22,27,25,30) %>% 
  ggplot(aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(range=c(3, 8))+
  labs(y=NULL) +
  ggtitle("")+
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black",
                                   size = font.size, vjust =1 ),
        axis.text.y = element_text(colour = "black",
                                   size = font.size, hjust =1 ),
        axis.title = element_text(margin=margin(10, 5, 0, 0),
                                  color = "black",size = font.size),
        axis.title.y = element_text(angle=90))
####################### END #######################################
