###########
# Install #
###########
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("tximportData", "tximport", "ensembldb", "DESeq2", 
                       "IHW", "EnsDb.Hsapiens.v86", "ReportingTools", "apeglm", "clusterProfiler"))
BiocManager::install(c("clusterProfiler"))

library(tximportData)
library(tximport)
library(readr)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(dplyr)
library(forcats)
library(stringr)
library(gridExtra)
library(DESeq2)
library(genefilter)
library(RColorBrewer)
library(gplots)
library(clusterProfiler)
library(tibble)

###############
# IMPORT DATA #
###############
# You will need to change this to your current folder
setwd("~/honors")
dir <- getwd()
list.files(dir)

# Tx 2 Gene
edb <- EnsDb.Hsapiens.v86
# tx2gene <- ensembldb::transcripts(edb, return.type="DataFrame")

k <- keys(edb, keytype = "GENEID")
df <- select(edb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]  # tx ID, then gene ID

# Runs
samples <- read.table(file.path(dir, "acc.txt"), header = FALSE)

# Combine TSVs
files <- file.path(dir, samples$V1, "abundance.tsv")
names(files) <- paste0("sample", 1:54)

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
head(txi$counts)
length(txi$counts)

##############
# CONDITIONS #
##############

CONDITIONS = c("primary colorectal cancer", 
               "normal-looking surrounding colonic epithelium", 
               "metastatic colorectal cancer to the liver")

##########
# DeSeq2 #
##########

sampleTable <- data.frame(condition = factor(rep(CONDITIONS, each = 18)))
rownames(sampleTable) <- colnames(txi$counts)
ddsData <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
ddsData$condition <- relevel(ddsData$condition, ref = "normal-looking surrounding colonic epithelium")

dds <- DESeq(ddsData)
resultsNames(dds) # lists the coefficients

# Normal vs Primary
primary_vs_normal <- lfcShrink(dds,  type="apeglm", coef="condition_primary.colorectal.cancer_vs_normal.looking.surrounding.colonic.epithelium" )
primary_vs_normal$group <- "upregulated"
primary_vs_normal$group[primary_vs_normal$log2FoldChange < 0] <- "downregulated"
summary(primary_vs_normal)

up11 = length(which(primary_vs_normal$log2FoldChange > 1))
up12 = length(which(primary_vs_normal$log2FoldChange < -1))
(up11 / 35193) * 100
(up12 / 35193) * 100

# Normal vs Metastatic
metastatic_vs_normal <- lfcShrink(dds, type="apeglm", coef="condition_metastatic.colorectal.cancer.to.the.liver_vs_normal.looking.surrounding.colonic.epithelium")
metastatic_vs_normal$group <- "upregulated"
metastatic_vs_normal$group[metastatic_vs_normal$log2FoldChange < 0] <- "downregulated"
length(metastatic_vs_normal$group[metastatic_vs_normal$log2FoldChange > 0])
summary(metastatic_vs_normal)

up21 = length(which(metastatic_vs_normal$log2FoldChange > 1))
up22 = length(which(metastatic_vs_normal$log2FoldChange < -1))
(up21 / 35193) * 100
(up22 / 35193) * 100

#### FO
######R TESTING
#######
sampleTable2 <- data.frame(condition = factor(rep(CONDITIONS, each = 18)))
rownames(sampleTable2) <- colnames(txi$counts)
ddsData2 <- DESeqDataSetFromTximport(txi, sampleTable2, ~condition)
ddsData2$condition <- relevel(ddsData2$condition, ref = "metastatic colorectal cancer to the liver")

dds2 <- DESeq(ddsData2)
resultsNames(dds2) # lists the coefficients

# Primary vs Metastatic
primary_vs_metastatic <- lfcShrink(dds2, type="apeglm", coef="condition_primary.colorectal.cancer_vs_metastatic.colorectal.cancer.to.the.liver")
primary_vs_metastatic$group <- "upregulated"
primary_vs_metastatic$group[primary_vs_metastatic$log2FoldChange < 0] <- "downregulated"
summary(primary_vs_metastatic)

up31 = length(which(primary_vs_metastatic$log2FoldChange > 1))
up32 = length(which(primary_vs_metastatic$log2FoldChange < -1))
(up31 / 35193) * 100
(up32 / 35193) * 100


##################################
# PCA Analysis and Gene clustering
###################################

# PCA Plot
rld <- vst(dds)
plotPCA(rld, intgroup = c("condition"))

# Take 25 genes from each coefficient
topGenes = function (data, number) {
  resOrdered <- data[order(-abs(data$log2FoldChange)),]
  resOrdered = na.omit(resOrdered)
  resOrdered$ID <- mapIds(org.Hs.eg.db, rownames(resOrdered), "SYMBOL", "ENSEMBL")
  topVarGenes = na.omit(topVarGenes)
  topVarGenes <- head(topVarGenes, number)
  topVarGenes
}

topVarGenes = topGenes(metastatic_vs_normal, 25)
topVarGenes2 = topGenes(primary_vs_normal, 25)
topVarGenes3 = topGenes(primary_vs_metastatic, 25)

# Gene clustering
colors <- colorRampPalette( rev(brewer.pal(9, "PuOr")) )(255)
sidecols <- c("grey", "dodgerblue")[ rld$condition ]
mat <- assay(rld)[ c(rownames(topVarGenes), rownames(topVarGenes2), rownames(topVarGenes3)), ]
mat <- mat - rowMeans(mat)
colnames(mat) <- c(rep("primary CRC", 18), rep("normal colonic epithelium", 18), rep("metastatic liver CRC", 18))
heatmap.2(mat, 
          trace="none", 
          col=colors, 
          ColSideColors=sidecols,
          labRow=TRUE,
          mar=c(10,2), 
          scale="row")


#################
# GO Enrichment #
#################

generateGeneList = function (data) {
  d = data %>% 
    data.frame %>% 
    tibble::rownames_to_column("ID") %>% 
    dplyr::rename(FC=log2FoldChange) %>%
    dplyr::select(ID, FC)
  
  geneList = d[,2]
  names(geneList) = as.character(d[,1])
  geneList = sort(geneList, decreasing = TRUE)
  
  # LOG2 FOLD HAS TO BE GREATER THAN 1
  geneList <- names(geneList)[abs(geneList) > 1]
  
  return(geneList)
}

primary_vs_normal_geneList = generateGeneList(primary_vs_normal)
primary_vs_normal_go = enrichGO(gene          = primary_vs_normal_geneList,
                                OrgDb         = org.Hs.eg.db,
                                ont           = "ALL",
                                keyType       = "ENSEMBL",
                                pAdjustMethod = "holm",
                                pvalueCutoff  = 0.05,
                                readable      = TRUE)

metastatic_vs_normal_geneList = generateGeneList(metastatic_vs_normal)
metastatic_vs_normal_go = enrichGO(gene           = metastatic_vs_normal_geneList,
                                   OrgDb         = org.Hs.eg.db,
                                   ont           = "ALL",
                                   pAdjustMethod = "holm",
                                   keyType       = "ENSEMBL",
                                   pvalueCutoff  = 0.05,
                                   readable      = TRUE)

primary_vs_metastatic_geneList = generateGeneList(primary_vs_metastatic)
primary_vs_metastatic_go = enrichGO(gene           = primary_vs_metastatic_geneList,
                                    OrgDb         = org.Hs.eg.db,
                                    ont           = "ALL",
                                    pAdjustMethod = "holm",
                                    keyType       = "ENSEMBL",
                                    pvalueCutoff  = 0.05,
                                    readable      = TRUE)

primary_vs_normal_go
metastatic_vs_normal_go
primary_vs_metastatic_go

###############
## Visualize ##
###############

primary_vs_normal_go@result = primary_vs_normal_go@result %>% dplyr::arrange(desc(Count))
metastatic_vs_normal_go@result = metastatic_vs_normal_go@result %>% dplyr::arrange(desc(Count))
primary_vs_metastatic_go@result = primary_vs_metastatic_go@result %>% dplyr::arrange(desc(Count))

generateBarPlot = function (data) {
  barplot(data, showCategory=10, font.size=10)
}

generateDotPlot = function (data, title) {
  plot = dotplot(data, x=~ONTOLOGY, font.size=11, showCategory=10) +  
    scale_y_discrete(labels = function(y) str_wrap(y, width = 20)) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  return(plot)
}

barplot(primary_vs_normal_go)
barplot(metastatic_vs_normal_go)
barplot(primary_vs_metastatic_go)

plot1 = generateDotPlot(primary_vs_normal_go, "Primary vs Normal")
plot2 = generateDotPlot(metastatic_vs_normal_go, "Metastatic vs Normal")
plot3 = generateDotPlot(primary_vs_metastatic_go, "Metastatic vs Primary")
grid.arrange(plot1, plot2, plot3, nrow = 1)

########
# KEGG #
########

generateKegg = function (data) {
  result  <- mapIds(org.Hs.eg.db, data, "ENTREZID", "ENSEMBL") %>% 
    unname %>%
    enrichKEGG(gene          = .,
               organism      = 'hsa',
               pAdjustMethod = "holm",
               pvalueCutoff  = 0.05)
  return(result)
}

kk1 <- generateKegg(primary_vs_normal_geneList) %>% data.frame
kk2 <- generateKegg(metastatic_vs_normal_geneList) %>% data.frame
kk3 <- generateKegg(primary_vs_metastatic_geneList) %>% data.frame

rownames(kk1) <- NULL
rownames(kk2) <- NULL
rownames(kk3) <- NULL

head(kk1 %>% dplyr::select(ID, Description, GeneRatio, p.adjust))
head(kk2 %>% dplyr::select(ID, Description, GeneRatio, p.adjust))
head(kk3 %>% dplyr::select(ID, Description, GeneRatio, p.adjust))