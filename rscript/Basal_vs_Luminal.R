library(DESeq2)
library(data.table)
source('/Users/nickharper/Documents/random_scripts/R-files/RNA-Seq-functions.R')

#CD45 fraction:

countsTable <- read.csv('/Volumes/storage-1/ngs/analyzed_data/Carlos/190213_CGDA6432/combined-gene-level-counts.csv', header=TRUE, row.names=1, check.names = FALSE, comment.char = '', stringsAsFactors = FALSE)
drops <- c("CD45_6R_C", "DN_10L_B", "CD45_sample_36")
countsTable <- countsTable[,!names(countsTable) %in% drops]
countsTable <- countsTable[rowSums(countsTable) > 10,]
countsTable <- countsTable[,order(colnames(countsTable))]
infoTable <- read.table('/Volumes/storage-1/ngs/analyzed_data/Carlos/190213_CGDA6432/analysis/Differential_expression/sample-info.tsv', header=TRUE, sep='\t', row.names=1, check.names = FALSE, comment.char = '')
infoTable <- infoTable[!rownames(infoTable) %in% drops,]
infoTable <- infoTable[order(rownames(infoTable)),]

infoTable <- infoTable[c(1:16),]
infoTable$type <- c("Basal", "Luminal", "Basal", "Luminal", "Nothing", "Luminal", "Luminal", "Basal", "Basal", "Nothing", "Basal", "Luminal", "Nothing", "Basal", "Basal", "Basal")
infoTable <- infoTable[!infoTable$type == "Nothing",]
countsTable <- countsTable[,names(countsTable) %in% rownames(infoTable)]
countsTable <- countsTable[,order(colnames(countsTable))]
infoTable <- infoTable[order(rownames(infoTable)),]

dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = infoTable, design = ~ type)
dds <- estimateSizeFactors(dds)

padjFilter <- 0.05
absFCMin <- 2

dds <- DESeq(dds)

setwd("/Volumes/storage-1/ngs/analyzed_data/Carlos/190213_CGDA6432/analysis/basal_vs_luminal/")
CD45_Basal_vs_Luminal.DF <- DESeq2WriteDiff(dds, 'type', 'Basal', 'Luminal', outputFolder = './', outputFilePrefix = 'CD45_Basal_vs_Luminal', absLog2FCMin = log2(absFCMin), padjFilter = padjFilter)

rm(list=ls())

#DN fraction:
source('/Users/nickharper/Documents/random_scripts/R-files/RNA-Seq-functions.R')

countsTable <- read.csv('/Volumes/storage-1/ngs/analyzed_data/Carlos/190213_CGDA6432/combined-gene-level-counts.csv', header=TRUE, row.names=1, check.names = FALSE, comment.char = '', stringsAsFactors = FALSE)
drops <- c("CD45_6R_C", "DN_10L_B", "CD45_sample_36")
countsTable <- countsTable[,!names(countsTable) %in% drops]
countsTable <- countsTable[rowSums(countsTable) > 10,]
countsTable <- countsTable[,order(colnames(countsTable))]
infoTable <- read.table('/Volumes/storage-1/ngs/analyzed_data/Carlos/190213_CGDA6432/analysis/Differential_expression/sample-info.tsv', header=TRUE, sep='\t', row.names=1, check.names = FALSE, comment.char = '')
infoTable <- infoTable[!rownames(infoTable) %in% drops,]
infoTable <- infoTable[order(rownames(infoTable)),]

infoTable <- infoTable[c(17:31),]
infoTable$type <- c("Luminal", "Luminal", "Basal", "Luminal", "Luminal", "luminal", "Basal", "Luminal", "Luminal", "Basal", "Nothing", "Basal", "Luminal", "Basal", "Basal")
infoTable <- infoTable[!infoTable$type == "Nothing",]
countsTable <- countsTable[,names(countsTable) %in% rownames(infoTable)]
countsTable <- countsTable[,order(colnames(countsTable))]
infoTable <- infoTable[order(rownames(infoTable)),]

dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = infoTable, design = ~ type)
dds <- estimateSizeFactors(dds)

padjFilter <- 0.05
absFCMin <- 2

dds <- DESeq(dds)

DN_Basal_vs_Luminal.DF <- DESeq2WriteDiff(dds, 'type', 'Basal', 'Luminal', outputFolder = './', outputFilePrefix = 'DN_Basal_vs_Luminal', absLog2FCMin = log2(absFCMin), padjFilter = padjFilter)

rm(list=ls())

#EpCAM fraction:
source('/Users/nickharper/Documents/random_scripts/R-files/RNA-Seq-functions.R')

countsTable <- read.csv('/Volumes/storage-1/ngs/analyzed_data/Carlos/190213_CGDA6432/combined-gene-level-counts.csv', header=TRUE, row.names=1, check.names = FALSE, comment.char = '', stringsAsFactors = FALSE)
drops <- c("CD45_6R_C", "DN_10L_B", "CD45_sample_36")
countsTable <- countsTable[,!names(countsTable) %in% drops]
countsTable <- countsTable[rowSums(countsTable) > 10,]
countsTable <- countsTable[,order(colnames(countsTable))]
infoTable <- read.table('/Volumes/storage-1/ngs/analyzed_data/Carlos/190213_CGDA6432/analysis/Differential_expression/sample-info.tsv', header=TRUE, sep='\t', row.names=1, check.names = FALSE, comment.char = '')
infoTable <- infoTable[!rownames(infoTable) %in% drops,]
infoTable <- infoTable[order(rownames(infoTable)),]

infoTable <- infoTable[c(32:nrow(infoTable)),]
infoTable$type <- c("Basal", "Luminal", "Basal", "Luminal", "Basal", "Luminal", "Luminal", "Luminal", "Basal", "Luminal", "Luminal", "Basal", "Basal", "Luminal", "Basal", "Basal", "Basal", "Basal")
countsTable <- countsTable[,names(countsTable) %in% rownames(infoTable)]
countsTable <- countsTable[,order(colnames(countsTable))]
infoTable <- infoTable[order(rownames(infoTable)),]

dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = infoTable, design = ~ type)
dds <- estimateSizeFactors(dds)

padjFilter <- 0.05
absFCMin <- 2

dds <- DESeq(dds)

EpCAM_Basal_vs_Luminal.DF <- DESeq2WriteDiff(dds, 'type', 'Basal', 'Luminal', outputFolder = './', outputFilePrefix = 'EpCAM_Basal_vs_Luminal', absLog2FCMin = log2(absFCMin), padjFilter = padjFilter)

