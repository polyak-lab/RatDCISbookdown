# Angiogenesis signature is from table S1 from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3743050/#mmc2, BC overexpressed angiogenesis signature.

gene_sig <- read.table("/Volumes/storage/ngs/analyzed_data/Carlos/190213_CGDA6432/analysis/angiogenesis_markers/BC_angiogenesis_overexpressed.txt", sep = '\t')
gene_sig <- gene_sig$V1

library(biomaRt)

Hmart = useMart("ensembl",dataset="hsapiens_gene_ensembl")
Rmart = useMart("ensembl",dataset="rnorvegicus_gene_ensembl")

orth = getBM(c("external_gene_name","rnorvegicus_homolog_associated_gene_name"), filters="external_gene_name", values = gene_sig, mart = Hmart)

library(DESeq2)
library(pheatmap)

countsTable <- read.csv('/Volumes/storage/ngs/analyzed_data/Carlos/190213_CGDA6432/combined-gene-level-counts.csv', header=TRUE, row.names=1, check.names = FALSE, comment.char = '', stringsAsFactors = FALSE)
drops <- c("CD45_6R_C", "DN_10L_B", "CD45_sample_36")
countsTable <- countsTable[,!names(countsTable) %in% drops]
countsTable <- countsTable[rowSums(countsTable) > 10,]
countsTable <- countsTable[,order(colnames(countsTable))]
infoTable <- read.table('/Volumes/storage/ngs/analyzed_data/Carlos/190213_CGDA6432/analysis/Differential_expression/sample-info.tsv', header=TRUE, sep='\t', row.names=1, check.names = FALSE, comment.char = '')
infoTable <- infoTable[!rownames(infoTable) %in% drops,]
infoTable <- infoTable[order(rownames(infoTable)),]

dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = infoTable, design = ~ 1)
dds <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds)

mat  <- assay(vsd)
mat <- mat[rownames(mat) %in% orth$rnorvegicus_homolog_associated_gene_name, ]
mat  <- mat - rowMeans(mat)

out <- pheatmap(mat, show_rownames = F, filename = "/Volumes/storage/ngs/analyzed_data/Carlos/190213_CGDA6432/analysis/angiogenesis_markers/angiogenesis_overexpressed_signature_clustering_all_samples.pdf")

#EP only
EP_samples <- rownames(infoTable[c(32:49),])
mat_EP <- assay(vsd)
mat_EP <- mat_EP[rownames(mat_EP) %in% orth$rnorvegicus_homolog_associated_gene_name, ]
mat_EP <- mat_EP[, colnames(mat_EP) %in% EP_samples]
mat_EP  <- mat_EP - rowMeans(mat_EP)
out <- pheatmap(mat_EP, show_rownames = F, filename = "/Volumes/storage/ngs/analyzed_data/Carlos/190213_CGDA6432/analysis/angiogenesis_markers/angiogenesis_overexpressed_signature_clustering_EP_samples.pdf")

DN_samples <- rownames(infoTable[c(17:31),])
mat_DN <- assay(vsd)
mat_DN <- mat_DN[rownames(mat_DN) %in% orth$rnorvegicus_homolog_associated_gene_name, ]
mat_DN <- mat_DN[, colnames(mat_DN) %in% DN_samples]
mat_DN  <- mat_DN - rowMeans(mat_DN)
out <- pheatmap(mat_DN, show_rownames = F, filename = "/Volumes/storage/ngs/analyzed_data/Carlos/190213_CGDA6432/analysis/angiogenesis_markers/angiogenesis_overexpressed_signature_clustering_DN_samples.pdf")

CD45_samples <- rownames(infoTable[c(1:16),])
mat_CD45 <- assay(vsd)
mat_CD45 <- mat_CD45[rownames(mat_CD45) %in% orth$rnorvegicus_homolog_associated_gene_name, ]
mat_CD45 <- mat_CD45[, colnames(mat_CD45) %in% CD45_samples]
mat_CD45  <- mat_CD45 - rowMeans(mat_CD45)
out <- pheatmap(mat_CD45, show_rownames = F, filename = "/Volumes/storage/ngs/analyzed_data/Carlos/190213_CGDA6432/analysis/angiogenesis_markers/angiogenesis_overexpressed_signature_clustering_CD45_samples.pdf")

