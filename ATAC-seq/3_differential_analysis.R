# Script for differential analysis / anotation / functional analysis
# https://tobiasrausch.com/courses/atac/atac-seq-data-analysis.html
# set up
rm(list = ls())
library(DESeq2)
library(ggplot2)
library(ChIPseeker)
library(rtracklayer)
library(clusterProfiler)
library(chipenrich)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# control panel
lfc_threshold <- 1
sign_threshold <- 0.01

# load count matrix
count_matrix <- read.csv('./count_matrix.csv')
rownames(count_matrix) <- count_matrix$peak
count_matrix$peak <- NULL
head(count_matrix)

# PCA
# your code here!

# Differential analysis
# your code here!
# creating the pheno data
pheno <- data.frame(cell_type = factor(c(rep("HL60", 3), rep("Mac3h", 3))))
pheno$cell_type <- relevel(pheno$cell_type, "HL60")
# creating the DESeq2 object
dds <- DESeqDataSetFromMatrix(count_matrix, pheno, ~ cell_type)
# rlog transformed data
rlog_transformed <- rlog(dds)
# pca
plotPCA(rlog_transformed, intgroup = "cell_type")
# Differential analysis
# your code here!
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
head(res)


#### Chipseeker: annotation ####

# loading the reference bed file
peaks <- import('./refererence.bed')

# creating the peak name
peaks$name <- paste0(seqnames(peaks), '_', start(peaks), '_', end(peaks))

# annotating
peakAnno <- annotatePeak(peaks, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnno

# plotting
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno)

# selecting the differentially expressed peaks
sign_idx <- which(abs(res$log2FoldChange) > lfc_threshold &
                    +                     res$padj <= sign_threshold)
sign_peaks <- peaks[peaks$name %in% rownames(res)[sign_idx]]

# annotation for significant peaks
signPeakAnno <- annotatePeak(sign_peaks, tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")

# plotting
plotAnnoBar(signPeakAnno)
plotDistToTSS(signPeakAnno)

#### enrichment analysis (cluster profiler) ####

# which are the significant genes?
sign_genes <- seq2gene(sign_peaks, tssRegion = c(-1000, 1000), 
                       flankDistance = 3000, TxDb = txdb)
sign_genes

# all genes
all_genes <- seq2gene(peaks, tssRegion = c(-1000, 1000), 
                       flankDistance = 3000, TxDb = txdb)

# GO ORA 
sign_genes <- unique(sign_genes[!grepl('ENS', sign_genes)])
all_genes <- unique(all_genes[!grepl('ENS', all_genes)])
enrich_res <- enrichGO(gene = sign_genes, 
                OrgDb = org.Hs.eg.db, 
                keyType = 'ENTREZID',
                universe = all_genes, 
                qvalueCutoff = 1, 
                minGSSize = 20, 
                maxGSSize = 200, 
                readable = TRUE)

# GO ORA results
head(enrich_res@result)

# dot plot
dotplot(enrich_res)

#### enriching according to genomic intervals ####

# writing the significant peaks for GREAT
to_export <- sign_peaks
to_export$name <- NULL
to_export$score <- NULL
export.bed(to_export, con = './sign_peaks.bed')

# test for broad peaks with ChipEnrich
gs_path = system.file('extdata','./vignette_genesets.txt', package='chipenrich')
gs <- read.table(gs_path, header = TRUE)
table(gs$gs_id)
broadenrich_res <- broadenrich(as.data.frame(sign_peaks), 
                               genome = 'hg38', genesets = gs_path, #'GOBP', 
                               #max_geneset_size = 200,
                               n_cores = 8)
# results chip enrich
head(broadenrich_res$results)
