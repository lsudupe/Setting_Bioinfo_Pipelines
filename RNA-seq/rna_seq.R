

# libraries
library("NOISeq") #https://rstudio-pubs-static.s3.amazonaws.com/525119_64c1fe6e1a514b89a1ef26d23bf4aae3.html
library("tidyverse")
library("GEOquery")

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
GSE198256_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# read metadata
#meta <- read.table("./data/GSE198256/SraRunTable.txt" , sep=",", header = TRUE)
gds <- getGEO("GSE198256")
Meta_GSE198256 <- pData(gds$GSE198256_series_matrix.txt.gz@phenoData)
Meta_GSE198256 <- Meta_GSE198256[,c("title","source_name_ch1","characteristics_ch1",
                                    "characteristics_ch1.1","description","cell type:ch1","disease state:ch1")]

# extract rownames
genes <- rownames(tbl)
writeLines(genes, "./data/GSE198256/genes_names.txt")

# read biomart export
annotgene <- read.table("./data/GSE198256/mart_export_david.txt",sep = "\t", header = TRUE)
rownames(annotgene) <- annotgene$Entrezgene # some are repited that could be associated to different isophorms

# How many genes do I get annotated?
sum(rownames(GSE198256_count) %in% annotgene$Entrezgene)

# Filter the information
annotgene <- annotgene[annotgene$Chromosome %in% c(as.character(1:22) ,"X","Y"),]
sum(rownames(GSE198256_count) %in% annotgene$Entrezgene)

## Annotation... solving some issues...
rownames(annotgene) <- annotgene$Entrezgene
annotgene[annotgene$Entrezgene=="132989",]

annotgene_filt <- annotgene[!duplicated(annotgene$Entrezgene),]
sum(rownames(GSE198256_count) %in% annotgene$Entrezgene)
sum(annotgene_filt$Entrezgene %in% rownames(GSE198256_count))
annotgene_filt[annotgene_filt$Entrezgene=="132989",]

## Overlap between annotation and gnes
rownames(annotgene_filt) <- as.character(annotgene_filt$Entrezgene)
sum(as.character(rownames(annotgene_filt)) %in% rownames(GSE198256_count))

##  Work with the annotated genes!
GSE198256_count_filt <- GSE198256_count[rownames(GSE198256_count) %in% rownames(annotgene_filt),]
GSE198256_count_exc <-GSE198256_count[!(rownames(GSE198256_count) %in% rownames(annotgene_filt)),]
annotgene_ord <- annotgene_filt[rownames(GSE198256_count_filt ),]

sum(rownames(annotgene_ord)==rownames(GSE198256_count_filt))


######NOISeq
Factors_GSE198256 <- data.frame(Meta_GSE198256 [ colnames(GSE198256_count_filt),c("disease state:ch1")])
colnames(Factors_GSE198256)[1]<- "Group"

data_NOISEQ <- readData(data = GSE198256_count_filt,
                        length=abs(annotgene_ord$end-annotgene_ord$start),
                        gc=annotgene_ord$GC,
                        biotype= annotgene_ord$type ,
                        chromosome = annotgene_ord[,c("Chromosome","start","end")],
                        factors = Factors_GSE198256)
# problems?
myexplodata <- dat(data_NOISEQ, type = "biotype")
explo.plot(myexplodata, plottype = "persample")
mynicedata <- dat2save(myexplodata)
mybiodetection <- dat(data_NOISEQ, k = 0, type = "biodetection", factor = NULL)


lengthuse <- abs(annotgene_ord$end-annotgene_ord$start)
names(lengthuse) <- rownames(annotgene_ord)
gc <- annotgene_ord$GC
names(gc) <- rownames(annotgene_ord)
biotype <-annotgene_ord$type
names(biotype) <- rownames(annotgene_ord)

chromosome <- annotgene_ord[,c("Chromosome","start","end")]


data_NOISEQ <- readData(data = GSE198256_count_filt,
                        length=lengthuse,
                        gc=gc,
                        biotype= biotype ,
                        chromosome = annotgene_ord[,c("Chromosome","start","end")],
                        factors = Factors_GSE198256)

myexplodata <- dat(data_NOISEQ, type = "biodetection")
myexplodata <- dat(data_NOISEQ, type = "biodetection")

## explore our data
pdf("./results/RNAseq/GSE198256/GSE198256_type.pdf")
print(explo.plot(myexplodata, plottype = "persample"))
dev.off()

pdf("./results/RNAseq/GSE198256/GSE198256_biotype_detection.pdf", width = 12, height = 8)
par(mfrow = c(1, 2))
explo.plot(myexplodata, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")
dev.off()

mycountsbio = dat(data_NOISEQ, factor = NULL, type = "countsbio")
pdf("./results/RNAseq/GSE198256/GSE198256_each_expressionvalues.pdf", width = 10, height = 8)
par(mfrow = c(1, 2))
print(explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot"))
dev.off()

mysaturation = dat(data_NOISEQ, k = 0, ndepth = 7, type = "saturation")
pdf("./results/RNAseq/GSE198256/GSE198256_sequencingdepth.pdf", width = 10, height = 8)
par(mfrow = c(2, 1))
explo.plot(mysaturation, toplot = 1, samples = 1:8, yleftlim = NULL, yrightlim = NULL)
explo.plot(mysaturation, toplot = "protein_coding", samples = 1:8)
dev.off()

pdf("./results/RNAseq/GSE198256/GSE198256_cpm.pdf", width = 10, height = 8)
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")
dev.off()

mylengthbias = dat(data_NOISEQ, factor = "Group", type = "lengthbias")
pdf("./results/RNAseq/GSE198256/GSE198256_meanExpression_lengthBins.pdf", width = 10, height = 8)
explo.plot(mylengthbias, samples = NULL, toplot = "global")
dev.off()

myGCbias = dat(data_NOISEQ, factor = "Group", type = "GCbias")
pdf("./results/RNAseq/GSE198256/GSE198256_gc_content.pdf", width = 10, height = 8)
explo.plot(myGCbias, samples = NULL, toplot = "global")
dev.off()

mycd = dat(data_NOISEQ, type = "cd", norm = FALSE, refColumn = 1)
pdf("./results/RNAseq/GSE198256/GSE198256_cd_density.pdf", width = 10, height = 8)
explo.plot(mycd,samples = 1:8)
dev.off()

myPCA = dat(data_NOISEQ, type = "PCA")
pdf("./results/RNAseq/GSE198256/GSE198256_pca.pdf", width = 10, height = 8)
explo.plot(myPCA, factor = "Group")
dev.off()

pdf("./results/RNAseq/GSE198256/GSE198256_qa.pdf", width = 10, height = 8)
QCreport(data_NOISEQ, samples = NULL, factor = "Group", norm = FALSE)
dev.off()

## save 
save(data_NOISEQ,GSE198256_count_filt,annotgene_ord,file="./data/GSE198256/GSE198256_step1.Rda")

############
## STEP 3: NORMALIZATION & DIFF EXPRESSION
############

myRPKM = rpkm(assayData(data_NOISEQ)$exprs, long = lengthuse, k = 0, lc = 1)
myUQUA = uqua(assayData(data_NOISEQ)$exprs, long = lengthuse, lc = 0.5, k = 0)
myTMM = tmm(assayData(data_NOISEQ)$exprs, long = 1000, lc = 0)

############
## STEP 3.1: DESEQ2
############
library(DESeq2)
load("./data/GSE198256/GSE198256_step1.Rda")

############
# STEP 3.1.1: SET THE CLASS
############

GSE198256_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt,
                                           colData = pData(data_NOISEQ),
                                           design = ~ Group)
# Warning
# Warning
pDataUSE <- pData(data_NOISEQ)
pDataUSE[pDataUSE=="Covid19: Acute infection"] <- "Covid19AI"
pDataUSE[pDataUSE=="Covid19: Recovery 3Mo"] <- "Covid193Mo"
pDataUSE[pDataUSE=="Covid19: Recovery 6Mo"] <- "Covid196Mo"
pDataUSE[,1] <- as.factor(pDataUSE[,1])

GSE198256_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt,
                                           colData = pDataUSE,
                                           design = ~ -1 + Group)
resultsNames(GSE198256_DESeq2)
GSE198256_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt,
                                           colData = pDataUSE,
                                           design = ~ Group)

GSE198256_DESeq2

############
# STEP 3.1.2: WITH WHICH GENES TO WORK?
############

## Do we use all the genes?
## How do we select which ones?

smallestGroupSize <- 6
keep <- rowSums(counts(GSE198256_DESeq2) >= 10) >= smallestGroupSize
GSE198256_DESeq2_F <- GSE198256_DESeq2[keep,]


############
# STEP 3.1.3: DIFFERENTIAL EXPRESSION?
############

GSE198256_DESeq2_F<- DESeq(GSE198256_DESeq2_F)
GSE198256_res <- results(GSE198256_DESeq2_F)
GSE198256_res
resultsNames(GSE198256_DESeq2_F)

############
# STEP 3.1.4: WE NEED TO UNDERSTAND MORE...
############

## Questions in my mind:
# How do I define the question?
# How the differential expression is done?
# How to interpret the results?
# Technical replicates?

## STEP 3.1.4: plot MA

#Interpretation?
pdf("./results/RNAseq/GSE198256/GSE198256_MA.pdf", width = 10, height = 8)
plotMA(GSE198256_res, ylim=c(-2,2))
dev.off()

#lfcShrink(GSE198256_DESeq2_F,coef=c("Group_case_vs_control"))
#res_lfcShrink <- lfcShrink(GSE198256_DESeq2_F,coef=c("Group_case_vs_control"))

#pdf("./results/RNAseq/GSE198256/GSE198256_MA_lfc.pdf", width = 10, height = 8)
#plotMA(res_lfcShrink, ylim=c(-2,2))
#dev.off()

#GSE198256_res
#resOrdered <- GSE198256_res[order(GSE198256_res$pvalue),]
#summary(GSE198256_res)


## STEP 3.1.4: Define questions
GSE198256_DESeq2_F<- DESeq(GSE198256_DESeq2_F)
#res <- results(GSE198256_DESeq2_F, contrast=c('factorName','numeratorLevel','denominatorLevel'))
res <- results(GSE198256_DESeq2_F, contrast=c("Group","Healthy","Covid19AI"))
res
resultsNames(GSE198256_DESeq2_F)

GSE198256_DESeq2_F <- DESeq(GSE198256_DESeq2_F, test="LRT", reduced=~1)
GSE198256_DESeq2_res_LRT <- results(GSE198256_DESeq2_F)
GSE198256_DESeq2_res_LRT
res <- results(GSE198256_DESeq2_res_LRT)

# Technical replicates?

# How to interpret the results?

plotCounts(GSE198256_DESeq2_F, gene="100287102", intgroup="Group")

## STEP 3.1.4: How differential expression is conducted...

# DESeq2 offers two kinds of hypothesis tests: 
#   the Wald test, 
#        where we use the estimated standard error of a log2 fold 
#        change to test if it is equal to zero, 
#   the likelihood ratio test (LRT). 
#        The LRT examines two models for the counts, a full model 
#        with a certain number of terms and a reduced model, in 
#        which some of the terms of the full model are removed. 
#        The test determines if the increased likelihood of the 
#        data using the extra terms in the full model is more 
#        than expected if those extra terms are truly zero.

resNorm <- lfcShrink(GSE198256_DESeq2_F, coef=2, type="normal")
resAsh <- lfcShrink(GSE198256_DESeq2_F, coef=2, type="ashr")

pdf("./results/RNAseq/GSE198256/GSE198256_MA_plots.pdf", width = 14, height = 5)
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res_lfcShrink, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
dev.off()

## STEP 3.1.4: QC??

# How do visualize?
vsd <- vst(GSE198256_DESeq2_F, blind=FALSE)
rld <- rlog(GSE198256_DESeq2_F, blind=FALSE)
head(assay(vsd), 3)

# this gives log2(n + 1)
ntd <- normTransform(GSE198256_DESeq2_F)
library("vsn")
pdf("./results/RNAseq/GSE198256/GSE198256_effects_variation_transf1.pdf", width = 14, height = 5)
meanSdPlot(assay(ntd))
dev.off()

pdf("./results/RNAseq/GSE198256/GSE198256_effects_variation_transf2.pdf", width = 14, height = 5)
meanSdPlot(assay(vsd))
dev.off()


pdf("./results/RNAseq/GSE198256/GSE198256_effects_variation_transf3.pdf", width = 14, height = 5)
meanSdPlot(assay(rld))
dev.off()

############
# STEP 4: BIOLOGICAL INTERPRETATION
############

## heatmap
library("pheatmap")
select <- order(rowMeans(counts(GSE198256_DESeq2_F,normalized=TRUE)),
                decreasing=TRUE)[1:20]
ntd_data <- assay(ntd)[select,]

df <- as.data.frame(colData(GSE198256_DESeq2_F)[, "Group", drop = FALSE])
pdf("./results/RNAseq/GSE198256/GSE198256_heatmap.pdf", width = 10, height = 12)
print(pheatmap(assay(ntd)[select,], 
         cluster_rows=FALSE, 
         show_rownames=TRUE,
         cluster_cols=TRUE, 
         annotation_col=df))
dev.off()

#change names genes
gene_names <- setNames(annotgene$Entrezgene, annotgene$Entrezgene)
rownames(GSE198256_res) <- gene_names[rownames(GSE198256_res)]
head(rownames(GSE198256_res))

library(EnhancedVolcano)
pdf("./results/RNAseq/GSE198256/GSE198256_volcano.pdf", width = 10, height = 8)
EnhancedVolcano(GSE198256_res,
                lab = rownames(GSE198256_res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-10,10))
dev.off()

## Gene Set Enrichment Analysis.
#    a. ORA.
#    b. GSEA
# How do we bring all the information together at once?

BiocManager::install("topGO")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(geneAnnotations)
library(exprAnalysis)

results_anno <- GSE198256_res

## cluster profiler
OrgDb <- org.Hs.eg.db 

#geneList <- as.vector(results_anno$log2FoldChange)
#names(geneList) <- rownames(results_anno)
#gene <- as.vector(rownames(results_anno))

# First, convert gene symbols to Entrez IDs
#entrez_ids <- AnnotationDbi::mapIds(OrgDb,
#                                    keys = gene,
#                                    column = "ENTREZID",
#                                    keytype = "SYMBOL",
#                                    multiVals = "first")


#entrez_ids_vector <- unlist(valid_entrez_ids)
#entrez_ids_vector <- as.character(entrez_ids_vector)
#head(entrez_ids_vector)

# Group GO
ggo <- clusterProfiler::groupGO(gene     = gene,
                                OrgDb    = OrgDb,
                                ont      = "BP",
                                level    = 3,
                                readable = TRUE)
head(summary(ggo)[,-5])

pdf("./results/RNAseq/GSE198256/GSE198256_clusterprofiler_barplot.pdf", width = 10, height = 8)
barplot(ego, showCategory=25)
dev.off()

# GO over-representation test
ego <- clusterProfiler::enrichGO(gene          = gene,
                                 OrgDb         = OrgDb,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.03,
                                 qvalueCutoff  = 0.03, 
                                 readable      = TRUE)
head(summary(ego)[,-8])

pdf("./results/RNAseq/GSE214282/GSE198256_clusterprofiler_ego.pdf", width = 8, height = 10)
barplot(ego, showCategory=25)
dev.off()
