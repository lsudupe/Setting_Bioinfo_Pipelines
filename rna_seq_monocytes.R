#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214282
#https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# libraries
library("NOISeq") #https://rstudio-pubs-static.s3.amazonaws.com/525119_64c1fe6e1a514b89a1ef26d23bf4aae3.html
library("tidyverse")
library("GEOquery")

##############
## PART 1: LOAD THE DATA
##############

##################
#   FROM GEO samples
##################
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE214282&format=file&file=GSE214282%5FrawCounts%2Etxt%2Egz"
path <- paste(urld, "acc=GSE214282", "file=GSE214282_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
GSE214282_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
tbl <- GSE214282_count

# pre-filter low count genes
# keep genes with at least 2 counts > 10
dim(tbl)
keep <- rowSums( tbl >= 5 ) >= 2
tbl <- tbl[keep, ]
dim(tbl)

# log transform raw counts
# instead of raw counts can display vst(as.matrix(tbl)) i.e. variance stabilized counts
dat <- log10(tbl + 1)

# box-and-whisker plot
dev.new(width=3+ncol(tbl)/6, height=5)
pdf("./results/RNAseq/GSE214282/GSE214282_boxplot_logcounts.pdf")
par(mar=c(7,4,2,1))
print(boxplot(dat, boxwex=0.7, notch=T, main="GSE214282", ylab="lg(cnt + 1)", outline=F, las=2))
dev.off()

# read metadata
gds <- getGEO("GSE214282")
Meta_GSE214282 <- pData(gds$GSE214282_series_matrix.txt.gz@phenoData)
Meta_GSE214282 <- Meta_GSE214282[,c("title","source_name_ch1","characteristics_ch1",
                                    "characteristics_ch1.1","description","cell type:ch1","disease state:ch1")]

# biological information: GC, gene length, chromosome
# Write the names
write.table(rownames(GSE214282_count),"./data/GSE214282/GSE214282_gene_names.entrez.txt",
            col.names = FALSE,row.names = FALSE,quote=F)

# Import the information
annotgene <- read.csv("./data/GSE214282/mart_export.txt",sep="\t",header = T)

# How many genes do I get annotated?
sum(rownames(GSE214282_count) %in% annotgene$Gene.stable.ID)

# Filter the information
annotgene <- annotgene[annotgene$Chromosome.scaffold.name %in% c(as.character(1:22) ,"X","Y"),]
sum(rownames(GSE214282_count) %in% annotgene$Gene.stable.ID)

## Annotation... solving some issues...
rownames(annotgene) <- annotgene$Gene.stable.ID
annotgene[annotgene$Gene.stable.ID=="ENSG00000012817",]

annotgene_filt <- annotgene[!duplicated(annotgene$Gene.stable.ID),]
sum(rownames(GSE214282_count) %in% annotgene$Gene.stable.ID)
sum(annotgene_filt$Gene.stable.ID %in% rownames(GSE214282_count))
annotgene_filt[annotgene_filt$Gene.stable.ID=="ENSG00000012817",]

## Overlap between annotation and gnes
rownames(annotgene_filt) <- as.character(annotgene_filt$Gene.stable.ID)
sum(as.character(rownames(annotgene_filt)) %in% rownames(GSE214282_count))


##  Work with the annotated genes!
GSE214282_count_filt <- GSE214282_count[rownames(GSE214282_count) %in% rownames(annotgene_filt),] #dim 58402     8
GSE214282_count_exc <-GSE214282_count[!(rownames(GSE214282_count) %in% rownames(annotgene_filt)),] #dim 482   8
annotgene_ord <- annotgene_filt[rownames(GSE214282_count_filt ),]

sum(rownames(annotgene_ord)==rownames(GSE214282_count_filt))

## prepare data 
Meta_GSE214282$title <- gsub("-", ".", Meta_GSE214282$title)
rownames(Meta_GSE214282) <- Meta_GSE214282$title
Factors_GSE214282 <- data.frame(Meta_GSE214282[colnames(GSE214282_count_filt),c("disease state:ch1")])
colnames(Factors_GSE214282)[1]<- "Group"

GSE214282_count_filt
annotgene_ord
Factors_GSE214282

##############
## PART 2: QC
##############

######NOISeq
library(NOISeq)

## prepare data for NOISeq
data_NOISEQ <- readData(data = GSE214282_count_filt,
                        length=abs(annotgene_ord$Gene.end..bp.-annotgene_ord$Gene.start..bp.),
                        gc=annotgene_ord$Gene...GC.content,
                        biotype= annotgene_ord$type ,
                        chromosome = annotgene_ord[,c("Chromosome.scaffold.name","Gene.start..bp.","Gene.end..bp.")],
                        factors = Factors_GSE214282)
# problems?


myexplodata <- dat(data_NOISEQ, type = "biotype")
explo.plot(myexplodata, plottype = "persample")
mynicedata <- dat2save(myexplodata)
mybiodetection <- dat(data_NOISEQ, k = 0, type = "biodetection", factor = NULL)


lengthuse <- abs(annotgene_ord$Gene.end..bp.-annotgene_ord$Gene.start..bp.)
names(lengthuse) <- rownames(annotgene_ord)
gc <- annotgene_ord$Gene...GC.content
names(gc) <- rownames(annotgene_ord)
biotype <-annotgene_ord$Transcript.type
names(biotype) <- rownames(annotgene_ord)

chromosome <- annotgene_ord[,c("Chromosome.scaffold.name","Gene.start..bp.","Gene.end..bp.")]

## now yes
data_NOISEQ <- readData(data = GSE214282_count_filt,
                        length=lengthuse,
                        gc=gc,
                        biotype= biotype ,
                        chromosome = annotgene_ord[,c("Chromosome.scaffold.name","Gene.start..bp.","Gene.end..bp.")],
                        factors = Factors_GSE214282)

myexplodata <- dat(data_NOISEQ, type = "biodetection")

## explore our data
pdf("./results/RNAseq/GSE214282/GSE214282_type.pdf")
print(explo.plot(myexplodata, plottype = "persample"))
dev.off()

pdf("./results/RNAseq/GSE214282/GSE214282_biotype_detection.pdf", width = 12, height = 8)
par(mfrow = c(1, 2))
explo.plot(myexplodata, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")
dev.off()

mycountsbio = dat(data_NOISEQ, factor = NULL, type = "countsbio")
pdf("./results/RNAseq/GSE214282/GSE214282_each_expressionvalues.pdf", width = 10, height = 8)
par(mfrow = c(1, 2))
print(explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot"))
dev.off()

mysaturation = dat(data_NOISEQ, k = 0, ndepth = 7, type = "saturation")
pdf("./results/RNAseq/GSE214282/GSE214282_sequencingdepth.pdf", width = 10, height = 8)
par(mfrow = c(2, 1))
explo.plot(mysaturation, toplot = 1, samples = 1:8, yleftlim = NULL, yrightlim = NULL)
explo.plot(mysaturation, toplot = "protein_coding", samples = 1:8)
dev.off()

pdf("./results/RNAseq/GSE214282/GSE214282_cpm.pdf", width = 10, height = 8)
explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")
dev.off()

mylengthbias = dat(data_NOISEQ, factor = "Group", type = "lengthbias")
pdf("./results/RNAseq/GSE214282/GSE214282_meanExpression_lengthBins.pdf", width = 10, height = 8)
explo.plot(mylengthbias, samples = NULL, toplot = "global")
dev.off()

myGCbias = dat(data_NOISEQ, factor = "Group", type = "GCbias")
pdf("./results/RNAseq/GSE214282/GSE214282_gc_content.pdf", width = 10, height = 8)
explo.plot(myGCbias, samples = NULL, toplot = "global")
dev.off()

mycd = dat(data_NOISEQ, type = "cd", norm = FALSE, refColumn = 1)
pdf("./results/RNAseq/GSE214282/GSE214282_cd_density.pdf", width = 10, height = 8)
explo.plot(mycd,samples = 1:8)
dev.off()

myPCA = dat(data_NOISEQ, type = "PCA")
pdf("./results/RNAseq/GSE214282/GSE214282_pca.pdf", width = 10, height = 8)
explo.plot(myPCA, factor = "Group")
dev.off()

pdf("./results/RNAseq/GSE214282/GSE214282_qa.pdf", width = 10, height = 8)
QCreport(data_NOISEQ, samples = NULL, factor = "Group", norm = FALSE)
dev.off()

## save 
save(data_NOISEQ,GSE214282_count_filt,annotgene_ord,file="./data/GSE214282/GSE214282_step1.Rda")

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
load("./data/GSE214282/GSE214282_step1.Rda")

############
# STEP 3.1.1: SET THE CLASS
############

GSE214282_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE214282_count_filt,
                                           colData = pData(data_NOISEQ),
                                           design = ~ Group)
# Warning
pDataUSE <- pData(data_NOISEQ)
pDataUSE[pDataUSE=="Control"] <- "control"
pDataUSE[pDataUSE=="CFS"] <- "case"
pDataUSE[,1] <- as.factor(pDataUSE[,1])
pDataUSE$Group <- relevel(pDataUSE$Group, ref = "control")

GSE214282_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE214282_count_filt,
                                           colData = pDataUSE,
                                           design = ~ -1 + Group)
resultsNames(GSE214282_DESeq2)
GSE214282_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE214282_count_filt,
                                           colData = pDataUSE,
                                           design = ~ Group)

GSE214282_DESeq2

############
# STEP 3.1.2: WITH WHICH GENES TO WORK?
############

## Do we use all the genes?
## How do we select which ones?

smallestGroupSize <- 1
keep <- rowSums(counts(GSE214282_DESeq2) >= 2) >= smallestGroupSize
GSE214282_DESeq2_F <- GSE214282_DESeq2[keep,]


############
# STEP 3.1.3: DIFFERENTIAL EXPRESSION?
############

GSE214282_DESeq2_F<- DESeq(GSE214282_DESeq2_F)
GSE214282_res <- results(GSE214282_DESeq2_F)
GSE214282_res
resultsNames(GSE214282_DESeq2_F)

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
pdf("./results/RNAseq/GSE214282/GSE214282_MA.pdf", width = 10, height = 8)
plotMA(GSE214282_res, ylim=c(-2,2))
dev.off()

lfcShrink(GSE214282_DESeq2_F,coef=c("Group_case_vs_control"))
res_lfcShrink <- lfcShrink(GSE214282_DESeq2_F,coef=c("Group_case_vs_control"))

pdf("./results/RNAseq/GSE214282/GSE214282_MA_lfc.pdf", width = 10, height = 8)
plotMA(res_lfcShrink, ylim=c(-2,2))
dev.off()

GSE214282_res
resOrdered <- GSE214282_res[order(GSE214282_res$pvalue),]
summary(GSE214282_res)

## STEP 3.1.4: Define questions
 
resultsNames(GSE214282_DESeq2_F)

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

resNorm <- lfcShrink(GSE214282_DESeq2_F, coef=2, type="normal")
resAsh <- lfcShrink(GSE214282_DESeq2_F, coef=2, type="ashr")

pdf("./results/RNAseq/GSE214282/GSE214282_MA_plots.pdf", width = 14, height = 5)
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res_lfcShrink, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
dev.off()

## STEP 3.1.4: QC??

# How do visualize?
vsd <- vst(GSE214282_DESeq2_F, blind=FALSE)
rld <- rlog(GSE214282_DESeq2_F, blind=FALSE)
head(assay(vsd), 3)

# this gives log2(n + 1)
ntd <- normTransform(GSE214282_DESeq2_F)
library("vsn")
pdf("./results/RNAseq/GSE214282/GSE214282_effects_variation_transf1.pdf", width = 14, height = 5)
meanSdPlot(assay(ntd))
dev.off()

pdf("./results/RNAseq/GSE214282/GSE214282_effects_variation_transf2.pdf", width = 14, height = 5)
meanSdPlot(assay(vsd))
dev.off()


pdf("./results/RNAseq/GSE214282/GSE214282_effects_variation_transf3.pdf", width = 14, height = 5)
meanSdPlot(assay(rld))
dev.off()

############
# STEP 4: BIOLOGICAL INTERPRETATION
############

## heatmap
library("pheatmap")
select <- order(rowMeans(counts(GSE214282_DESeq2_F,normalized=TRUE)),
                decreasing=TRUE)[1:20]
ntd_data <- assay(ntd)[select,]

df <- as.data.frame(colData(GSE214282_DESeq2_F)[,c("control","case")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
df <- as.data.frame(colData(GSE214282_DESeq2_F)[, "Group", drop = FALSE])

# Plot the heatmap
pdf("./results/RNAseq/GSE214282/GSE214282_heatmap.pdf", width = 10, height = 12)
pheatmap(ntd_data, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()
#change names genes
gene_names <- setNames(annotgene$Gene.name, annotgene$Gene.stable.ID)
rownames(GSE214282_res) <- gene_names[rownames(GSE214282_res)]
head(rownames(GSE214282_res))

library(EnhancedVolcano)
pdf("./results/RNAseq/GSE214282/GSE214282_volcano.pdf", width = 10, height = 8)
EnhancedVolcano(GSE214282_res,
                lab = rownames(GSE214282_res),
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

results_anno <- GSE214282_res

## cluster profiler
OrgDb <- org.Hs.eg.db 

geneList <- as.vector(results_anno$log2FoldChange)
names(geneList) <- rownames(results_anno)
gene <- as.vector(rownames(results_anno))

# First, convert gene symbols to Entrez IDs
entrez_ids <- AnnotationDbi::mapIds(OrgDb,
                                    keys = gene,
                                    column = "ENTREZID",
                                    keytype = "SYMBOL",
                                    multiVals = "first")


entrez_ids_vector <- unlist(valid_entrez_ids)
entrez_ids_vector <- as.character(entrez_ids_vector)
head(entrez_ids_vector)

# Group GO
ggo <- clusterProfiler::groupGO(gene     = entrez_ids_vector,
                                OrgDb    = OrgDb,
                                ont      = "BP",
                                level    = 3,
                                readable = TRUE)
head(summary(ggo)[,-5])

pdf("./results/RNAseq/GSE214282/GSE214282_clusterprofiler_barplot.pdf", width = 10, height = 8)
barplot(ego, showCategory=25)
dev.off()

# GO over-representation test
ego <- clusterProfiler::enrichGO(gene          = entrez_ids_vector,
                                 OrgDb         = OrgDb,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.03,
                                 qvalueCutoff  = 0.03, 
                                 readable      = TRUE)
head(summary(ego)[,-8])

pdf("./results/RNAseq/GSE214282/GSE214282_clusterprofiler_ego.pdf", width = 8, height = 10)
barplot(ego, showCategory=25)
dev.off()
