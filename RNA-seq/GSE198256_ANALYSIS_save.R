setwd("~/Dropbox/WORKING_KAUST/LECTURING/B3XX_SettinngPipelines/DATA_Sets/GSE198256")

##############
## PART 1: LOAD THE DATA
##############

##################
#   FROM GEO samples
##################

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# pre-filter low count genes
# keep genes with at least 2 counts > 10
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]

# log transform raw counts
# instead of raw counts can display vst(as.matrix(tbl)) i.e. variance stabilized counts
dat <- log10(tbl + 1)

# box-and-whisker plot
dev.new(width=3+ncol(tbl)/6, height=5)
pdf("./results/RNAseq/GSE198256/GSE198256_boxplot_logcounts.pdf")
par(mar=c(7,4,2,1))
boxplot(dat, boxwex=0.7, notch=T, main="GSE198256", ylab="lg(cnt + 1)", outline=F, las=2)
dev.off()

# UMAP plot (dimensionality reduction)
library(umap)
dat <- dat[!duplicated(dat), ] # first remove duplicates
ump <- umap(t(dat), n_neighbors = 14, random_state = 123)
plot(ump$layout, main="GSE198256 UMAP plot, nbrs =14", xlab="", ylab="", pch=20, cex=1.5)
library(car)
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

##################
#   FROM DOWNLOAD 
##################

data_GSE198256 <- read.table("GSE198256_Monocytes_Covid19_Normalized_counts.tsv",sep="\t",
                             row.names = 1,header = T)

##################
#   COMPARING
##################

# Gene names
# Data processed.
# XLS vs TSV

##################
#   BETTER COUNT
##################

urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
GSE198256_count <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)


##############
## PART 2: QC
##############

######
## WHAT TO CHECK
######

######
## NOISEQ
######
###### https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("NOISeq")

## WE NEED THE METADATA

#BiocManager::install("GEOquery")
library(GEOquery)
## https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html
gds <- getGEO("GSE198256")
Meta_GSE198256 <- pData(gds$GSE198256_series_matrix.txt.gz@phenoData)
Meta_GSE198256 <- Meta_GSE198256[,c("title","source_name_ch1","characteristics_ch1","characteristics_ch1.1","description","cell type:ch1","disease state:ch1")]

Factors_GSE198256 <- Meta_GSE198256[,c("disease state:ch1")]
GSE198256_count

## WE NEED BIOLOGICAL INFORMATION: GC, GENE LENGTH, CHROMOSOME,...

# Write the names
write.table(rownames(GSE198256_count),"gene_names.entrez.txt",
            col.names = FALSE,row.names = FALSE,quote=F)

# Additional Biological information.
# https://www.ensembl.org/biomart/martview/7f2a95d66853c3b8aea7639401e47aba

# Import the information
annotgene <- read.csv("mart_export.txt",sep="\t",header = T)
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
  
# READY
GSE198256_count_filt
annotgene_ord
Factors_GSE198256
  
######
library(NOISeq)
#BiocManager::install("NOISeq",force = TRUE)

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

pdf("./results/RNAseq/GSE198256/GSE198256_type.pdf")
print(explo.plot(myexplodata, plottype = "persample"))
dev.off()


par(mfrow = c(1, 2))
explo.plot(myexplodata, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")


mycountsbio = dat(data_NOISEQ, factor = NULL, type = "countsbio")
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")

mysaturation = dat(data_NOISEQ, k = 0, ndepth = 7, type = "saturation")
explo.plot(mysaturation, toplot = 1, samples = 1:2, yleftlim = NULL, yrightlim = NULL)
explo.plot(mysaturation, toplot = "protein_coding", samples = 1:4)

explo.plot(mycountsbio, toplot = "protein_coding", samples = NULL, plottype = "boxplot")

explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")

mylengthbias = dat(data_NOISEQ, factor = "Group", type = "lengthbias")
explo.plot(mylengthbias, samples = NULL, toplot = "global")

myGCbias = dat(data_NOISEQ, factor = "Group", type = "GCbias")
explo.plot(myGCbias, samples = NULL, toplot = "global")

mycd = dat(data_NOISEQ, type = "cd", norm = FALSE, refColumn = 1)
explo.plot(mycd,samples = 1:12)

myPCA = dat(data_NOISEQ, type = "PCA")
explo.plot(myPCA, factor = "Group")

QCreport(data_NOISEQ, samples = NULL, factor = "Group", norm = FALSE)

########
##### SAVE JUST IN CASE
########

save(data_NOISEQ,GSE198256_count_filt,annotgene_ord,file="GSE198256_step1.Rda")

############
## STEP 3: NORMALIZATION & DIFF EXPRESSION
############

myRPKM = rpkm(assayData(data_NOISEQ)$exprs, long = lengthuse, k = 0, lc = 1)
myUQUA = uqua(assayData(data_NOISEQ)$exprs, long = lengthuse, lc = 0.5, k = 0)
myTMM = tmm(assayData(data_NOISEQ)$exprs, long = 1000, lc = 0)

############
## STEP 3.1: DESEQ2
############


#BiocManager::install("DESeq2")
# https://mac.r-project.org/tools/
# sudo xcode-select --install
library(DESeq2)
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
setwd("~/Dropbox/WORKING_KAUST/LECTURING/B3XX_SettinngPipelines/DATA_Sets/GSE198256")

load("GSE198256_step1.Rda") # data_NOISEQ,GSE198256_count_filt,annotgene_ord,file="GSE198256_step1.Rda")

############
# STEP 3.1.1: SET THE CLASS
############

GSE198256_DESeq2 <- DESeqDataSetFromMatrix(countData = GSE198256_count_filt,
                              colData = pData(data_NOISEQ),
                              design = ~ Group)
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

plotMA(GSE198256_res, ylim=c(-2,2))
#Interpretation?


lfcShrink(GSE198256_DESeq2_F,coef=c("Group_Healthy_vs_Covid193Mo"))
res_lfcShrink <- lfcShrink(GSE198256_DESeq2_F,coef=c("Group_Covid196Mo_vs_Covid193Mo"))

plotMA(res_lfcShrink, ylim=c(-2,2))
# Why to shrink: it looks at the largest fold changes that are not due 
# to low counts and uses these to inform a prior distribution. 
# So the large fold changes from genes with lots of statistical information are 
# not shrunk, while the imprecise fold changes are shrunk. This allows you to 
# compare all estimated LFC across experiments, for example, which is not really
# feasible without the use of a prior. 
# Michael Love https://support.bioconductor.org/p/77461/



## STEP 3.1.4: Define questions

GSE198256_DESeq2_F<- DESeq(GSE198256_DESeq2_F)
#res <- results(GSE198256_DESeq2_F, contrast=c('factorName','numeratorLevel','denominatorLevel'))
res <- results(GSE198256_DESeq2_F, contrast=c("Group","Healthy","Covid19AI"))
res
resultsNames(GSE198256_DESeq2_F)


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

GSE198256_DESeq2_F <- DESeq(GSE198256_DESeq2_F, test="LRT", reduced=~1)
GSE198256_DESeq2_res_LRT <- results(GSE198256_DESeq2_F)
GSE198256_DESeq2_res_LRT
res <- results(GSE198256_DESeq2_res_LRT)

# Technical replicates?

# How to interpret the results?

plotCounts(GSE198256_DESeq2_F, gene="100287102", intgroup="Group")

## STEP 3.1.4: QC??

# How do visualize?
vsd <- vst(GSE198256_DESeq2_F, blind=FALSE)
rld <- rlog(GSE198256_DESeq2_F, blind=FALSE)
head(assay(vsd), 3)

# heatmap
library("pheatmap")
select <- order(rowMeans(counts(GSE198256_DESeq2_F,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(GSE198256_DESeq2_F)[,c("Group")])
colnames(df) <- "Group"

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE)

# PCA
plotPCA(vsd, intgroup=c("Group"))


############
# STEP 3. BIS: NORMALIZATION AND DIFFERENTIAL EXPRESSION BASED ON LIMMA
############

# Starts at Page 70
# https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

# In the limma approach to RNA-seq, read counts are converted to log2-counts-per-million
# (logCPM) and the mean-variance relationship is modeled either with precision weights
# or with an empirical Bayes prior trend. The precision weights approach is called 
# “voom” and the prior trend approach is called “limma-trend”

require(limma)
require(edgeR)
dge <- DGEList(counts=GSE198256_count)
design <- model.matrix(~ pDataUSE[,1] )

# Filter
keep <- filterByExpr(dge, design=design)
dge <- dge[keep,,keep.lib.sizes=FALSE]

# Normalization
dge <- calcNormFactors(dge)

############
# STEP 3.1 LIMMA: TREND

logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))

############
# STEP 3.2 LIMMA: VOOM

v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))


###########
#### EXERCISE: COMPARE
###########

# How do we compare limma trend vs limma voom?


# How do we compare DESeq2 vs limma trend?


# How do we compare DESeq2 vs limma voom?



############
# STEP 4: BIOLOGICAL INTERPRETATION
############

# Gene Set Enrichment Analysis.
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

# Review Example:
# https://rpubs.com/jrgonzalezISGlobal/enrichment

