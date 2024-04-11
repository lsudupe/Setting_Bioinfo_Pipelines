setwd('/Users/sotoorda/Documents/bioinfo_course/methylation')

# load packages required for analysis
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)

# creating the annotation
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)


dataDirectory <- "/Users/medinils/Documents/bioinfo_course/methylation/dataDirectory"
list.files(dataDirectory , recursive = TRUE)

# Define el directorio donde tienes tus archivos IDAT
idatPath <- dataDirectory
# Obtén todos los archivos IDAT en la carpeta
idatFiles <- list.files(idatPath, pattern = 'idat.gz$', full.names = TRUE)
# Genera nombres base únicos (sin _Grn o _Red y sin la extensión .idat.gz)
basenames <- unique(sub('_Grn.idat.gz|_Red.idat.gz', '', basename(idatFiles)))
# Crea un data.frame con estos nombres base
targets <- data.frame(Basename = basenames)
# Ahora usa este data.frame en la función read.metharray.exp
rgSet <- read.metharray.exp(base = idatPath, targets = targets)
rgSet
# give the samples descriptive names
targets$Sample_Group <- c('HD', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 
                'HD', 'severe_COVID', 'severe_COVID',
                'HD', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID',
                'HD', 'severe_COVID', 'severe_COVID', 'severe_COVID',
                'HD', 'severe_COVID', 'severe_COVID', 'severe_COVID',
                'HD', 'severe_COVID',
                'HD', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID',
                'HD', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID',
                'HD', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID',
                'HD', 
                'HD', 'severe_COVID', 'severe_COVID', 'severe_COVID', 'severe_COVID')

sampleNames(rgSet) <- targets$Sample_Group



#######################
### QUALITY CONTROL ###
#######################



detP <- detectionP(rgSet)
head(detP)

# examine mean detection p-values across all samples to identify any failed samples

pdf("results/barplot.pdf", width = 10, height = 15,)

pal <- brewer.pal(8,"Dark2")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.6, ylim=c(0,0.0005), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white")

dev.off()


qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
        pdf="results/qcReport.pdf")

# check for poor quality samples
keep <- colMeans(detP) < 0.01
keep

# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet) 

# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)

# visualise what the data looks like before and after normalisation

pdf("results/normalization.pdf", width = 20, height = 15)

par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
      text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))

dev.off()


########
#Data exploration

# MDS plots to look at largest sources of variation

pdf("results/pca.pdf", width = 20, height = 15,)

plotMDS(getM(mSetSq), top=1000, gene.selection="common", labels = NULL,
        col=pal[factor(targets$Sample_Group)])

legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)


dev.off()

################
## FILTERING ###
################

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt



# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])

bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])

#Plotting m and b values
pdf("results/m-values.pdf", width = 20, height = 15)
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)), 
       text.col=brewer.pal(8,"Dark2"))
dev.off()



##################################################
## Probe-wise differential methylation analysis ##
##################################################

targets$Sample_Group <- factor(targets$Sample_Group, levels = c('HD', 'severe_COVID'))
# Crear la matriz de diseño basada en cellType
design <- model.matrix(~ cellType, data = targets)
# Ajustar el modelo lineal usando mVals como tus datos de expresión
fit <- lmFit(mVals, design)
# Aplicar el método eBayes para obtener estadísticas de los coeficientes ajustados
fit2 <- eBayes(fit)
head(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

# get the table of results for the first contrast (HD - severe COVID)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)


# extracting significantly methylated probes
deff.meth = topTable(fit2, coef=ncol(design), number = nrow(mVals), adjust.method = "BY")


# Visualization
# plot the top 10 most significantly differentially methylated CpGs 
par(mfrow=c(2,5))

pdf("results/DMP.pdf", width = 10, height = 15,)

sapply(rownames(deff.meth)[1:10], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})

dev.off()


