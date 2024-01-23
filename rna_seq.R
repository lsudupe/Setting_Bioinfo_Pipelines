

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
data_NOISEQ
