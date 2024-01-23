

# libraries
library("NOISeq") #https://rstudio-pubs-static.s3.amazonaws.com/525119_64c1fe6e1a514b89a1ef26d23bf4aae3.html
library("tidyverse")
library("GEOquery")

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE198256", "file=GSE198256_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# read metadata
#meta <- read.table("./data/GSE198256/SraRunTable.txt" , sep=",", header = TRUE)
gds <- getGEO("GSE198256")
Meta_GSE198256 <- pData(gds$GSE198256_series_matrix.txt.gz@phenoData)
Meta_GSE198256 <- Meta_GSE198256[,c("title","source_name_ch1","characteristics_ch1",
                                    "characteristics_ch1.1","description","cell type:ch1","disease state:ch1")]

