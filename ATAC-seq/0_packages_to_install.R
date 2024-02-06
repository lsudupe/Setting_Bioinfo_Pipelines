
# BiocManager package
if(!require('BiocManager')){
  install.packages("BiocManager")
}

# bioconductor
BiocManager::install()


## Then, run the following code to install the necessary packages. 
# This can take some time.
packages <- c(
  "universalmotif", "MotifDb", 
  "TFBSTools", "JASPAR2020", "motifmatchr", 
  "BSgenome.Hsapiens.UCSC.hg38", 
  "GenomicRanges", "ChIPseeker",
  "clusterProfiler", "chipenrich",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "rtracklayer", "DESeq2")

# avoid installing already present packaged
packages <- setdiff(packages, rownames(installed.packages()))
if(length(packages) > 0){
  BiocManager::install(packages)
}

# ggplot2
install.packages('ggplot2')
