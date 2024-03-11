# ATAC https://www.10xgenomics.com/datasets/8k-adult-mouse-cortex-cells-atac-v2-chromium-controller-2-standard

# Libraries

library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)
library(EnsDb.Mmusculus.v79)


#####################################
##### PRE-PROCESSING WORKFLOW #######
#####################################

#loading the data
counts <- Read10X_h5(filename = "./homework/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "./homework/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_singlecell.csv",
  header = TRUE,
  row.names = 1
)

#
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = './homework/8k_mouse_cortex_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

#creating the seurat object
atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "atac",
  meta.data = metadata
)

granges(atac)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to mm10
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(atac) <- annotations

#####################################
##### QUALITY CONTROL METRICS #######
#####################################

# compute nucleosome signal score per cell
atac <- NucleosomeSignal(object = atac)
atac$nucleosome_group <- ifelse(atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
pdf("./fragmenthistogram.pdf")
FragmentHistogram(object = atac, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
dev.off()

saveRDS(atac,"./objects/raw_atac.rds")
atac <- readRDS("./objects/raw_atac.rds")

#####################################
############# ANALYSIS ##############
#####################################

# add blacklist ratio and fraction of reads in peaks
atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$peak_region_fragments

pdf("./results/violin.pdf")
VlnPlot(
  object = atac,
  features = c( 'peak_region_fragments',
                 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()


## Normalization and linear dimensional reduction

atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 'q0')
atac <- RunSVD(object = atac)

pdf("./results/cor.pdf")
DepthCor(atac)
dev.off()

#visualization of significatives principal components
pdf("./results/elbow.pdf")
ElbowPlot(object = atac, 
          ndims = 40, reduction ="lsi")
dev.off()


## Non-linear dimension reduction and clustering

atac <- RunUMAP(
  object = atac,
  reduction = 'lsi',
  dims = 2:30
)
atac <- FindNeighbors(
  object = atac,
  reduction = 'lsi',
  dims = 2:30
)
atac <- FindClusters(
  object = atac,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

pdf("./results/umap.pdf")
DimPlot(object = atac, label = TRUE) + NoLegend()
dev.off()

saveRDS(atac,"./objects/atac.rds")


## Create a gene activity matrix
gene.activities <- GeneActivity(atac)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
atac <- NormalizeData(
  object = atac,
  assay = 'RNA'
)



#####################################
## INTEGRATING WITH scRNAseq DATA ###
#####################################

#loading scRNAseq data
sc <- readRDS('/home/sotoorda/single_cell/objects/sc.rds')
sc
DefaultAssay(atac) <- 'RNA'

#cross-modality integration and label transfer,
transfer.anchors <- FindTransferAnchors(
  reference = sc,
  query = atac,
  reduction = 'cca',
  dims = 1:40
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata =  sc@meta.data$predicted.subclass,
  weight.reduction = atac[['lsi']],
  dims = 2:30
)

brain <- AddMetaData(object = atac, metadata = predicted.labels)

saveRDS(brain,"./objects/atac_predicted.rds")

#visualizing integrated data

plot1 <- DimPlot(sc, group.by = 'predicted.subclass', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(brain, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

pdf("./results/prediction.pdf", width = 10)
print(plot1 + plot2)
dev.off()

#####################################
######### JACKARD COMPARISON ########
#####################################

library(dplyr)

## read data
gex_atac <- readRDS("./objects/atac_predicted.rds")
multiome <- readRDS("./objects/mouse_brain_multiome_Integrated.rds")

Idents(multiome) <- "predicted.subclass"
Idents(gex_atac) <- "predicted.id"

calculate_jaccard <- function(set1, set2) {
  return(length(intersect(set1, set2)) / length(union(set1, set2)))
}

## Extracting cell types from both datasets
gex_atac_celltypes <- Idents(gex_atac)
multiome_celltypes <- Idents(multiome)

## Calculate Jaccard Index for each cell type present in both datasets
cell_types <- intersect(unique(gex_atac_celltypes), unique(multiome_celltypes))
jaccard_indices <- data.frame(CellType = cell_types, JaccardIndex = NA_real_)

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  jaccard_indices$JaccardIndex[i] <- calculate_jaccard(
    which(gex_atac_celltypes == cell_type),
    which(multiome_celltypes == cell_type)
  )
}

## Review the Jaccard indices
jaccard_indices

## Plot the Jaccard indices
pdf("./results/jackard.pdf")
ggplot(jaccard_indices, aes(x = CellType, y = JaccardIndex, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Jaccard Index by Cell Type",
    x = "Cell Type",
    y = "Jaccard Index",
    fill = "Cell Type"
  )
dev.off()