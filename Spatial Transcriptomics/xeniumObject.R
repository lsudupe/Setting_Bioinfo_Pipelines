#XENIUM ANALYSIS SPATIAL TRANSCRIPTOMICS


## Libraries
library('Seurat')
library(ggplot2)
library("ggpubr")
library("scran")
library("tibble")
library(dplyr)

getwd()


## Load the Xenium data
path <- "/ibex/scratch/medinils/spatialProyect/data/xenium/"
xenium.obj_x <- LoadXenium(path, fov = "fov", assay = "Xenium")

## Save object
xenium <- xenium.obj_x
saveRDS(xenium, './objects/xenium.rds')
xenium <- readRDS('./objects/xenium.rds')

######QC########
## Plot features
VlnPlot(xenium, features = c("nFeature_Xenium", "nCount_Xenium"), pt.size = 0)
ImageDimPlot(xenium, fov="fov", molecules =c("Gad1", "Sst"), nmols=20000)
ImageFeaturePlot(xenium, features = c("nFeature_Xenium", "nCount_Xenium"), size = 0.75, cols = c("white", "red"))

## Extract metadata
metadata <- xenium@meta.data


## nFeature_Xenium statistics
features_stats <- summary(metadata$nFeature_Xenium)
features_sd <- sd(metadata$nFeature_Xenium)
features_mean <- mean(metadata$nFeature_Xenium)

## nCount_Xenium statistics
counts_stats <- summary(metadata$nCount_Xenium)
counts_sd <- sd(metadata$nCount_Xenium)
counts_mean <- mean(metadata$nCount_Xenium)

## Histograms with quantiles nFeature_Xenium
a <- ggplot(metadata, aes(x = nFeature_Xenium)) +
  geom_histogram(bins = 50, fill = 'blue', alpha = 0.7) +
  geom_vline(xintercept = features_mean, color = "darkblue", linetype = "dashed", size = 1) +
  geom_vline(xintercept = features_mean + features_sd, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = features_mean - features_sd, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = quantile(metadata$nFeature_Xenium, probs = c(0.01, 0.99)), color = "green", linetype = "dashed", size = 1) +
  labs(title = "Distribución de nFeature_Xenium", x = "Número de Features", y = "Frecuencia")

## Histograms with quantiles nCount_Xenium
b <-ggplot(metadata, aes(x = nCount_Xenium)) +
  geom_histogram(bins = 50, fill = 'red', alpha = 0.7) +
  geom_vline(xintercept = counts_mean, color = "darkred", linetype = "dashed", size = 1) +
  geom_vline(xintercept = counts_mean + counts_sd, color = "blue", linetype = "dashed", size = 1) +
  geom_vline(xintercept = counts_mean - counts_sd, color = "blue", linetype = "dashed", size = 1) +
  geom_vline(xintercept = quantile(metadata$nCount_Xenium, probs = c(0.03, 0.97)), color = "green", linetype = "dashed", size = 1) +
  labs(title = "Distribución de nCount_Xenium", x = "Número de Counts", y = "Frecuencia")

fig1 <- ggarrange(a, b,
                  ncol = 2, nrow = 1)
fig1

## Spatial plots before and after

ImageFeaturePlot(xenium, features = c("nFeature_Xenium_log10", "nFeature_Xenium" ), size = 0.75, cols = c("white", "red"))
table(is.finite(xenium@meta.data$nFeature_Xenium_log10)) #FALSE 49  TRUE  36553 

xenium_subset <- subset(xenium, subset = (nCount_Xenium > 46 & nFeature_Xenium > 18))

table(is.finite(xenium_subset@meta.data$nFeature_Xenium_log10)) # TRUE 35440 
ImageFeaturePlot(xenium_subset, features = c("nFeature_Xenium" , "nCount_Xenium"), size = 0.75, cols = c("blue", "yellow"))


######NORMALIZATION########
xenium_normalize <- xenium_subset

xenium_normalize <- NormalizeData(xenium_normalize)
xenium_normalize <- ScaleData(xenium_normalize)
VariableFeatures(xenium_normalize) <- rownames(xenium_normalize)

a <- ImageFeaturePlot(xenium_subset, features = c("Bhlhe40"), size = 0.75, cols = c("blue", "yellow"))
b <- ImageFeaturePlot(xenium_normalize, features = c("Bhlhe40"), size = 0.75, cols = c("blue", "yellow"))
c <- ImageFeaturePlot(xenium_subset, features = c("Gpr17"), size = 0.75, cols = c("blue", "yellow"))
d <- ImageFeaturePlot(xenium_normalize, features = c("Gpr17"), size = 0.75, cols = c("blue", "yellow"))


fig2 <- ggarrange(a,b,c, d,
                  ncol = 2, nrow = 2)
fig2


######ANALYSIS########
xenium_normalize_processed <- FindVariableFeatures(xenium_normalize)
xenium_normalize_processed <- RunPCA(xenium_normalize_processed, npcs = 30, features = rownames(xenium_normalize_processed))
xenium_normalize_processed <- RunUMAP(xenium_normalize_processed, dims = 1:30)
xenium_normalize_processed <- FindNeighbors(xenium_normalize_processed, reduction = "pca", dims = 1:30)
xenium_normalize_processed <- FindClusters(xenium_normalize_processed, resolution = 0.3)

## Clusters
DimPlot(xenium_normalize_processed)
FeaturePlot(xenium_normalize_processed, features = c("Gad1", "Bdnf"))
ImageDimPlot(xenium_normalize_processed, cols = "polychrome", size = 0.6)


######DE########
xenium_16vs17 <- FindMarkers(xenium_normalize_processed, ident.1 = "16", ident.2 = "17")
head(xenium_16vs17, n = 10)

ImageFeaturePlot(xenium_normalize, features = c("Bcl11b", "Cpne4"), size = 0.75, cols = c("blue", "yellow"))

######SEGMENTATION########
cropped.coords <- Crop(xenium_normalize_processed[["fov"]], x = c(1200, 2900), y = c(3750, 4550), coords = "plot")
xenium_normalize_processed[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium_normalize_processed[["zoom"]]) <- "segmentation"
ImageDimPlot(xenium_normalize_processed, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("Bcl11b", "Cpne4"), nmols = 10000)


