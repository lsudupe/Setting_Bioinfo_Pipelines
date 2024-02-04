packages <- c("limma","stats","methods","RColorBrewer","clusterProfiler","GGally",
              "network","clustree","readxl","org.Hs.eg.db",
              "org.Mm.eg.db","cluster","factoextra","STRINGdb","WebGestaltR","stringr",
              "AnnotationDbi","ComplexHeatmap","GO.db","GetoptLong","bigstatsr","colorRamp2",
              "cowplot","doParallel","dplyr","foreach","ggdendro","ggnewscale","ggplot2",
              "ggtree","ggwordcloud","grid","httr","jsonlite","parallel","patchwork","pbapply",
              "reshape2","rgl","seriation","simplifyEnrichment","slam","tidyverse","umap",
              "utils","grDevices")

new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]

install.packages(new.packages)