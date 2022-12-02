library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library (ggplot2)

wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat/'
med <- readRDS(file = paste(wdir, 'med_msphate.RDS', sep = ''))

DimPlot(med, reduction = 'msphate')
Idents(med) <- 'msphate_clusters'
clus <- subset(med, subset = msphate_clusters %in% c('1'))
# OR
clus <- subset(med, subset = jansky_idents %in% c('Bridge'))

clus <- FindVariableFeatures(object = clus)
clus <- ScaleData(clus)
clus <- RunPCA(object = clus)
clus <- FindNeighbors(object = clus)
clus <- FindClusters(object = clus, resolution = 0.3)
clus <- RunTSNE(object = clus, features = VariableFeatures(clus))
DimPlot(object = clus, reduction = "tsne", pt.size=1.5)

msphate_varmarkers_clus_zoom <- FindAllMarkers(clus, features = VariableFeatures(clus))
msphate_varmarkers_clus_zoom %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC) -> top_de
DoHeatmap(clus, features = top_de$gene) + NoLegend()

DimPlot(object = clus, reduction = "umap", pt.size=1.5)


## msPHATE of all data, some manual merging based on heatmap and umap
## Zoom into one msPHATE cluster, recluster with msPHATE or seurat, visualise with tSNE
# Test if zooming into umap and tsne clusters is always trash
# Check difference between reclustering a global msPHATE cluster with msPHATE or Seurat
# Check if tSNE is best vis method for reclustered global msPHATE cluster

# Hypothesis: good workflow is to do global clustering with msPHATE, zoom into a cluster, recluster (with msPHATE OR seurat??), vis with tSNE
