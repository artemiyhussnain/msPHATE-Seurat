library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library (ggplot2)

# Basic visualisation
DimPlot(med, reduction = 'msphate', group.by = 'msphate_clusters')

# Comparing how well Jansky labelling corresponds to msphate
DimPlot(med, reduction = 'msphate', group.by = 'msphate_clusters') + 
  DimPlot(med, reduction = 'msphate') +
  NoLegend()

# Splitting both by timepoint - alternative to pseudotime?
DimPlot(med, reduction = 'msphate', group.by = 'msphate_clusters',
        split.by = 'timepoint2')
DimPlot(med, reduction = 'umap', split.by = 'timepoint2')

# Plotting cell cycle over msphate, but need better definition of 'cycling' - perhaps sum of G2M and S scores?
#DimPlot(med, reduction = 'msphate', group.by = 'msphate_clusters', pt.size = 1) + 
#  DimPlot(med, reduction = 'msphate', group.by = 'Phase', pt.size = 1) +
#  NoLegend()
#DimPlot(med, reduction = 'umap') + 
#  DimPlot(med, reduction = 'umap', group.by = 'Phase') +
#  NoLegend()

# Key marker gene expression on msphate
marker_plots <- FeaturePlot(med, reduction = 'msphate', 
            features = c('SOX10', 'DBH', 'ISL1'), combine=TRUE)
marker_plots[[1]] + coord_cartesian(xlim = c(-0.01, 0.02), ylim = c(-0.01, 0.02)) + NoLegend() +
marker_plots[[2]] + coord_cartesian(xlim = c(-0.01, 0.02), ylim = c(-0.01, 0.02)) + NoLegend() +
marker_plots[[3]] + coord_cartesian(xlim = c(-0.01, 0.02), ylim = c(-0.01, 0.02))

# Same for umap - marker genes cluster equally well
marker_plots <- FeaturePlot(med, reduction = 'umap', features = c('SOX10', 'DBH', 'ISL1'))
marker_plots[[1]] + NoLegend() + marker_plots[[2]] + NoLegend() + marker_plots[[3]]

# Should have done this on top
jansky_idents <- levels(med)
med$jasky_idents <- Idents(med) # Old idents stored under jansky_idents, can be switched into as in next line
Idents(med) <- 'msphate_clusters'
DimPlot(med, reduction = 'msphate') # Now msphate clusters are Idents
jansky_markers <- msphate_markers
msphate_markers <- FindAllMarkers(med)

# Seems that cluster 1 has no clear markers, and clusters 22 and 5472 are the same
msphate_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top_de
DoHeatmap(med, features = top_de$gene) + NoLegend()
DoHeatmap(med, features = top_de$gene[!(top_de$gene %in% c(cc.genes.updated.2019$s.genes,
                                                      cc.genes.updated.2019$g2m.genes))], 
          label=FALSE) # Excluding cc genes

# Merging 22 and 5427, renaming clusters 1-5
# B+lS = Bridge and some late SCPs?
# N+B+clC = Neuroblasts and some Bridge plus connecting and some late Chromaffin cells?
# clS = SCPs including cycling and most late
# C = Chromaffin cells
# lC = Some late chromaffin cells
new_cluster_ids <- c('B+lS', 'N+B+clC', 'clS', 'C', 'C', 'lC')
names(new_cluster_ids) <- levels(med)
med <- RenameIdents(med, new_cluster_ids)
DimPlot(med, reduction = "msphate", pt.size = 1)
