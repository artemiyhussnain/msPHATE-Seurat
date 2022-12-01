library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library (ggplot2)

### RESET ###
wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat/'
med <- readRDS(file = paste(wdir, 'med_msphate.RDS', sep = ''))

FeaturePlot(med, reduction = 'msphate', features = 'CC.Difference', cols = c('red', 'green')) + 
  coord_cartesian(xlim = c(-0.01, 0.02), ylim = c(-0.01, 0.02)) + 
  FeaturePlot(med, reduction = 'umap', features = 'CC.Difference', cols = c('red', 'green'))
# Looks like cycling cells cluster much more clearly in umap

FeaturePlot(med, reduction = 'msphate', features = 'CC.Difference', cols = c('red', 'green')) + 
  coord_cartesian(xlim = c(-0.01, 0.02), ylim = c(-0.01, 0.02)) + DimPlot(med, reduction = 'msphate')
# Unclear, rerunning msPHATE clustering with finer clusters
# How will I then recombine them objectively, without reference to Jansky naming?
# I think, manually look through top 10 DE for each cluster for marker genes

# Kind of a standard workflow:
HoverLocator(DimPlot(med, reduction = 'msphate', pt.size = 1), 
             information = FetchData(med, vars = 'msphate_clusters'))

med <- FindVariableFeatures(med, nfeatures = 5000)
top_de <- head(VariableFeatures(med), 10)
LabelPoints(VariableFeaturePlot(med), points = top_de, repel = TRUE)

Idents(med) <- 'msphate_clusters'
msphate_varmarkers <- FindAllMarkers(med, features = VariableFeatures(med))
msphate_varmarkers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top_de
DoHeatmap(med, features = top_de$gene, label=FALSE)

# Resetting clusters to jansky ones to compare heatmap
Idents(med) <- 'jansky_idents'
jansky_varmarkers <- FindAllMarkers(med, features = VariableFeatures(med))
jansky_varmarkers %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC) -> top_de
DoHeatmap(med, features = top_de$gene, label=FALSE)
# I can see that a lot of genes co-express. Is that not a better representation? I'll try swapping in msphate

## Based on DimPlot and Heatmap:
## 0 - unsure, will keep as is
## 1 - good
## 2 - good
## 3 + 5 merge + small subset of 2
## 6 - good
## 22 - good
## 10 + 29 + 45 + 72 merge
## 82 - good
## 336 + 637 + 912 merge
## 1814 + 1838 merge
## 1893 + 1914 + 5472 merge
## 2701 - good
## 3784 + 3894 + 3927 + 4607 + 8171 merge
## 4829 - good
## 5598 - good
## 7097 - good
Idents(med) <- 'msphate_clusters'

msphate_orig_idents <- levels(med)
med$msphate_orig_idents <- Idents(med)
msphate_new_idents <- c('1', '2', '3', '4', '4', '5', '6', '7', '6', '6', '6', '8', '9', '9', '9', '10', 
                     '10', '11', '11', '12', '13', '13', '13', '13', '14', '11', '15', '16', '13')
names(msphate_new_idents) <- levels(med)
med <- RenameIdents(med, msphate_new_idents)
med$msphate_new_clusters <- Idents(med)

HoverLocator(DimPlot(med, reduction = 'msphate', pt.size = 1), 
             information = FetchData(med, vars = 'msphate_new_clusters'))
msphate_varmarkers2 <- FindAllMarkers(med, features = VariableFeatures(med))
msphate_varmarkers2 %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC) -> top_de
DoHeatmap(med, features = top_de$gene) + NoLegend()
# Merging is taking me nowhere. Maybe I'll cluster more coarsely and separate

# Now for another round of merging
## 1 - 1
## 2 - 2
## 3 - 3
## 4 - 3
## 5 - 3
## 6 - 1
## 7 - 4
## 8 - 5
## 9 - 6
## 10 - 7
## 11 - 4
## 12 - 4
## 13 - 1
## 14 - 3
## 15 - 6
## 16 - 6
Idents(med) <- 'msphate_new_clusters'

msphate_orig_idents2 <- levels(med)
med$msphate_orig_idents2 <- Idents(med)
msphate_new_idents2 <- c('1', '2', '3', '3', '3', '1', '4', '5', '6', '7', '4', '4', '1', '3', '6', '6')
names(msphate_new_idents2) <- levels(med)
med <- RenameIdents(med, msphate_new_idents2)
med$msphate_new_clusters2 <- Idents(med)

HoverLocator(DimPlot(med, reduction = 'msphate', pt.size = 1), 
             information = FetchData(med, vars = 'msphate_new_clusters2'))
msphate_varmarkers3 <- FindAllMarkers(med, features = VariableFeatures(med))
msphate_varmarkers3 %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC) -> top_de
DoHeatmap(med, features = top_de$gene) + NoLegend()