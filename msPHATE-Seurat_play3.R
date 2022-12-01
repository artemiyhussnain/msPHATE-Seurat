library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library (ggplot2)

### RESET ###
wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat/'
med <- readRDS(file = paste(wdir, 'med_msphate.RDS', sep = ''))
### RESET ###
Idents(med) <- 'msphate_clusters'
med$msphate_orig1_clusters <- Idents(med)

msphate_new_idents <- c('1', '2', '3', '4', '5', '6', '7')
names(msphate_new_idents) <- levels(med)
med <- RenameIdents(med, msphate_new_idents)
med$msphate_clusters <- Idents(med)

HoverLocator(DimPlot(med, reduction = 'msphate', pt.size = 1), 
             information = FetchData(med, vars = 'msphate_clusters'))

msphate_varmarkers <- FindAllMarkers(med, features = VariableFeatures(med))
msphate_varmarkers %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC) -> top_de
DoHeatmap(med, features = top_de$gene) + NoLegend()

### Ok now I'll try to invert the data (specifically mix up the cell and gene row and column in msPHATE)
# I knew it wouldn't change anything - it didn't