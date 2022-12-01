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

msphate_orig_idents2 <- levels(med)
med$msphate_orig_idents <- Idents(med)
msphate_new_idents2 <- c('1', '2', '1', '3', '3', '4')
names(msphate_new_idents2) <- levels(med)
med <- RenameIdents(med, msphate_new_idents2)
med$msphate_new_clusters2 <- Idents(med)
Idents(med) <- 'msphate_new_clusters2' # Necessary to run findallmarkers correctly
HoverLocator(DimPlot(med, reduction = 'msphate', pt.size = 1), 
             information = FetchData(med, vars = 'msphate_new_clusters2'))

msphate_varmarkers2 <- FindAllMarkers(med, features = VariableFeatures(med))
msphate_varmarkers2 %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC) -> top_de
DoHeatmap(med, features = top_de$gene) + NoLegend()

# Let's look at the cell cycle distribution again
FeaturePlot(med, reduction = 'msphate', features = 'CC.Difference', cols = c('red', 'green')) + 
  coord_cartesian(xlim = c(-0.01, 0.02), ylim = c(-0.01, 0.02))
# Better in umap
FeaturePlot(med, reduction = 'umap', features = 'CC.Difference', cols = c('red', 'green'))

# Let's look more closely at cluster 2, the biggest one (THIS WORKSS!!!)
clus <- subset(med, subset = msphate_new_clusters2 == '2')
scaled_data <- GetAssayData(clus, 'scale.data')
wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat/'
write.table(scaled_data, 
            file = paste(wdir, 'counts.csv', sep = ''), 
            sep = ',', row.names = FALSE, col.names = FALSE)
write.table(rownames(scaled_data), file = paste(wdir, 'gene_names.txt', sep = ''), 
            sep = ',', row.names = FALSE, col.names = FALSE)
write.table(colnames(scaled_data), file = paste(wdir, 'cell_names.txt', sep = ''), 
            sep = ',', row.names = FALSE, col.names = FALSE)

msPHATE_embedding <- read.csv(paste(wdir, 'msPHATE_embedding.csv', sep = ''), 
                              header=FALSE,
                              row.names = 1, 
                              col.names = c('', 'msPHATE_1', 'msPHATE_2'))
msPHATE_embedding <- as.matrix(msPHATE_embedding)
clus[['msphate']] <- CreateDimReducObject(embeddings = msPHATE_embedding, key = 'msPHATE_')
msPHATE_clusters <- as.matrix(read.csv(paste(wdir, 'msPHATE_clusters.csv', sep = ''), header=FALSE))
msPHATE_clusters <- factor(msPHATE_clusters)
clus$msphate_clusters <- msPHATE_clusters
clus$jansky_idents <- Idents(clus)
Idents(clus) <- 'msphate_clusters'
DimPlot(clus, reduction = 'msphate')

msphate_new_idents <- c('1', '2', '3', '4', '5')
names(msphate_new_idents) <- levels(clus)
clus <- RenameIdents(clus, msphate_new_idents)
clus$msphate_clusters <- Idents(clus)
Idents(clus) <- 'msphate_clusters'
DimPlot(clus, reduction = 'msphate')

msphate_varmarkers_clus2 <- FindAllMarkers(clus, features = VariableFeatures(clus))
msphate_varmarkers_clus2 %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC) -> top_de
DoHeatmap(clus, features = top_de$gene) + NoLegend()

DimPlot(clus, reduction = 'msphate', split.by = 'timepoint2')
