library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

data_file = '~/Jansky/data/adrenal_medulla_Seurat.RDS.gz'
wdir = '~/Jansky/ms_phate/msphate_jansky/timepoint_phateplot/'

medulla_seurat <- readRDS(file = data_file)
writeMM(medulla_seurat@assays$RNA@counts, paste(wdir, 'counts.mtx', sep = ""))
write.table(colnames(medulla_seurat), file = paste(wdir, 'cell_names.txt', sep = "") , sep = ',', row.names = FALSE, col.names = FALSE)
write.table(rownames(medulla_seurat), file = paste(wdir, 'gene_names.txt', sep = ""), sep = ',', row.names = FALSE, col.names = FALSE)

c18 <- colnames(subset(medulla_seurat, subset = timepoint2 == 'CS18'))
cs19 <- colnames(subset(medulla_seurat, subset = timepoint2 == 'CS19'))
cs20 <- colnames(subset(medulla_seurat, subset = timepoint2 == 'CS20'))
cs23 <- colnames(subset(medulla_seurat, subset = timepoint2 == 'CS23'))
pcw11 <- colnames(subset(medulla_seurat, subset = timepoint2 == '11pcw'))
pcw14 <- colnames(subset(medulla_seurat, subset = timepoint2 == '14pcw'))
pcw17 <- colnames(subset(medulla_seurat, subset = timepoint2 == '17pcw'))
timepoint_plot <- DimPlot(medulla_seurat, reduction = 'umap', split.by = 'timepoint2', combine = TRUE)
