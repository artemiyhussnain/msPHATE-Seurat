library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

data_file = '~/Jansky/data/adrenal_medulla_Seurat.RDS.gz'
wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat'
use_scaled = TRUE

medulla_seurat <- readRDS(file = data_file)
if (use_scaled) {
  write.table(GetAssayData(medulla_seurat, 'scale.data'), 
              file = paste(wdir, 'counts_scaled.csv', sep = "/"), 
              sep = ',', row.names = FALSE, col.names = FALSE)
} else {
  writeMM(medulla_seurat@assays$RNA@counts, paste(wdir, 'counts.mtx', sep = ""))
}

write.table(colnames(medulla_seurat), file = paste(wdir, 'cell_names.txt', sep = "/"),
            sep = ',', row.names = FALSE, col.names = FALSE)
write.table(rownames(medulla_seurat), file = paste(wdir, 'gene_names.txt', sep = "/"),
            sep = ',', row.names = FALSE, col.names = FALSE)