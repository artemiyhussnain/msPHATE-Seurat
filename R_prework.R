library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

use_scaled = TRUE

med <- readRDS('adrenal_medulla_Seurat.RDS.gz')

if (use_scaled) {
  write.table(GetAssayData(med, 'scale.data'), 
              file = 'counts.csv', 
              sep = ',', row.names = FALSE, col.names = FALSE)
} else {
  writeMM(medulla_seurat@assays$RNA@counts, 'counts.mtx')
}

write.table(colnames(GetAssayData(med, 'scale.data')), file = 'cell_names.txt',
            sep = ',', row.names = FALSE, col.names = FALSE)
write.table(rownames(GetAssayData(med, 'scale.data')), file = 'gene_names.txt',
            sep = ',', row.names = FALSE, col.names = FALSE)