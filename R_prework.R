library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

use_scaled = TRUE
do_subset <- FALSE
# Still working on extracting data based on criteria


wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat/'
med <- readRDS('~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat/adrenal_medulla_Seurat.RDS.gz')

# if (do_subset == TRUE) {
#   use_scaled <- FALSE
#   med <- readRDS(paste(wdir, 'med_msphate.RDS', sep = ''))
# clus <- colnames(subset(med, subset = msphate_clusters == '1'))
# scale_data <- GetAssayData(med, 'scale.data')
# scale_subset <- scaled_data[colnames(scaled_data) %in% clus]
# write.table(scale_subset, file = paste(wdir, 'counts.csv', sep = ''), 
#             sep = ',', row.names = FALSE, col.names = FALSE)
# write.table(clus, file = paste(wdir, 'cell_names.txt', sep = ''),
#             sep = ',', row.names = FALSE, col.names = FALSE)
# write.table(rownames(scale_subset), file = paste(wdir,
#                                                                     'gene_names.txt', sep = ''),
#             sep = ',', row.names = FALSE, col.names = FALSE)
# }

if (use_scaled) {
  write.table(GetAssayData(med, 'scale.data'), 
              file = paste(wdir, 'counts.csv', sep = ''), 
              sep = ',', row.names = FALSE, col.names = FALSE)
} else {
  writeMM(medulla_seurat@assays$RNA@counts, paste(wdir, 'counts.mtx', sep = ''))
}

write.table(colnames(GetAssayData(med, 'scale.data')), file = paste(wdir,
                                                                    'cell_names.txt', sep = ''),
            sep = ',', row.names = FALSE, col.names = FALSE)
write.table(rownames(GetAssayData(med, 'scale.data')), file = paste(wdir,
                                                                    'gene_names.txt', sep = ''),
            sep = ',', row.names = FALSE, col.names = FALSE)