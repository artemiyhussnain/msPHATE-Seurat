library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(tibble)
library(cowplot)
library(purrr)
library(SeuratDisk)
library(RColorBrewer)

wdir <- '~/path/to/data/folder'
data_name <- 'data_file.RDS(.gz)'
data <- readRDS(paste(wdir, data_name, sep = ''))

# Generate scale.data first using https://satijalab.org/seurat/articles/pbmc3k_tutorial.html (https://satijalab.org/seurat/articles/essential_commands.html for summary)

write.table(GetAssayData(data, 'scale.data'), file = paste(wdir, 'counts.csv', sep = ''), sep = ',', row.names = FALSE, col.names = FALSE)
write.table(colnames(GetAssayData(data, 'scale.data')), file = paste(wdir, 'cell_names.txt', sep = ''), sep = ',', row.names = FALSE, col.names = FALSE)
write.table(rownames(GetAssayData(data, 'scale.data')), file = paste(wdir, 'gene_names.txt', sep = ''), sep = ',', row.names = FALSE, col.names = FALSE)

# Now run msphate.py