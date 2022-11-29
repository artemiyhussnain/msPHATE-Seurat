library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat'
medulla_seurat <- readRDS(file = data_file)
msPHATE_embedding <- read.csv(paste(wdir, 'msPHATE_embedding.csv', sep = '/'), 
                              header=FALSE,
                              row.names = 1, 
                              col.names = c('', 'msPHATE_1', 'msPHATE_2'))
msPHATE_embedding <- as.matrix(msPHATE_embedding)
medulla_seurat[['msphate']] <- CreateDimReducObject(embeddings = msPHATE_embedding, 
                                                    key = 'msPHATE_')

msPHATE_clusters <- as.matrix(read.csv(paste(wdir, 'msPHATE_clusters.csv', sep = '/'), 
                             header=FALSE))
msPHATE_clusters <- factor(msPHATE_clusters)
medulla_seurat$msphate_clusters <- msPHATE_clusters