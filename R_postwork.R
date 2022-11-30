library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

msPHATE_embedding <- read.csv('msPHATE_embedding.csv', 
                              header=FALSE,
                              row.names = 1, 
                              col.names = c('', 'msPHATE_1', 'msPHATE_2'))
msPHATE_embedding <- as.matrix(msPHATE_embedding)
med[['msphate']] <- CreateDimReducObject(embeddings = msPHATE_embedding, 
                                                    key = 'msPHATE_')

msPHATE_clusters <- as.matrix(read.csv('msPHATE_clusters.csv', 
                             header=FALSE))
msPHATE_clusters <- factor(msPHATE_clusters)
med$msphate_clusters <- msPHATE_clusters