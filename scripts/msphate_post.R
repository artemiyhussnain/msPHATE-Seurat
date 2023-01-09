msPHATE_embedding <- read.csv(paste(wdir, 'msPHATE_embedding.csv', sep = ''), header=FALSE, row.names = 1, col.names = c('', 'msPHATE_1', 'msPHATE_2'))
msPHATE_embedding <- as.matrix(msPHATE_embedding)
data[['msphate']] <- CreateDimReducObject(embeddings = msPHATE_embedding, key = 'msPHATE_')

msPHATE_clusters <- as.matrix(read.csv(paste(wdir, 'msPHATE_clusters.csv', sep = ''), header=FALSE))
msPHATE_clusters <- factor(msPHATE_clusters)
data$msphate_clusters <- msPHATE_clusters

DimPlot(data, reduction = 'msphate', group.by = 'msphate_clusters', pt.size = 1)
