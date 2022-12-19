### USAGE ###
# Export counts.csv. gene_names.txt and cell.names.txt from R
# Now you can load in data by restarting
# Find msPHATE clusters at different levels with clus
### USAGE ###

import multiscale_phate as mp
import numpy as np
import pandas as pd
import scprep
import os

wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat'
npca = 20
gran = 0.1

wdir = os.path.expanduser(wdir)
gene_path = os.path.join(wdir, 'gene_names.txt')
cell_path = os.path.join(wdir, 'cell_names.txt')

print('Loading scaled data...')
data_path = os.path.join(wdir, 'counts.csv')
data = scprep.io.load_csv(data_path, cell_axis='column',
                        gene_names=gene_path,
                        cell_names=cell_path)

mp_op = mp.Multiscale_PHATE(n_pca=npca, granularity=gran, random_state=0)
levels = mp_op.fit(data)

def clus(clustering_level, vis_level = 0):
    clus_level = clustering_level
    print('Generating embedding and assigning clusters...')
    embedding, clusters, sizes = mp_op.transform(visualization_level = levels[vis_level],
                                                               cluster_level = levels[-1*clus_level])

    print('Exporting embedding and clusters...')
    msPHATE_embedding = pd.DataFrame(embedding,
                          index = pd.read_csv(cell_path, header=None).iloc[:, 0].to_list())
    msPHATE_embedding.to_csv('msPHATE_embedding.csv', header=False)
    msPHATE_clusters = pd.DataFrame(clusters)
    msPHATE_clusters.to_csv('msPHATE_clusters.csv', header=False, index=False)

    print('\n')
    print(data)
    print('\n')
    print(levels)
    print('\n')
    print(embedding)
    print('\n')
    print(set(clusters))
    print ('\n')
    print(' Clustering level ' + str(clus_level))
    print('\n')
    print('Done')
