import multiscale_phate as mp
import numpy as np
import pandas as pd
import scprep
import os
import matplotlib.pyplot as plt
import pickle

# Git test
# Seurat: unsure how to find out QC from seurat object
# Assuming tutorial defaults
# Seurat: FindNeighbours used 20 PC's
# Seurat: FindClusters (resolution = 1.2)
# Shunya's report: min 500 genes, min 1000 counts, max 2.5% mt - where are QC params from?
min_reads = 200
wdir = '~/Jansky/ms_phate/msphate_jansky/timepoint_phateplot'
mincells = 3
npca = 20
gran = 0.1
run_multiple_embeddings=False
multiple_spread = 1
vis_level = 0
clus_level = 3
marker_dict = {'SCPs': ['SOX10', 'PLP1', 'ERBB3'],
               'Neuroblasts': ['ISL1', 'ALK'],
               'Chromaffin cells': ['DBH', 'TH', 'CHGA']}
names = list(marker_dict.keys())
genes = list(marker_dict.values())


mtx_dir = os.path.expanduser(wdir)
mtx_path = os.path.join(mtx_dir, 'counts.mtx')
gene_dir = os.path.expanduser(wdir)
gene_path = os.path.join(gene_dir, 'gene_names.txt')
cell_dir = os.path.expanduser(wdir)
cell_path = os.path.join(cell_dir, 'cell_names.txt')

print('Loading counts matrix...')
data = scprep.io.load_mtx(mtx_path, cell_axis='column',
                          gene_names=gene_path,
                          cell_names=cell_path,
                          sparse=True)

print('Filtering and normalising...')
data = scprep.filter.filter_library_size(data, cutoff=min_reads,
                                         keep_cells='above')
data = scprep.filter.filter_rare_genes(data, min_cells = mincells)
data = scprep.normalize.library_size_normalize(data)
data = np.sqrt(data)
# Can QC'd and normalised counts matrix be imported from Seurat?

print('Running multiscale PHATE...')
mp_op = mp.Multiscale_PHATE(n_pca=npca, granularity=gran, random_state=0)

print('Finding levels...')
levels = mp_op.fit(data)
ax = plt.plot(mp_op.gradient)
ax = plt.scatter(levels, mp_op.gradient[levels], c = 'r', s=100)
plt.savefig(fname='levels.png')
plt.close()

tree = mp_op.build_tree()
tree_clusters = mp_op.get_tree_clusters(levels[-1*clus_level])
scprep.plot.scatter3d(tree, c = tree_clusters, s= 50,
                      fontsize=16, ticks=False, figsize=(10,10))
plt.savefig(fname='tree.png')
plt.close()

print('Generating embedding(s)...')
if run_multiple_embeddings:
    for i in range(vis_level, round(len(levels)/2)+multiple_spread, 2):
        for j in [-1*clus_level, -1*clus_level+multiple_spread, -1*clus_level-multiple_spread]:
            embedding, clusters, sizes = mp_op.transform(visualization_level = levels[i],
                                                         cluster_level = levels[j])
            scprep.plot.scatter2d(embedding, s = 100*np.sqrt(sizes), c = clusters,
                                  fontsize=16, ticks=False,
                                  label_prefix="Multiscale PHATE", figsize=(10,8))
            plt.savefig(fname='embedding_vis' + str(j) + '_clus' + str(-i)+ '.png')
            plt.close()
else:
    embedding, clusters, sizes = mp_op.transform(visualization_level = levels[vis_level],
                                                         cluster_level = levels[-1*clus_level])
    scprep.plot.scatter2d(embedding, s = 100*np.sqrt(sizes), c = clusters,
                                  fontsize=16, ticks=False,
                                  label_prefix="Multiscale PHATE", figsize=(10,8))
    plt.savefig(fname='embedding_vis' + str(vis_level) + '_clus' +
                        str(-1*clus_level)+ '.png')
    plt.close()


print('Finding expression')
expression = pd.DataFrame()
for i in range(len(list(marker_dict.values()))):
        expression[list(marker_dict.values())[i]] = mp_op.get_expression(data[list(marker_dict.values())[i]],
                                                   visualization_level =  levels[vis_level])
# Something's broken here, reintroduce names and genes lists and only do dict kerfuffle at the start?

print('Pickling important files...')
keep = [embedding, clusters, sizes, data]
keep_names = ['embedding', 'clusters', 'sizes', 'data', 'expression']
for i in range(len(keep)):
    pickle_out = open(keep_names[i] + '.pickle', 'wb')
    pickle.dump(keep[i], pickle_out)
    pickle_out.close()

print(data)
print(levels)
print('Done')
