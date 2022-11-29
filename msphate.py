import multiscale_phate as mp
import numpy as np
import pandas as pd
import scprep
import os
import matplotlib.pyplot as plt
import pickle

# Seurat: unsure how to find out QC from seurat object
# Assuming tutorial defaults
# Seurat: FindNeighbours used 20 PC's
# Seurat: FindClusters (resolution = 1.2)
# Shunya's report: min 500 genes, min 1000 counts, max 2.5% mt - where are QC params from?
wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat'
do_scaling = False
min_reads = 200
mincells = 3
npca = 20
gran = 0.1
run_multiple_embeddings=False
generate_tree=True
multiple_spread = 1
vis_level = 0
clus_level = 3
marker_dict = {'SCPs': ['SOX10', 'PLP1', 'ERBB3', 'MPZ', 'FOXD3'],
               'Neuroblasts': ['ISL1', 'STMN2', 'NEFM'],
               'Chromaffin cells': ['DBH', 'PHOX2B']}

names = list(marker_dict.keys())
genes = []
for i in marker_dict:
    for j in marker_dict[i]:
        genes.append(j)


wdir = os.path.expanduser(wdir)
gene_path = os.path.join(wdir, 'gene_names.txt')
cell_path = os.path.join(wdir, 'cell_names.txt')

if do_scaling:
    print('Loading counts matrix...')
    data_path = os.path.join(wdir, 'counts.mtx')
    data = scprep.io.load_mtx(data_path, cell_axis='column',
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
else:
    print('Loading scaled data...')
    data_path = os.path.join(wdir, 'counts.mtx')
    data = scprep.io.load_mtx(data_path, cell_axis='column',
                              gene_names=gene_path,
                              cell_names=cell_path,
                              sparse=True)

print('Finding levels...')
mp_op = mp.Multiscale_PHATE(n_pca=npca, granularity=gran, random_state=0)
levels = mp_op.fit(data)
ax = plt.plot(mp_op.gradient)
ax = plt.scatter(levels, mp_op.gradient[levels], c = 'r', s=100)
f_dir = os.path.expanduser(wdir)
f_name = os.path.join(f_dir, 'levels.png')
plt.savefig(fname=f_name)
plt.close()

if generate_tree:
    tree = mp_op.build_tree()
    tree_clusters = mp_op.get_tree_clusters(levels[-1*clus_level])
    scprep.plot.scatter3d(tree, c = tree_clusters, s= 50,
                          fontsize=16, ticks=False, figsize=(10,10))
    f_dir = os.path.expanduser(wdir)
    f_name = os.path.join(f_dir, 'tree.png')
    plt.savefig(fname=f_name)
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
            f_dir = os.path.expanduser(wdir)
            custom_name = 'embedding_vis' + str(i) + '_clus' + str(-1*j)+ '.png'
            f_name = os.path.join(f_dir, custom_name)
            plt.savefig(fname=f_name)
            plt.close()
else:
    embedding, clusters, sizes = mp_op.transform(visualization_level = levels[vis_level],
                                                         cluster_level = levels[-1*clus_level])
    scprep.plot.scatter2d(embedding, s = 100*np.sqrt(sizes), c = clusters,
                                  fontsize=16, ticks=False,
                                  label_prefix="Multiscale PHATE", figsize=(10,8))
    f_dir = os.path.expanduser(wdir)
    custom_name = 'embedding_vis' + str(vis_level) + '_clus' + str(clus_level)+ '.png'
    f_name = os.path.join(f_dir, custom_name)
    plt.savefig(fname=f_name)
    plt.close()


print('Finding expression')
expression = pd.DataFrame()
for i in range(len(genes)):
    expression[genes[i]] = mp_op.get_expression(data[genes[i]].values,
    visualization_level = levels[vis_level])

print('Pickling important files...')
keep = [embedding, clusters, sizes, data]
keep_names = ['embedding', 'clusters', 'sizes', 'data', 'expression']
for i in range(len(keep)):
    pickle_out = open(keep_names[i] + '.pickle', 'wb')
    pickle.dump(keep[i], pickle_out)
    pickle_out.close()

# Pickling into folder is not working due to TypeError, don't know why

print('Exporting embedding and clusters...')
msPHATE_embedding = pd.DataFrame(embedding,
                        index = pd.read_csv(cell_path, header=None).iloc[:, 0].to_list())
msPHATE_embedding.to_csv('msPHATE_embedding.csv', header=False)
msPHATE_clusters = pd.DataFrame(clusters)
msPHATE_clusters.to_csv('msPHATE_clusters.csv', header=False, index=False)

print(data)
print(levels)
print('Done')
