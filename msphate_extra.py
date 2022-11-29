import multiscale_phate as mp
import numpy as np
import pandas as pd
import scprep
import os
import matplotlib.pyplot as plt
import pickle

wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat'
min_reads = 200
mincells = 3
npca=20
gran = 0.1
vis_level = 0
clus_level = -3
marker_dict = {'SCPs': ['SOX10', 'PLP1', 'ERBB3'],
               'Neuroblasts': ['ISL1', 'ALK'],
               'Chromaffin cells': ['DBH', 'TH', 'CHGA']}
names = list(marker_dict.keys())
genes = list(marker_dict.values())


print('Unpickling important files...')
kept_names = ['embedding', 'clusters', 'sizes', 'data', 'expression']
pickle_path = os.path.join(wdir, 'pickle')

custom_name = kept_names[0]+'.pickle'
f_name = os.path.join(pickle_path, custom_name)
pickle_in = open(f_name, 'rb')
embedding = pickle.load(pickle_in)
pickle_in.close()

custom_name = kept_names[1]+'.pickle'
f_name = os.path.join(pickle_path, custom_name)
pickle_in = open(f_name, 'rb')
clusters = pickle.load(pickle_in)
pickle_in.close()

custom_name = kept_names[2]+'.pickle'
f_name = os.path.join(pickle_path, custom_name)
pickle_in = open(f_name, 'rb')
sizes = pickle.load(pickle_in)
pickle_in.close()

custom_name = kept_names[3]+'.pickle'
f_name = os.path.join(pickle_path, custom_name)
pickle_in = open(f_name, 'rb')
data = pickle.load(pickle_in)
pickle_in.close()

custom_name = kept_names[4]+'.pickle'
f_name = os.path.join(pickle_path, custom_name)
pickle_in = open(f_name, 'rb')
data = pickle.load(pickle_in)
pickle_in.close()




##print('Generating marker expression plot for msPHATE...')
##row_list=[]
##for i in marker_dict:
##        row.list.append(len(marker_dict[i]))
##for i in range(0, len(row_list)-1):
##        for j in marker_dict:
##                if len(marker_dict[j])==i:
##                        cols = i
##                        break
##        else:
##                continue
##        break
##fig, axes = plt.subplots(nrows=len(list(marker_dict.keys())), ncols=cols,
##                         figsize=(14, 4))
##for i, ax in enumerate(axes.flatten()):
##        for j in marker_dict:
##                scprep.plot.scatter2d(embedding, title=marker_dict[j]+list(marker_dict.values())[i],
##                                      s = 25*np.sqrt(sizes),
##                                      c=expression[list(marker_dict.values())[i]],
##                                      legend_anchor=(1,1), ax=ax,
##                                      xticks=False, yticks=False,
##                                      label_prefix="PHATE", fontsize=16, cmap = 'RdBu_r')
##fig.tight_layout()
##plt.savefig(fname='clusters_markers.png')
##plt.close()
# Assuming it goes through the grid of subplots in simple order and skips empty ones
# Have not tested yet

##print('Generating marker expression plot by cluster...')
##scprep.plot.marker_plot(expression, clusters, markers = marker_dict,
##                        gene_names = list(marker_dict.values()), cmap = 'RdBu_r')
##plt.savefig(fname='markers.png')
##plt.close()

print(levels)
print(expression)
