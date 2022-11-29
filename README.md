# msPHATE-Seurat
## Key variables

### R_prework.R
data_file = '~/Jansky/data/adrenal_medulla_Seurat.RDS.gz'
wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat'
use_scaled = TRUE

### msphate.py
Can input raw counts mtx or scaled data into python

wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat'
min_reads = 200
mincells = 3
npca = 20
gran = 0.1
run_multiple_embeddings=False
generate_tree=True
multiple_spread = 1
vis_level = 0
clus_level = 3
marker_dict = {'SCPs': ['SOX10', 'PLP1', 'ERBB3'],
               'Neuroblasts': ['ISL1', 'ALK'],
               'Chromaffin cells': ['DBH', 'TH', 'CHGA']}

### R_postwork.R
wdir = '~/Jansky/ms_phate/msphate_jansky/msPHATE-Seurat'

## Usage
`DimPlot(medulla_seurat, reduction = 'msphate', group.by = 'msphate_clusters')`
Please ignore the '"" not responding' when running Python script!