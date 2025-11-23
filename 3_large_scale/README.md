# dynchro clusters similar and dissimilar samples based on trajectories

This usecase shows how dynchro clusters datasets based on the trajectories contained therein.

## Workflow

1. `0_functions.py`: utility functions
2. `1_split_data.ipynb`: Analyse data, split in replicates and lineages.
3. `2_hvgs.ipynb`: Calculate HVGs per replicate, and find the union of them.
4. `3_subsampling.ipynb`: Subsample datasets, to generate 30 datasets.
5. `4_dynchro_ss.ipynb`: Calculate dynchro distances between each dataset and cluster.

# dynchro identifies and tracks similarity and dissimilarity between trajectory pseudotimes

1. `5_dynchro_examples.ipynb` showcases a good and bad alignment example.
2. `0_figures.ipynb` reproduces the paper figures exactly.
