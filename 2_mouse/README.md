# dynchro identifies correspondences between branching trajectories

This usecase demonstrates how dynchro identifies matching trajectories, using a mouse hematopoeisis dataset with 3 knockouts, see: [Olsson, A. et al. Single-cell analysis of mixed-lineage states leading to a binary cell fate choice. Nature (2016)](https://doi.org/10.1038/nature19348)

This data was used by trajectory inference tools, such as [STREAM](https://doi.org/10.1038/s41467-019-09670-4).

## Workflow

1. `utils.py` contains functions lifted from [the STREAM TI tool](https://github.com/pinellolab/STREAM), in order to be able to run the analysis performed by the now failing tool.
2. `1_prep_data.ipynb` details the STREAM preprocessing the data undergoes
3. `2_dynchro.ipynb` details the dynchro analysis
4. `3_figures.ipynb` reproduces the paper figures exactly.
