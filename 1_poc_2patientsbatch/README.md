# dynchro shows the effect of batch removal on trajectories and identifies premature stops

This usecase demonstrates how to align two patients with a batch effect.
We use generated data to simulate two patients, using `dyngen`. Afterwards, we consider batch effect correction methods, and how dynchro can pinpoint issues with the batch effect correction methods.

## Workflow

1. `1_generate_data.R` uses `dyngen` and `0_diverging_kinetics.R` to generate the datasets.
2. `slingshot.R` runs slingshot, at different points in the process. Run the corresponding part of this script during preprocessing, or specific analysis.
3. `2_preprocess.ipynb` preprocesses the two synthetic datasets.
4. `3x_method` runs the corresponding method on the datasets. For fastmnn, run `3c_fastmnn.ipynb` first and switch to `3c_fastmnn.R` when asked to.
5. `4_figures.ipynb` reproduces the paper figures exactly.
