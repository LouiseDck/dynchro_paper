# Aligning two patients with a batch effect

This example demonstrates how to align two patients with a batch effect.
We use generated data to simulate two patients, using `dyngen`. Afterwards, we align the patients and consider the differences between a classical analysis and one that takes the batch effect into account, up to the trajectory inference step.

`1_generate_data.R` uses `dyngen` to generate the datasets.
`2_preprocess/py` preprocesses the data of the two patients.
`3_analysis.py` performs the analysis and visualises the results.
