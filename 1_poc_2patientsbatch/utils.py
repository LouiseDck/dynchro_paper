import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def plot_iroot(data):
    colormap_iroot = np.array(["#C0C0C0"] * data.n_obs)
    colormap_iroot[data.uns["iroot"]] = "#FF0000"

    plt.scatter(
        data.obsm["X_scanorama"][:, 0],
        data.obsm["X_scanorama"][:, 1],
        c=colormap_iroot,
        s=10,
        alpha=0.1,
    )


def comp_kde(dataset, milestones, batches = None, milestone_column = "milestones", batch_column = "batch", time_column = "dpt_pseudotime", colors = None, title = None):
    """
    Function to plot the kernel density estimation of the data
    """
    if batches is None:
        batches = dataset.obs[batch_column].cat.categories.values
    if isinstance(milestones, str) or isinstance(milestones, int):
        milestones = [milestones] * len(batches)
    if colors is None:
        colors = [None] * len(batches)

    for batch, milestone, color in zip(batches, milestones, colors):
        dataset_batch = dataset[dataset.obs[batch_column] == batch]
        dataset_batch = dataset_batch[dataset_batch.obs[milestone_column] == milestone]

        print(dataset_batch)

        sns.kdeplot(dataset_batch.obs[time_column], label=f"{batch} {milestone}", color=color, fill = True, alpha=0.5)

    plt.legend()
    plt.figure()


def kdeplot2(adata, batch, milestone, milestone_column = "milestones", time = "dpt_pseudotime", color=None):
    """
    Plot the kernel density estimate of the pseudotime for a given batch and milestone.
    """
    adata = adata[adata.obs[milestone_column] == milestone]
    adata = adata[adata.obs["batch"] == batch]
    sns.kdeplot(adata.obs[time], color=color, label=f"{batch} {milestone}", fill=True, alpha=0.5)

def kdeplot(adata, batch, milestone, time = "dpt_pseudotime", color=None):
    adata_mst = adata[adata.obs["milestones"] == milestone]
    adata_btc = adata_mst[adata_mst.obs["batch"] == batch]

    sns.kdeplot(adata_btc.obs[time], label=f"{batch} {milestone}", fill=True, alpha=0.5, color= color)