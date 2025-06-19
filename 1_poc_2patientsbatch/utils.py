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

def plot_datasets(dataset, milestones = None, batches = None, milestone_column = "milestones", time_column = "dpt_pseudotime"):
    milestones = ["sB", "sBmid", "sC", "sC_batch", "sEndC"] if milestones is None else milestones
    batches = ["control", "premature stop"] if batches is None else batches


    comp_kde(dataset, milestones[0], time_column = "slingPseudotime_1", colors = ["#ecf39e", "#caf0f8"])
    comp_kde(dataset, milestones[1],time_column = "slingPseudotime_1", colors = ["#90a955", "#90e0ef"])
    
    kdeplot2(dataset, batches[0], milestones[2],time = "slingPseudotime_1", milestone_column = "milestones", color="#4f772d")
    kdeplot2(dataset, batches[1], milestones[3],time = "slingPseudotime_1", milestone_column = "milestones", color="#00b4d8")
    kdeplot2(dataset, batches[0], milestones[4], time = "slingPseudotime_1",milestone_column = "milestones", color="#31572c")
    plt.legend()

def compare(dataset, milestones, batch = None, time1 = "sling_Pseudotime_1", time2 = "sim_time", colors = None):
    if batch is not None:
        dataset = dataset[dataset.obs["batch"] == batch]
    
    if isinstance(milestones, str) or isinstance(milestones, int):
        milestones = [milestones]

    dataset = dataset[dataset.obs["milestones"].isin(milestones)]
    
    sns.kdeplot(dataset.obs[time1], label=f"{batch} {milestones} {time1}", fill=True, alpha=0.5, color=colors[0])
    sns.kdeplot(dataset.obs[time2], label=f"{batch} {milestones} {time2}", fill=True, alpha=0.5, color=colors[1])
    plt.legend()

def norm(x, max_value=1):
    """
    Normalize the input array to the range [0, max_value].
    """
    return (x - np.min(x)) / (np.max(x) - np.min(x)) * max_value

def get_kde_eval(vector, bandwith = 0.1, x = 100):
    """
    Get the kernel density estimation of a vector.
    """
    from scipy.stats import gaussian_kde

    kde = gaussian_kde(vector, bw_method=bandwith)
    x_values = np.linspace(0, 1, x)
    y_values = kde(x_values)

    return x_values, y_values, kde

