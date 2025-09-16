import anndata as ad
import numpy as np
import pandas as pd

import dynchro


def write_with_pc(adata, name):
    for pckey in adata.uns["pseudocells"]:
        if isinstance(adata.uns[f"{pckey}_pseudotime"], pd.core.series.Series):
            pseudocells = adata.uns[f"{pckey}_pseudotime"]
            pseudocells_index = pseudocells.index.values
            pseudocells_values = pseudocells.values

            adata.uns[f"{pckey}_pseudotime_index"] = pseudocells_index
            adata.uns[f"{pckey}_pseudotime_values"] = pseudocells_values

            adata.uns[f"{pckey}_pseudotime"] = "removed"

    adata.write_h5ad(f"{name}.h5ad")


def read_with_pc(file):
    adata = ad.read_h5ad(file)

    for pckey in adata.uns["pseudocells"]:
        if f"{pckey}_pseudotime" in adata.uns and adata.uns[f"{pckey}_pseudotime"] == "removed":
            pseudocells_index = adata.uns[f"{pckey}_pseudotime_index"]
            pseudocells_values = adata.uns[f"{pckey}_pseudotime_values"]

            pseudocells = pd.Series(pseudocells_values, index=pseudocells_index)
            adata.uns[f"{pckey}_pseudotime"] = pseudocells

    return adata


def get_config(trajectories, sort=True):
    config = {
        "compare_lineages_pseudocells": False,
        "compare_lineages_pseudocells_label": "pseudocells_50",
        "compare_trajectories_pseudocells": False,
        "compare_trajectories_pseudocells_label": "pseudocells_50",
        "align_pseudocells": False,
    }

    common_vars = list(set.intersection(*(set(trajectory.var_names) for trajectory in trajectories)))

    lineages = [
        [
            dynchro.tl.get_counts_common_vars(
                trajectory, config, "compare_trajectories_pseudocells", linlabel, common_vars, get_x=False
            )
            for linlabel in trajectory.uns["lineage_labels"]
        ]
        for trajectory in trajectories
    ]

    if sort:
        lineages = [[lin[np.argsort(lin.obs.pseudotime)] for lin in lineages1] for lineages1 in lineages]

    return lineages


import matplotlib.pyplot as plt
import seaborn as sns

def flatten(values : list) -> list:
    """
    Flatten a list of lists into a single list.
    """
    flat_list = []

    for sublist in values:
        if isinstance(sublist, list):
            # If the item is a list, extend the result with its contents
            flat_list.extend(sublist)
        else:
            # If the item is not a list, append it directly if it is not None
            if sublist is not None:
                flat_list.append(sublist)

    return flat_list


def plot_cost_distances(cost, distances, path1, path2):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(12, 4)
    sns.heatmap(cost, ax=ax1)
    distances = distances[1:, 1:]
    sns.heatmap(distances, ax=ax2)
    path1 = flatten(path1.values())
    path2 = flatten(path2.values())
    print(len(path1), len(path2))
    plt.scatter(path2, path1, c="white", s=3)
    plt.show()


class DTWResult:
    def __init__(self, lin1, lin2, cost, distances, path1, path2, dist, dist_norm, dist_norm_path):
        self._lin1 = lin1
        self._lin2 = lin2
        self._cost = cost
        self._distances = distances
        self._path1 = path1
        self._path2 = path2
        self._dist = dist
        self._dist_norm = dist_norm
        self._dist_norm_path = dist_norm_path

    def cost(self):
        return self._cost

    def distances(self):
        return self._distances

    def path1(self):
        return self._path1

    def path2(self):
        return self._path2

    def dist(self):
        return self._dist

    def dist_norm(self):
        return self._dist_norm

    def dist_norm_path(self):
        return self._dist_norm_path

    def get_warping(self):
        return self._path1, self._path2

    def plot_cost_distances(self):
        plot_cost_distances(self._cost, self._distances, self._path1, self._path2)


def dynchro_wrapper(ad1, ad2, distance_measure="correlation"):
    trajectories = [ad1, ad2]
    # lineages = get_config(trajectories)

    common_vars = list(set(ad1.var_names).intersection(set(ad2.var_names)))

    lineages1 = [
        ad1[ad1.obs[label]][:, common_vars]
        for label in ad1.uns["lineage_labels"]
    ]

    lineages2 = [
        ad2[ad2.obs[label]][:, common_vars]
        for label in ad2.uns["lineage_labels"]
    ]

    # sort according to pseudotime_key
    lineages1 = [lin[lin.obs["pseudotime"].argsort()] for lin in lineages1]
    lineages2 = [lin[lin.obs["pseudotime"].argsort()] for lin in lineages2]

    results = []

    for lineage1, label1 in zip(lineages1, ad1.uns["lineage_labels"]):
        for lineage2, label2 in zip(lineages2, ad2.uns["lineage_labels"]):
            # print(lineage2.X.shape, lineage1.X.shape)
            total_dist, cost, distances_matrix = dynchro.tl.dtw(
                lineage1, lineage2, distance="euclidean", mode="only_results"
            )
            path1, path2 = dynchro.tl.traceback(D = distances_matrix)
            norm_distance = total_dist / (cost.shape[0] * cost.shape[1])

            # lin1, lin2 = lineage1.X, lineage2.X
            # total_dist, cost, distances = dynchro.tl.dtw(lin1, lin2, distance=distance_measure)
            # path1, path2 = dynchro.tl.traceback(distances)

            total_dist_norm = total_dist / (lineage1.X.shape[0] + lineage2.X.shape[0])
            # total_dist_norm2 = total_dist / len(path1)

            res = DTWResult(
                label1, label2, cost, distances_matrix, path1, path2, total_dist, total_dist_norm, norm_distance
            )
            results.append(res)

            # plot_cost_distances(cost, distances, path1, path2)

    return results
