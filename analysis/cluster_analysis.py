from typing import Literal

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import cdist, squareform

from datastructure import SPECase
from spe11_io import nonlinear_transform


def _cluster_analysis(
    distance_matrix: np.ndarray,
    labels: list,
    linkage_type: Literal[
        "single", "complete", "average", "weighted", "centroid", "median", "ward"
    ] = "average",
):
    assert len(labels) == distance_matrix.shape[0]
    dm = squareform(distance_matrix)
    Z = linkage(dm, method=linkage_type)
    return Z


# ! ---- PLOTTING FUNCTIONS ---- ! #


def plot_distance_matrix_values(distance_matrix: dict[str, np.ndarray]):
    flat_distance_matrix = {
        key: squareform(matrix) for key, matrix in distance_matrix.items()
    }

    plt.figure("values")
    for key, values in flat_distance_matrix.items():
        plt.scatter(values, np.ones_like(values) * (-1), label=key)
    plt.legend()

    plt.figure("raw - srg")
    for i, (key, data) in enumerate(srg_flat_distance_matrix.items()):
        plt.scatter(data, np.ones_like(data) * (-i), label=key)
    plt.legend()
    plt.show()


def plot_subgroup_distance_matrix_values(
    distance_matrix: dict[tuple[str, str], np.ndarray],
):
    flat_distance_matrix = {
        key: squareform(matrix) for key, matrix in distance_matrix.items()
    }

    subgroups = set([key[0] for key in distance_matrix.keys()])
    for subgroup in subgroups:
        subgroup_values = {
            key[1]: values
            for key, values in flat_distance_matrix.items()
            if key[0] == subgroup
        }

        plt.figure(f"{subgroup} values")
        for level, (key, values) in enumerate(subgroup_values.items()):
            plt.scatter(values, np.ones_like(values) * (-level), label=key)
        plt.legend()

    plt.show()


def plot_linkage_clustering(
    distance_matrix: np.ndarray,
    labels: list,
    linkage_type: Literal[
        "single", "complete", "average", "weighted", "centroid", "median", "ward"
    ] = "average",
) -> None:
    Z = _cluster_analysis(distance_matrix, labels, linkage_type)
    plt.figure("dendrogram")
    plt.title(f"Dendrogram - {linkage_type}")
    dendrogram(
        Z,
        labels=labels,
        orientation="left",
        # distance_sort="descending",
        count_sort="ascending",
    )
    plt.show()


def plot_linkage_clustering_for_subgroups(
    subgroups_distance_matrix: dict[tuple[str, str], np.ndarray],
    subgroups_labels: list,
    linkage_type: Literal[
        "single", "complete", "average", "weighted", "centroid", "median", "ward"
    ] = "average",
) -> None:
    for key, matrix in subgroups_distance_matrix.items():
        labels = subgroups_labels[key[0]]
        Z = _cluster_analysis(matrix, labels, linkage_type)
        plt.figure(f"dendrogram - {key}")
        plt.title(f"Dendrogram - {key} - {linkage_type}")
        dendrogram(
            Z,
            labels=labels,
            orientation="left",
            # distance_sort="descending",
            count_sort="ascending",
        )

    plt.show()


def plot_linkage_clustering_with_colors(
    spe_case: SPECase,
    distance_matrix: np.ndarray,
    labels: list,
    linkage_type: Literal[
        "single", "complete", "average", "weighted", "centroid", "median", "ward"
    ] = "average",
    nonlinear_transformation: Literal["linear", "symlog"] = "linear",
    path: str = None,
):
    plt.figure("dendrogram with colors")
    ax = plt.gca()

    Z = _cluster_analysis(distance_matrix, labels, linkage_type)

    def color_func(k):
        return "black"

    dendrogram(
        Z,
        ax=ax,
        labels=labels,
        orientation="left",
        count_sort=True,
        distance_sort=True,
        link_color_func=lambda k: color_func(k),
    )

    if nonlinear_transformation == "symlog":
        Z[:, 2] = nonlinear_transform(Z[:, 2], "symlog")

    xlbls = ax.get_ymajorticklabels()
    for lbl in xlbls:
        result = lbl.get_text()
        if result in labels:
            team = lbl.get_text()
            if team[-1].isdigit():
                team = team[:-1]
            team_color = spe_case.groups_and_colors[team.lower()]
            category = spe_case.results_and_categories[result.lower()]
            category_color = spe_case.categories_and_colors[category.lower()]
            # Set two boxes, one surrounding the text with pad just to the left, and the other with a smaller pad, but surrounding all
            lbl.set_bbox(
                dict(
                    facecolor=category_color,
                    edgecolor="k",  # team_color,
                    boxstyle="square, pad=0.3",
                ),
            )
            # Add empty text at the same location as the lbl
            txt = ax.text(
                # lbl.get_position()[0],
                0.0,
                lbl.get_position()[1],
                ".",
                color=team_color,
            )
            # Set bbox for the text
            txt.set_bbox(
                dict(
                    facecolor=team_color,
                    edgecolor=team_color,
                    boxstyle="square,pad=0.0",
                ),
            )

    plt.title(
        f"{spe_case.variant.upper()} - Hierarchical clustering of submissions based on SPE11 distance"
    )
    if nonlinear_transformation == "symlog":
        plt.xlabel("Symlog(SPE11 distance)")
    else:
        plt.xlabel("SPE11 distance")

    # Set the size of the plot exactly to be the size of the dendrogram plus the labels
    plt.gcf().set_size_inches(9, 8)
    if path:
        plt.tight_layout()
        plt.savefig(path, dpi=1000)
    plt.show()


# ! ---- MEDIAN CLUSTERS ---- ! #


def determine_median_cluster(
    distance_matrix: np.ndarray,
    labels: list,
    linkage_type: Literal[
        "single", "complete", "average", "weighted", "centroid", "median", "ward"
    ] = "average",
    mean_type: Literal["arithmetic", "geometric", "ag"] = "ag",
):
    Z = _cluster_analysis(distance_matrix, labels, linkage_type)

    # Form (a single) flat cluster
    clusters = fcluster(Z, t=1, criterion="maxclust")

    # Find the centroid - the one with smallest sum of distances to all other points
    for cluster_id in np.unique(clusters):
        cluster_indices = np.where(clusters == cluster_id)[0]
        cluster_distances = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
        means = []
        for i in range(len(cluster_indices)):
            distances = cluster_distances[i, :].tolist()
            distances.pop(i)
            distances = np.array(distances)
            means += [mean(distances, mean_type)]
        centroid_index = cluster_indices[np.argmin(means)]
        print(
            f"Cluster {cluster_id} centroid index: {centroid_index}, label: {labels[centroid_index]}"
        )

    # Find the partner of the identified cluster
    if centroid_index in Z[:, 0].astype(int).tolist():
        row = np.where(Z[:, 0].astype(int) == centroid_index)[0][0]
        partner_index = Z[row, 1].astype(int)
    elif centroid_index in Z[:, 1].astype(int).tolist():
        row = np.where(Z[:, 1].astype(int) == centroid_index)[0][0]
        partner_index = Z[row, 0].astype(int)
    else:
        raise ValueError("Centroid index not found in linkage matrix.")

    if partner_index < len(labels):
        median_cluster = [labels[centroid_index], labels[partner_index]]
    else:
        print("Partner index not found in labels.")
        median_cluster = [labels[centroid_index]]

    return median_cluster


def mean(
    values: np.ndarray,
    mean_type: Literal["arithmetic", "geometric", "ag"] = "ag",
) -> float:
    if mean_type == "arithmetic":
        transformation_type = "linear"
    elif mean_type == "geometric":
        transformation_type = "log"
    elif mean_type == "ag":
        transformation_type = "symlog"
    else:
        raise ValueError("Mean type must be arithmetic, geometric, or ag.")

    return nonlinear_transform(
        np.mean(nonlinear_transform(values, transformation_type)),
        transformation_type,
        inverse=True,
    )


def std(
    values: np.ndarray,
    mean_type: Literal["arithmetic", "geometric", "ag"] = "ag",
) -> float:
    if mean_type == "arithmetic":
        transformation_type = "linear"
    elif mean_type == "geometric":
        transformation_type = "log"
    elif mean_type == "ag":
        transformation_type = "symlog"
    else:
        raise ValueError("Mean type must be arithmetic, geometric, or ag.")

    transformed_values = nonlinear_transform(values, transformation_type)
    transformed_mean = np.mean(transformed_values)
    transformed_std = np.std(transformed_values)
    mean = nonlinear_transform(
        transformed_mean,
        transformation_type,
        inverse=True,
    )
    std = max(
        nonlinear_transform(
            transformed_mean + transformed_std,
            transformation_type,
            inverse=True,
        )
        - mean,
        mean
        - nonlinear_transform(
            transformed_mean - transformed_std,
            transformation_type,
            inverse=True,
        ),
    )
    return std


def mean_distance_to_group(
    distance_matrix: np.ndarray,
    probe: str,
    group: list,
    mean_type: Literal["arithmetic", "geometric", "ag"] = "ag",
):
    # Determine the mean distance to the group
    index = group.index(probe)
    # Extract the 0-distance to the probe itself
    distances = distance_matrix[index, :].tolist()
    distances.pop(index)
    distances = np.array(distances)
    group_distance = mean(distances, mean_type)

    return group_distance


def std_distance_to_group(
    distance_matrix: np.ndarray,
    probe: str,
    group: list,
    mean_type: Literal["arithmetic", "geometric", "ag"] = "ag",
):
    # Determine the mean distance to the group
    index = group.index(probe)
    # Extract the 0-distance to the probe itself
    distances = distance_matrix[index, :].tolist()
    distances.pop(index)
    distances = np.array(distances)
    standard_deviation = std(distances, mean_type)

    return standard_deviation


def argmin_distance(distances: dict[str, float]) -> str:
    return min(distances, key=distances.get)


def centroid_analysis(
    distance_matrix: np.ndarray,
    labels: list,
    linkage_type: Literal[
        "single", "complete", "average", "weighted", "ward"
    ] = "average",
    mean_type: Literal["arithmetic", "geometric", "ag"] = "ag",
):
    # Perform hierarchical clustering
    Z = _cluster_analysis(distance_matrix, labels, linkage_type)

    # Form flat clusters
    clusters = fcluster(Z, t=1, criterion="maxclust")

    # Find the centroid of each cluster
    for cluster_id in np.unique(clusters):
        cluster_indices = np.where(clusters == cluster_id)[0]
        cluster_distances = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
        means = []
        for i in range(len(cluster_indices)):
            distances = cluster_distances[i, :].tolist()
            distances.pop(i)
            distances = np.array(distances)
            means += [mean(distances, mean_type)]
        centroid_index = cluster_indices[np.argmin(means)]
        print(
            f"Cluster {cluster_id} centroid index: {centroid_index}, label: {labels[centroid_index]}"
        )
