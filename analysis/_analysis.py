"""Analysis functionality."""

from typing import Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

from utilities import nonlinear_transform
from datastructure import SPECase


def weighted_time_integral(data: np.ndarray, spe_case: SPECase):
    """Compute the weighted time integral of a data field."""

    return np.sqrt(
        np.sum(np.diag(spe_case.time_weight) @ np.square(data), axis=0)
        / np.sum(spe_case.time_weight)
    )


def weighted_time_integral_diff(
    data1: np.ndarray, data2: np.ndarray, spe_case: SPECase
):
    """Compute the weighted time integral of a difference of data fields."""

    return np.sqrt(
        np.sum(np.diag(spe_case.time_weight) @ np.square(data1 - data2))
        / np.sum(spe_case.time_weight)
    )


def scalar_data_series_to_distance_matrix(
    data: dict, index_map: dict, key: str, spe_case: SPECase, order: int = 2
) -> np.ndarray:
    """Distance for scalar data integrated over time.

    Apply the SPE11 time weight to the data to distinguish between the early and late
    effects.

    Parameters:
    data (dict): Data dictionary
    index_map (dict): Index map, mapping the keys to the index in the distance matrix
    key (str): Key of the data field to compare
    spe_case (SPECase): SPE case
    order (int): Order of the distance metric

    Returns:
    np.ndarray: Distance matrix

    """

    # Initialize distance matrix for the cross-comparison of all particiants in the data set
    distance_matrix = np.zeros((len(data), len(data)))

    # Fetch the index of the key in the data format
    key_index = spe_case.data_format[key]

    # Cross comparisons
    for key1, result1 in data.items():
        for key2, result2 in data.items():
            if order == 1:
                distance_matrix[index_map[key1], index_map[key2]] = np.sum(
                    spe_case.time_weight_sparse
                    * spe_case.reporting_times_sparse_delta
                    * np.absolute(
                        result1["scalar_data"][:, key_index]
                        - result2["scalar_data"][:, key_index]
                    )
                )
            elif order == 2:
                distance_matrix[index_map[key1], index_map[key2]] = np.sqrt(
                    np.sum(
                        spe_case.time_weight_sparse
                        * spe_case.reporting_times_sparse_delta
                        * np.square(
                            result1["scalar_data"][:, key_index]
                            - result2["scalar_data"][:, key_index]
                        )
                    )
                )
            else:
                raise ValueError("Order not supported.")

    return distance_matrix


def field_data_distance_matrix(
    data,
    snapshots,
    key,
    timing: Literal["early", "late", "all"],
    spe_case: SPECase,
    order: int = 2,
) -> np.ndarray:
    """..."""

    # Initialize distance matrix for the cross-comparison of all particiants in the data set
    distance_matrix = np.zeros((len(data), len(data)))

    # Apply weighted time integration using a rectangular rule
    if timing == "early":
        container = zip(
            spe_case.early_reporting_times_dense,
            spe_case.early_time_weight_dense,
            spe_case.early_reporting_times_dense_delta,
        )
    elif timing == "late":
        container = zip(
            spe_case.late_reporting_times_dense,
            spe_case.late_time_weight_dense,
            spe_case.late_reporting_times_dense_delta,
        )
    elif timing == "full":
        container = zip(
            spe_case.reporting_times_dense,
            spe_case.time_weight_dense,
            spe_case.reporting_times_dense_delta,
        )
    else:
        raise ValueError("Timing not supported.")

    if order not in [1, 2]:
        raise ValueError("Order not supported.")

    for time, weight, delta in container:
        if order == 1:
            distance_matrix += weight * delta * snapshots[time][key]
        elif order == 2:
            distance_matrix += weight * delta * np.square(snapshots[time][key])

    if order == 1:
        ...
    elif order == 2:
        distance_matrix = np.sqrt(distance_matrix)

    return distance_matrix


def mean_matrix(
    data: dict,
    keys: Optional[list] = None,
    mean_type: Literal["arithmetic", "geometric", "ag"] = "ag",
    order: int = 1,
):
    """Compute a global metric for a data field - collecting selected sparse and dense data."""

    if keys is None:
        keys = data.keys()
    if mean_type == "ag":
        return nonlinear_transform(
            sum([nonlinear_transform(data[key], "symlog") for key in keys]) / len(keys),
            "symlog",
            inverse=True,
        )
    # elif mean_type == "arithmetic":
    #    return sum([data[key] for key in keys]) / len(keys)
    # elif mean_type == "geometric":
    #    return np.sqrt(
    #        sum([np.square(data[key]) for key in keys]) / len(keys)
    #    )
    # elif mean_type == "rms":
    #    return nonlinear_transform(
    #        np.sqrt(
    #            sum(
    #                [
    #                    np.square(nonlinear_transform(data[key], transformation_type))
    #                    for key in keys
    #                ]
    #            )
    #            / len(keys)
    #        ),
    #        transformation_type,
    #        inverse=True,
    #    )
    else:
        raise ValueError("Order not supported.")


# ! ---- PLOTTING ----


def plot_heatmap_distance_matrix(
    distance_matrix: np.ndarray,
    labels: list,
    spe_case: SPECase,
    path=None,
):
    assert distance_matrix.shape[0] == distance_matrix.shape[1]
    assert distance_matrix.shape[0] == len(labels)

    # Plot the distance matrix using a heatmap
    fig, ax = plt.subplots()
    # Adjust figure size
    fig.set_size_inches(14, 14)

    def crop_cmap(cmap, min_val=0.0, max_val=1.0):
        cmap = plt.get_cmap(cmap)
        colors = cmap(np.linspace(min_val, max_val, cmap.N))
        return LinearSegmentedColormap.from_list("cropped_cmap", colors)

    cropped_cmap = crop_cmap("binary", 0.0, 0.4)  # Using 20% of the 'viridis' colormap
    im = ax.matshow(distance_matrix, cmap=cropped_cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    # Add label to color bar
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("SPE11 distance")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    # Add labels as x labels, but rotate them
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=90)
    # Remove x-ticks on the bottom
    ax.tick_params(axis="x", bottom=False)
    # Add numbers to the heatmap - for each box add the value in white
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(
                i,
                j,
                f"{distance_matrix[i, j]:.2f}",
                ha="center",
                va="center",
                color="k",
                fontsize=7,
            )

    # Add raster to around the diagonal
    for i in range(len(labels)):
        ax.add_patch(
            plt.Rectangle(
                (i - 0.5, i - 0.5),
                1,
                1,
                fill=False,
                edgecolor="k",
                lw=0.5,
            )
        )

    for orientation, lbls in zip(
        ["x", "y"], [ax.get_xticklabels(), ax.get_yticklabels()]
    ):
        for lbl in lbls:
            result = lbl.get_text()
            if result in labels:
                team = lbl.get_text()
                if team[-1].isdigit():
                    team = team[:-1]
                team_color = spe_case.groups_and_colors[team.lower()]
                team_text_color = spe_case.groups_and_text_colors[team.lower()]
                category = spe_case.results_and_categories[result.lower()]
                category_color = spe_case.categories_and_colors[category.lower()]
                lbl.set_bbox(
                    dict(
                        facecolor=category_color if orientation == "x" else team_color,
                        edgecolor="k",  # team_color,
                        boxstyle="square, pad=0.3",
                    ),
                )
                if orientation == "y":
                    lbl.set_color(team_text_color)

    if path is not None:
        plt.tight_layout()
        plt.savefig(path, dpi=1200)
        print(f"Saved heatmap to {path}.")

    plt.show()
