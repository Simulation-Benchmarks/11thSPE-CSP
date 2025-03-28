"""Collected I/O functionality."""

import logging
from pathlib import Path
from typing import Literal, Optional
from warnings import warn

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform

from datastructure import SPECase

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

logger = logging.getLogger(__name__)

# ! ---- NAME MANAGEMENT ----


def split_result_name(result_name: str):
    """Split the result name into its components."""

    # Check if the last component is a number - corresponds to the result number
    if result_name[-1].isdigit():
        team = result_name[:-1]
        result = result_name[-1]
    else:
        team = result_name
        result = None
    return team, result


def identify_data(
    data_folder: Path,
    spe_case: SPECase,
    participants,
):
    """Identify participants based on the folder structure and existence of files.

    The structure of the folder is assumed to be:
        - team1
            - {speID}_spatial_map_<time>.csv
            - speID_time_series.csv
        ...
        - team<N>
            - speID_time_series.csv
        - ...

    Multiple submissions by the same team are marked with a number at the end of the team name.
    If a team only has one submission, the number is anyhow used.

    """
    # Collect the data
    team_submissions = {}
    for submission in data_folder.iterdir():
        participants[submission.name] = {
            "sparse": submission / f"{spe_case.variant}_time_series.csv",
        }
        team, result = split_result_name(submission.name)
        if team not in team_submissions:
            team_submissions[team] = []
        team_submissions[team].append(result)

    # Clean the keys - omit number if the only submission by a team is "team1"
    for team, submissions in team_submissions.items():
        if submissions == ["1"]:
            participants[team] = participants.pop(f"{team}1")

    return participants


# ! ---- DATA MANAGEMENT ----


def read_time_series_data(path, spe_case: SPECase):
    """Read raw data from a CSV file."""

    # Check if any header is non-numeric
    df = pd.read_csv(path)
    non_numeric_headers = [
        str(col).lower()
        for col in df.columns
        if not str(col).replace(".", "", 1).isdigit()
    ]
    if "nan" in non_numeric_headers:
        non_numeric_headers.remove("nan")
    if "inf" in non_numeric_headers:
        non_numeric_headers.remove("inf")

    # Read the data (either without header or with)
    if non_numeric_headers:
        # print("Non-numeric headers found:", non_numeric_headers)
        df = pd.read_csv(path, skiprows=1, usecols=range(len(spe_case.data_format)))
    else:
        # print(f"All headers are numeric for {path}.")
        df = pd.read_csv(path, header=None, usecols=range(len(spe_case.data_format)))

    data = df.to_numpy().astype(float)

    # Make sure the data format matches
    assert data.shape[1] == len(spe_case.data_format), (
        f"Data format does not match {data.shape[1]} vs {len(spe_case.data_format)}."
    )

    return data


def read_extra_data(path, data_types):
    extra_data = {}

    df = pd.read_csv(path)
    header = df.columns.to_list()
    header = [x.strip() for x in header]

    # Strip all whitespaces from all cells
    data = df.to_numpy().astype(float)

    for i, key in enumerate(header):
        if "time" in key:
            extra_data["t"] = data[:, i]
        elif key in data_types:
            extra_data[key] = data[:, i]

    return extra_data


def read_field_data_distance_matrix_snapshots(folder, participants, spe_case: SPECase):
    distance_matrix_snapshots = {
        time: {
            "pressure_l2": read_dense_distance_data(
                folder
                / f"{spe_case.variant}_pressure_l2_diff_{time}{spe_case.reporting_time_unit}.csv",
                participants,
            ),
            "pressure_l2s": read_dense_distance_data(
                folder
                / f"{spe_case.variant}_pressure_l2semi_diff_{time}{spe_case.reporting_time_unit}.csv",
                participants,
            ),
            "mass_w1": read_dense_distance_data(
                folder
                / f"{spe_case.variant}_co2mass_w1_diff_{time}{spe_case.reporting_time_unit}.csv",
                participants,
            ),
        }
        for time in spe_case.reporting_times_dense
    }
    if spe_case.non_isothermal:
        for time in spe_case.reporting_times_dense:
            distance_matrix_snapshots[time].update(
                {
                    "temperature_l2s": read_dense_distance_data(
                        folder
                        / f"{spe_case.variant}_temperature_l2semi_diff_{time}{spe_case.reporting_time_unit}.csv",
                        participants,
                    )
                }
            )

    return distance_matrix_snapshots


def read_dense_distance_data(path, participants):
    """Read dense distance data from a CSV file.

    The CSV file contains a header and a first column with group labels.

    The distance matrix is then extracted from the remaining data.

    """
    df = pd.read_csv(path, skiprows=1, header=None)

    # Make sure the groups are identical to the participants
    groups = df.iloc[:, 0].values.tolist()
    # groups = clean_names(groups)
    if not groups == list(participants.keys()):
        warn(f"{groups} vs {list(participants.keys())}")

    # Extract data
    distance_data = df.iloc[:, 1:].values
    return distance_data


# def clean_data(data):
#    # Clean up the data:
#    # - Replace close-to zero values with zero
#    # - Replace negative values with NaN
#    # - Replace +/- inf values with NaN
#
#    clean_data = data.copy()
#    assert not np.any(np.isclose(data, -1)), "-1 values detected"
#    if np.any(data == -np.inf):
#        warn("-inf values detected")
#    clean_data[data == -np.inf] = np.nan
#    if np.any(data == np.inf):
#        warn("+inf values detected")
#    clean_data[data == np.inf] = np.nan
#
#    return clean_data

# ! ---- DATA SANITY CHECKS ----


def check_sanity(data, key: str):
    if np.any(np.isclose(data, -1)):
        logging.info(f"Missing values detected for submission {key}.")
        return False
    if np.any(np.isclose(data, -999)):
        logging.info(f"Missing values detected for submission {key}.")
        return False
    if np.isnan(data).any():
        logging.info(f"NaN values detected for submission {key}.")
        return False
    if np.any(data < 0):
        logging.info(
            f"Negative values (min value {np.min(data)}) detected for submission {key}."
        )
    if np.any(data == np.inf):
        logging.info(f"Positive infinity values detected for submission {key}.")
        return False
    if np.any(data == -np.inf):
        logging.info(f"Negative infinity values detected for submission {key}.")
        return False
    return True


# ! ---- DATA MODIFICATION ----


def replace_mC_values(data, extra_time, extra_data, spe_case: SPECase, key):
    """Replace M values in the data with the corresponding values from the extra data."""

    logger.info(f"Replace M-values for submission {key}.")
    clean_data = data.copy()
    interp_extra_data = np.interp(
        data[:, spe_case.data_format["t"]], extra_time, extra_data
    )

    clean_data[:, spe_case.data_format["M_C"]] = interp_extra_data
    return clean_data


def replace_sparse_values(
    data, extra_time, extra_data, sparse_id, spe_case: SPECase, key
):
    """Replace M values in the data with the corresponding values from the extra data."""

    logger.info(f"Replace {sparse_id} for submission {key}.")
    clean_data = data.copy()
    interp_extra_data = np.interp(
        data[:, spe_case.data_format["t"]], extra_time, extra_data
    )

    clean_data[:, spe_case.data_format[sparse_id]] = interp_extra_data
    return clean_data


def set_zero_boundaryCO2_values(data, spe_case: SPECase):
    """Set boundaryCO2 values to zero."""

    clean_data = data.copy()
    clean_data[:, spe_case.data_format["boundaryCO2"]] = 0.0
    return clean_data


def interpolate_data_reporting_times(data, spe_case: SPECase):
    """Interpolate data (in table format) to reporting times."""

    # Allocate memory for the interpolated data
    interp_data = np.zeros(
        (len(spe_case.reporting_times_sparse), len(spe_case.data_format))
    )

    # Define the time column
    assert spe_case.data_format["t"] == 0, "Time must be the first column."
    interp_data[:, 0] = spe_case.reporting_times_sparse

    # Interpolate the remaining columns
    for i in range(1, len(spe_case.data_format)):
        interp_data[:, i] = np.interp(
            spe_case.reporting_times_sparse,
            data[:, spe_case.data_format["t"]],
            data[:, i],
        )

    return interp_data


def reduce_to_increase_over_time(data):
    """Mostly unchanged for most data points (besides pressures)."""
    reduced_data = data.copy()
    reduced_data -= data[0]
    return reduced_data


# ! ---- DATA INTEGRATION ----


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


# ! ---- DATA TRANSFORMATION ----


def rescale_clean_sparse_data(
    participants,
    scaling_type: Literal["max", "median", "non_trivial_median"],
    spe_case: SPECase,
):
    """Rescale the sparse data for each participant under "rescaled_data".

    Need to take care of the type of the data. In particular. Phase data
    (immA, mobA, dissA, sealA) as well as (immB, mobB, dissB, sealB) are
    linearly scaled should be seen in combination.

    """
    # Collect the integrals of all participants
    collected_integrals = []
    for _, participant in participants.items():
        collected_integrals.append(participant["weighted_time_integral"])
    collected_integrals = np.array(collected_integrals)

    # Replace (mobA, immA, dissA) columns and (mobB, immB, dissB) rows with their respective sums
    indices_A = [spe_case.data_format[key] for key in ["mobA", "immA", "dissA"]]
    indices_B = [spe_case.data_format[key] for key in ["mobB", "immB", "dissB"]]
    sum_A = np.sum(collected_integrals[:, indices_A], axis=1)
    sum_B = np.sum(collected_integrals[:, indices_B], axis=1)
    collected_integrals[:, indices_A] = sum_A[:, None]
    collected_integrals[:, indices_B] = sum_B[:, None]

    # Determine max values and median values
    num_columns = len(spe_case.data_format)
    max_values = np.max(collected_integrals, axis=0)
    median_values = np.median(collected_integrals, axis=0)
    tol = 1e-10 * max_values
    non_trivial_median_values = np.array(
        [
            np.median(collected_integrals[:, i][collected_integrals[:, i] > tol[i]])
            for i in range(num_columns)
        ]
    )

    # Rescale data
    assert scaling_type in ["max", "median", "non_trivial_median"]
    if scaling_type == "max":
        scaling_values = np.divide(1.0, max_values)
    elif scaling_type == "median":
        scaling_values = np.divide(1.0, median_values)
    elif scaling_type == "non_trivial_median":
        scaling_values = np.divide(1.0, non_trivial_median_values)

    for _, participant in participants.items():
        participant["rescaled_data"] = participant["data"] * scaling_values


def determine_reference_value_distance_matrix(
    distance_matrix: np.ndarray,
    scaling_type: Literal["max", "median", "non_trivial_median"],
):
    """Rescale the distance matrices based on the upper triangle values.

    Assume that the distance matrix is symmetric and that the diagonal is zero.
    And that all values are non-negative.

    """

    # Collect all values of the upper triangle as flat vector
    values = squareform(distance_matrix)

    # Determine max values and median values
    max_value = np.max(values)
    median_value = np.median(values)

    # Consider the median among non-trivial values
    assert np.all(values >= -1e-20)
    tol = 1e-10 * np.max(values)
    if np.count_nonzero(values > tol) > 0:
        non_trivial_median_value = np.median(values[values > tol])
    else:
        non_trivial_median_value = 1.0

    # Rescale data
    assert scaling_type in ["max", "median", "non_trivial_median"]
    if scaling_type == "max":
        scaling_value = max_value
    elif scaling_type == "median":
        scaling_value = median_value
    elif scaling_type == "non_trivial_median":
        scaling_value = non_trivial_median_value

    if scaling_value is np.nan:
        scaling_value = 1.0

    return scaling_value


def rescale_distance_matrix(
    distance_matrix: np.ndarray,
    scaling_type: Literal["max", "median", "non_trivial_median"],
):
    """Rescale the distance matrices based on the upper triangle values.

    Assume that the distance matrix is symmetric and that the diagonal is zero.
    And that all values are non-negative.

    """

    # Collect all values of the upper triangle as flat vector
    values = squareform(distance_matrix)

    # Determine max values and median values
    max_value = np.max(values)
    median_value = np.median(values)

    # Consider the median among non-trivial values
    assert np.all(values >= -1e-20)
    tol = 1e-10 * np.max(values)
    if np.count_nonzero(values > tol) > 0:
        non_trivial_median_value = np.median(values[values > tol])
    else:
        non_trivial_median_value = 1.0

    # Rescale data
    assert scaling_type in ["max", "median", "non_trivial_median"]
    if scaling_type == "max":
        scaling_value = 1.0 / max_value
    elif scaling_type == "median":
        scaling_value = 1.0 / median_value
    elif scaling_type == "non_trivial_median":
        scaling_value = 1.0 / non_trivial_median_value

    if scaling_value is np.nan:
        scaling_value = 1.0

    return distance_matrix * scaling_value


def rescale_dense_distance_matrices(
    distance_matrix,
    scaling_type: Literal["max", "median", "non_trivial_median"],
    keys: list[str],
):
    for key in keys:
        # Collect all values as
        values = np.ravel(distance_matrix[key])

        # Determine max values and median values
        max_values = np.max(values)
        median_values = np.median(values)
        tol = 1e-10 * np.max(values)
        non_trivial_median_values = np.median(values[values > tol])

        # Rescale data
        assert scaling_type in ["max", "median", "non_trivial_median"]
        if scaling_type == "max":
            scaling_value = 1.0 / max_values
        elif scaling_type == "median":
            scaling_value = 1.0 / median_values
        elif scaling_type == "non_trivial_median":
            scaling_value = 1.0 / non_trivial_median_values

        distance_matrix[key] = distance_matrix[key] * scaling_value


def nonlinear_transform(
    distance_matrix: np.ndarray,
    transform_type: Literal["linear", "log10", "symlog10"],
    inverse: bool = False,
):
    """Component-wise transform of the distance matrix."""

    if transform_type == "linear":
        return distance_matrix
    elif transform_type == "log10":
        if inverse:
            return 10**distance_matrix
        else:
            return np.log10(distance_matrix)
    elif transform_type == "log":
        if inverse:
            return np.exp(distance_matrix)
        else:
            return np.log(distance_matrix)
    elif transform_type in ["symlog10", "symlog"]:
        dm = distance_matrix.copy()
        # Check if dm is a scalar
        is_scalar = isinstance(dm, int) or isinstance(dm, float)
        is_list = isinstance(dm, list)
        if is_scalar:
            dm = np.array([dm])
        elif is_list:
            dm = np.array(dm)
        ind_gt_1 = distance_matrix > 1
        if inverse:
            dm[ind_gt_1] = (
                10 ** (dm[ind_gt_1] - 1)
                if transform_type == "symlog10"
                else np.exp(dm[ind_gt_1] - 1)
            )
        else:
            dm[ind_gt_1] = (
                1 + np.log10(dm[ind_gt_1])
                if transform_type == "symlog10"
                else 1 + np.log(dm[ind_gt_1])
            )
        if is_scalar:
            return dm[0]
        elif is_list:
            return dm.tolist()
        else:
            return dm
    else:
        raise ValueError("Transform type not supported.")


# ! ---- DATA ANALYSIS ----


def find_distance(distance_matrix, key, labels):
    """Find the distance between two participants."""
    i = labels.index(key[0])
    j = labels.index(key[1])
    return distance_matrix[i, j]


def reduce_distance_matrix_to_subset(
    distance_matrix: np.ndarray,
    subset_participating_groups: list,
    all_participating_groups: list,
):
    """Reduce square matrix to rows and columns/participants of interest."""
    group_indices = np.array(
        [all_participating_groups.index(g) for g in subset_participating_groups]
    ).astype(int)
    assert len(group_indices) == len(subset_participating_groups), (
        "expected length mismatch"
    )
    return distance_matrix[group_indices, :][:, group_indices]


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
        plt.savefig(path, dpi=300)
        print(f"Saved heatmap to {path}.")

    plt.show()
