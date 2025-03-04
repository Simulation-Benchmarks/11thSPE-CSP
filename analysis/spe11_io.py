"""Collected I/O functionality."""

import logging
from pathlib import Path
from typing import Literal
from warnings import warn

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform

from datastructure import SPECase

logger = logging.getLogger(__name__)


def get_result_name_from_multiple(folder):
    # folder has the structure .../team/speID/result<N> and is a Path object.
    # Reduce to team<N>.
    team = Path(folder).parts[-3]
    result = Path(folder).parts[-1]
    result = result.split("result")[-1]

    return f"{team}{result}"


def get_result_name_from_unique(folder):
    # folder has the structure .../team/speID/result<N> and is a Path object.
    # Reduce to team.
    team = Path(folder).parts[-2]

    return f"{team}"


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


def identify_sparse_data(
    sparse_data_folder: Path,
    spe_case: SPECase,
    participants,
):
    """Identify participants based on the folder structure and existence of files."""

    # Find paths to relevant sparse data
    for team in sparse_data_folder.iterdir():
        spe_input_folder = team / spe_case.variant

        # Check if the folder "result1" is contained in spe_input_folder.
        # Then multiple results are present. Otherwise, only one result is present,
        # which is directly in spe_input_folder.
        if (spe_input_folder / "result1").exists():
            # Multiple results are present.
            for result in spe_input_folder.iterdir():
                if "result" not in result.name:
                    continue
                if (result / f"{spe_case.variant}_time_series.csv").exists():
                    result_name = get_result_name_from_multiple(result)
                    participants[result_name] = {
                        "sparse": result / f"{spe_case.variant}_time_series.csv",
                    }
                else:
                    warn(f"Skipping {team} {result.name}")
        elif (spe_input_folder / f"{spe_case.variant}_time_series.csv").exists():
            # Only one result is present
            result_name = get_result_name_from_unique(spe_input_folder)
            participants[result_name] = {
                "sparse": spe_input_folder / f"{spe_case.variant}_time_series.csv",
            }
        else:
            warn(f"Results missing - skipping {team}.")

    return participants


def identify_dense_data(
    dense_data_folder: Path,
    spe_case: SPECase,
    participants,
):
    """Identify participants based on the folder structure and existence of files."""
    ...


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


# def rescale_data(data, spe_case: SPECase):
#     """Rescale the data to the correct units."""
#     # Start with unit scaling
#     scaling = spe_case.scaling
#     if False:
#         ...
#     elif True:
#         # Scale data based on the weighted time integral
#         scaling.update({})
#     scaling_matrix = np.diag([scaling[key] for key in scaling.keys()])
#     rescaled_data = data @ scaling_matrix
#     return rescaled_data


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


def find_distance(distance_matrix, key, labels):
    """Find the distance between two participants."""
    i = labels.index(key[0])
    j = labels.index(key[1])
    return distance_matrix[i, j]


def reduce_distance_matrix_to_participating_groups(
    distance_matrix: np.ndarray, participating_groups: list, participant_index: dict
):
    """Reduce square matrix to rows and columns/participants of interest."""
    group_indices = np.array(
        [
            participant_index[g]
            for g in sorted(participating_groups)
            if g in participant_index
        ]
    )
    return distance_matrix[group_indices, :][:, group_indices]


# def read_spe11b_sparse_data(participants: dict, spe_case: SPECase):
#     # timeunit = 31536000
#     # IS = 50 * timeunit
#     # report_times = np.arange(0, 1000 + 0.1, 0.1) * timeunit
#     # clustering_type = [
#     #     linear_analysis,
#     #     linear_analysis,
#     #     linear_analysis,
#     #     linear_analysis,
#     #     linear_analysis,
#     # ]
#     # rel_tols_cluster = [1.2, 1.2, 1.2, 1.5, 1.4]
#     # plot_type = ["semilogx", "semilogx", "semilogx", "loglog", "semilogx"]
#     # plot_min_exp = [5, 3, -2, 2, 0]
#     # plot_max_exp = [7.7, 7, 4.2, 7.9179, 5]
#
#     for participant in participants.keys():
#         data = participants[participant]["data"]
#         timestamps = data[:, spe_case.data_format["t"]]
#         mass_boxA = (
#             data[:, data_format["mobA"]]
#             + data[:, data_format["immA"]]
#             + data[:, data_format["dissA"]]
#         )
#         mass_boxB = (
#             data[:, data_format["mobB"]]
#             + data[:, data_format["immB"]]
#             + data[:, data_format["dissB"]]
#         )
#         M_C = data[:, data_format["M_C"]]
#         sealTot = data[:, data_format["sealTot"]]
#         boundaryCO2 = data[:, data_format["boundaryCO2"]]
#     for key, value in data_groups.items():
#         data = participants[key]["data"]
#         timestamps = data[iter1].iloc[:, data_groups["time"]].values
#         resultstemp = data[iter1].iloc[:, data_groups[key]].values
