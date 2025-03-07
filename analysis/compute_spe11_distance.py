"""Standardized workflow for data analysis - under development."""

import argparse
import logging
from pathlib import Path
import pandas as pd

import numpy as np

from datastructure import SPECase
from utilities import (
    field_data_distance_matrix,
    mean_matrix,
    scalar_data_series_to_distance_matrix,
    determine_reference_value_distance_matrix,
    identify_sparse_data,
    interpolate_data_reporting_times,
    read_extra_data,
    read_field_data_distance_matrix_snapshots,
    read_time_series_data,
    reduce_to_increase_over_time,
    replace_sparse_values,
    rescale_distance_matrix,
    set_zero_boundaryCO2_values,
    split_result_name,
    check_sanity,
)
from sensitivity_analysis import (
    pearson_correlation_analysis,
)
from cluster_analysis import mean

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    # ! ---- SETUP ----

    parser = argparse.ArgumentParser(
        description="This script performs the SPE11 data analysis for the submitted data."
    )
    parser.add_argument(
        "-v",
        "--variant",
        help="The variant of the SPE11 testcase.",
    )
    parser.add_argument(
        "-f",
        "--folder",
        default="../shared_folder/data",
        help="path to folder containing group subfolders",
    )

    parser.add_argument(
        "-t",
        "--tablefolder",
        default="../shared_folder/evaluation",
        help="path to folder containing calculated tables",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="output",
        help="path to folder for saving output",
    )

    parser.add_argument(
        "-g",
        "--groups",
        nargs="+",
        help="names of groups, taking reported numbers",
        required=True,
    )

    parser.add_argument(
        "-gf",
        "--groupfolders",
        nargs="+",
        help="paths to group folders",
        required=False,
    )

    parser.add_argument(
        "-cAB",
        "--calculatedAB",
        nargs="+",
        help="names of groups, taking calculated numbers for Boxes A and B",
    )

    parser.add_argument(
        "-cC",
        "--calculatedC",
        nargs="+",
        help="names of groups, taking calculated numbers for Box C",
    )

    parser.add_argument(
        "-cBCO2",
        "--calculated_zero_boundaryCO2",
        nargs="+",
        help="names of groups, where boundary CO2 values are set to zero",
    )

    parser.add_argument(
        "-verbose",
        "--verbose",
        action="store_true",
        help="increase output verbosity",
    )

    args = parser.parse_args()

    # Define SPECase
    variant = args.variant
    spe_case = SPECase(variant)

    # Data management
    data_folder = Path(args.folder)
    tablefolder = Path(args.tablefolder)
    sparse_tablefolder = tablefolder / spe_case.variant / "sparse"
    dense_tablefolder = tablefolder / spe_case.variant / "dense"
    save_folder = Path(args.output)
    cache_folder = save_folder / "cache"
    save_folder.mkdir(parents=True, exist_ok=True)
    cache_folder.mkdir(parents=True, exist_ok=True)

    # Group management
    groups = args.groups
    groups_lower = [g.lower() for g in args.groups]
    # groupfolders = args.groupfolders
    calculatedAB = [g.lower() for g in args.calculatedAB]
    calculatedC = [g.lower() for g in args.calculatedC]
    if args.calculated_zero_boundaryCO2:
        zero_boundaryCO2 = [g.lower() for g in args.calculated_zero_boundaryCO2]
    else:
        zero_boundaryCO2 = []

    # Verbosity
    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    # ! ---- READ SPARSE DATA ----

    # Set up team structure - identify teams and submitted results based on the available sparse data
    participants = {}
    participants = identify_sparse_data(data_folder, spe_case, participants)

    # Read post-processed data - ImmA, ImmB, mC, sealA, sealB
    extra_values = {}
    for extra_key in [
        "immA",
        "immB",
        "mobA",
        "mobB",
        "dissA",
        "dissB",
        "mC",
        "sealA",
        "sealB",
    ]:
        extra_values_path = (
            sparse_tablefolder / f"{spe_case.variant}_{extra_key}_from_spatial_maps.csv"
        )
        extra_values[extra_key] = read_extra_data(extra_values_path, participants)

    # Post-process seal data - compute total seal data
    extra_values["sealTot"] = {}
    extra_values["sealTot"]["t"] = extra_values["sealA"]["t"].copy()
    extra_values["sealTot"].update(
        {
            key: extra_values["sealA"][key] + extra_values["sealB"][key]
            for key in participants.keys()
            if key != "t"
        }
    )

    # ! ---- READ DENSE DATA ----

    # Read data for each participant and clean the data
    errors = []
    for key, participant in participants.items():
        logging.info(f"Processing submission {key}")

        team, resultID = split_result_name(key)

        # Read raw data
        participant["scalar_data"] = read_time_series_data(
            participant["sparse"], spe_case
        )

        # Replace inconsistent data in boxes A and B
        if key in calculatedAB:
            for extra_key in [
                "immA",
                "immB",
                "dissA",
                "dissB",
                "mobA",
                "mobB",
                "sealA",
                "sealB",
                "sealTot",
            ]:
                participant["scalar_data"] = replace_sparse_values(
                    participant["scalar_data"],
                    extra_values[extra_key]["t"],
                    extra_values[extra_key][key],
                    extra_key,
                    spe_case,
                    key,
                )
        # Replace inconsistent data for the convection
        if key in calculatedC:
            for extra_key in ["mC"]:
                participant["scalar_data"] = replace_sparse_values(
                    participant["scalar_data"],
                    extra_values[extra_key]["t"],
                    extra_values[extra_key][key],
                    extra_key,
                    spe_case,
                    key,
                )

        # Set zero data for the boundary CO2.
        if key in zero_boundaryCO2:
            participant["scalar_data"] = set_zero_boundaryCO2_values(
                participant["scalar_data"], spe_case
            )

        # Check sanity of the data
        for data_type, index in spe_case.data_format.items():
            if data_type not in [
                "mobA",
                "mobB",
                "immA",
                "immB",
                "dissA",
                "dissB",
                "sealA",
                "sealB",
                "mC",
                "sealTot",
                "boundaryCO2",
            ]:
                continue
            if not check_sanity(participant["scalar_data"][:, index], key):
                raise ValueError(
                    f"Error in the data for {key} in {index}-th column - stop the analysis."
                )

    # Remove non-admissible participants
    for key in errors:
        participants.pop(key)

    # Sort the participants
    participants = dict(sorted(participants.items()))

    # Define participant -> index mapping
    participant_index = {key: i for i, key in enumerate(participants.keys())}

    # ! ---- SPARSE DATA DISTANCES ----

    # Conform data to reporting times and reduce to increase over time
    for key, participant in participants.items():
        # Interpolate data to reporting times
        participant["scalar_data"] = interpolate_data_reporting_times(
            participant["scalar_data"], spe_case
        )

        # Semi-norm view - restrict to increase over initial value
        participant["scalar_data"] = reduce_to_increase_over_time(
            participant["scalar_data"]
        )

    # Compute the sparse data distance tables for all cross combinations of participants, integrated over time.
    distance_matrix = {
        quantity: scalar_data_series_to_distance_matrix(
            participants,
            participant_index,
            quantity,
            spe_case,
        )
        for quantity in spe_case.data_format.keys()
    }

    # ! ---- DENSE DATA DISTANCES ----

    # Read evaluated field data distance matrices for selected snapshots
    distance_matrix_snapshots = read_field_data_distance_matrix_snapshots(
        dense_tablefolder, participants, spe_case
    )

    # Integrate in time
    distance_matrix.update(
        {
            timing + "_" + quantity: field_data_distance_matrix(
                participants,
                distance_matrix_snapshots,
                quantity,
                timing,
                spe_case,
            )
            for timing in ["early", "late", "full"]
            for quantity in ["pressure_l2", "pressure_l2s", "mass_w1"]
        }
    )
    if spe_case.non_isothermal:
        distance_matrix.update(
            {
                timing + "_" + quantity: field_data_distance_matrix(
                    participants,
                    distance_matrix_snapshots,
                    quantity,
                    timing,
                    spe_case,
                )
                for timing in ["early", "late", "full"]
                for quantity in [
                    "temperature_l2s",
                ]
            }
        )

    # ! ---- COMBINED DATA DISTANCES / GLOBAL METRIC ----

    # Define criteria for global metrics
    sparse_data_criteria = [
        # Sparse data
        "mobA",
        "immA",
        "dissA",
        "mobB",
        "immB",
        "dissB",
        "mC",
        "sealTot",
    ]
    # BoundaryCO2 data was not a reporting quantity for variant A
    if spe_case.variant in ["spe11b", "spe11c"]:
        sparse_data_criteria += ["boundaryCO2"]

    dense_data_criteria = [
        # Dense data
        "early_pressure_l2s",
        "late_pressure_l2",
        "early_mass_w1",
        "late_mass_w1",
    ]
    if spe_case.non_isothermal:
        dense_data_criteria += [
            "early_temperature_l2s",
            "late_temperature_l2s",
        ]
    criteria = sparse_data_criteria + dense_data_criteria

    # Compute global metrics for the selected criteria
    assert len(set([matrix.shape for matrix in distance_matrix.values()])) == 1, (
        "Matrices have different shapes. Need to reduce matrices to participants"
    )

    # ! ---- MEDIAN-BASED RESCALING ----

    median_values_distance_matrix = {
        key: determine_reference_value_distance_matrix(matrix, "non_trivial_median")
        for key, matrix in distance_matrix.items()
    }

    # Safety check for the median values - compare with spe_case.spe11_distance_median values
    for key, median_value in median_values_distance_matrix.items():
        assert np.isclose(median_value, spe_case.spe11_distance_median[key]), (
            f"Median value for {key} is {median_value}, expected value is {spe_case.spe11_distance_median[key]}"
        )

    median_scaled_distance_matrix = {
        key: rescale_distance_matrix(matrix, "non_trivial_median")
        for key, matrix in distance_matrix.items()
    }

    # ! ---- NONLINEAR TRANSFORM ----

    global_distance_matrix = mean_matrix(
        median_scaled_distance_matrix,
        keys=criteria,
        mean_type="ag",
    )

    # ! ---- ADD MEAN TO DIAGONAL ----

    for i in range(global_distance_matrix.shape[0]):
        distances = global_distance_matrix[i, :].tolist()
        distances.pop(i)
        distances = np.array(distances)
        global_distance_matrix[i, i] = mean(distances, mean_type="ag")

    # ! ---- STORE RESULTS -----

    # Store the distance matrix as csv file with indices as row and column headers using pandas
    pd.DataFrame(
        global_distance_matrix,
        index=groups,
        columns=groups,
    ).to_csv(
        save_folder / f"{spe_case.variant}_distance_matrix.csv",
        sep=",",
        header=True,
        index=True,
    )
    print()
    print(
        f"Distance matrix stored in {save_folder / f'{spe_case.variant}_distance_matrix.csv'}"
    )

    # ! ---- SENSIITIVITY ANALYSIS ----

    # Cross-Correlation analysis
    # for key in criteria:
    #    plot_correlation_test(
    #        median_scaled_distance_matrix[key],
    #        global_distance_matrix,
    #        key,
    #        "SPE11 distance",
    #    )

    pearson_correlation_result = pearson_correlation_analysis(
        {key: median_scaled_distance_matrix[key] for key in criteria},
        global_distance_matrix,
    )

    # Export a csv file with the criteria as labels in the rows, and two columns for the median values and the correlation values
    statistics_values = [
        [median_values_distance_matrix[key], pearson_correlation_result[key]]
        for key in criteria
    ]
    statistics_df = pd.DataFrame(
        statistics_values, index=criteria, columns=["Median", "PCC"]
    )
    statistics_df.to_csv(
        save_folder / f"{spe_case.variant}_statistics.csv",
        sep=",",
        header=True,
        index=True,
    )
    print(f"Statistics stored in {save_folder / f'{spe_case.variant}_statistics.csv'}")
