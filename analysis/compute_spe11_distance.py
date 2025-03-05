"""Standardized workflow for data analysis - under development."""

import argparse
import logging
from pathlib import Path
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
from icecream import ic
from scipy.spatial.distance import squareform

from analysis import (
    field_data_distance_matrix,
    mean_matrix,
    plot_heatmap_distance_matrix,
    scalar_data_series_to_distance_matrix,
)
from cluster_analysis import (
    argmin_distance,
    centroid_analysis,
    determine_median_cluster,
    mean_distance_to_group,
    plot_linkage_clustering_with_colors,
    std_distance_to_group,
)
from datastructure import SPECase, convert_result_name, update_team_name
from spe11_io import (
    determine_reference_value_distance_matrix,
    find_distance,
    identify_sparse_data,
    interpolate_data_reporting_times,
    read_extra_data,
    read_field_data_distance_matrix_snapshots,
    read_time_series_data,
    reduce_distance_matrix_to_subset,
    reduce_to_increase_over_time,
    replace_mC_values,
    replace_sparse_values,
    rescale_distance_matrix,
    set_zero_boundaryCO2_values,
    split_result_name,
)
from safety import (
    check_boundaryCO2_values,
    check_mC_values,
    check_sanity,
    clean_data,
)
from sensitivity_analysis import (
    pearson_correlation_analysis,
    plot_correlation_test,
)
from variability_analysis import (
    cummulative_distribution_function,
    left_tailed_test_for_alternative_hypothesis,
    two_tailed_test_for_alternative_hypothesis,
    variability_analysis,
)

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    # ! ---- SETUP ----

    # Parse the variant and the use of petsc
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
    groups = [g.lower() for g in args.groups]
    # groupfolders = args.groupfolders
    calculatedAB = [g.lower() for g in args.calculatedAB]
    calculatedC = [g.lower() for g in args.calculatedC]
    zero_boundaryCO2 = [g.lower() for g in args.calculated_zero_boundaryCO2]

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

        # # Replace nan values with zeros
        # participant["scalar_data"] = clean_data(
        #     participant["scalar_data"], replace_nan=True
        # )

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
    criteria = [
        # Sparse data
        "mobA",
        "immA",
        "dissA",
        "mobB",
        "immB",
        "dissB",
        "mC",
        "sealTot",
        # Dense data
        "late_pressure_l2",
        "early_pressure_l2s",
        "early_mass_w1",
        "late_mass_w1",
    ]
    # Add SPE11B and SPE11C specific data
    if spe_case.variant in ["spe11b", "spe11c"]:
        criteria += [
            "boundaryCO2",
            "early_temperature_l2s",
            "late_temperature_l2s",
        ]

    # Compute global metrics for the selected criteria
    assert len(set([matrix.shape for matrix in distance_matrix.values()])) == 1, (
        "Matrices have different shapes. Need to reduce matrices to participants"
    )

    # ! ---- MEDIAN-BASED RESCALING ----

    # Need to group some of the data, to build meaningful medians?
    # By using non-trivial values for the median, this is not true anymore.

    median_scaled_distance_matrix = {
        key: rescale_distance_matrix(matrix, "non_trivial_median")
        for key, matrix in distance_matrix.items()
    }

    median_values_distance_matrix = {
        key: determine_reference_value_distance_matrix(matrix, "non_trivial_median")
        for key, matrix in distance_matrix.items()
    }

    # ! ---- NONLINEAR TRANSFORM ----

    global_distance_matrix = mean_matrix(
        median_scaled_distance_matrix,
        keys=criteria,
        mean_type="ag",
    )

    # Restrict data to subgroups
    subgroups_global_distance_matrix = {
        subgroup: reduce_distance_matrix_to_subset(
            global_distance_matrix,
            spe_case.subgroups[subgroup],
            spe_case.subgroups["all"],
            # participant_index,
        )
        for subgroup in spe_case.subgroups.keys()
    }

    # ! ---- DISTANCE MATRIX VALUES -----
    plot_heatmap_distance_matrix(
        subgroups_global_distance_matrix["all"],
        [convert_result_name(r) for r in spe_case.subgroups["all"]],
        spe_case,
        add_mean_to_diagonal=True,
        mean_type="ag",
        path=save_folder / f"{spe_case.variant}_heatmap_distance_matrix_all.png",
    )

    # Visual inspection in form of a dendrogram
    plot_linkage_clustering_with_colors(
        spe_case,
        subgroups_global_distance_matrix["all"],
        [convert_result_name(r) for r in spe_case.subgroups["all"]],
        path=save_folder / f"{spe_case.variant}_dendrogram",
    )

    # Store the main distance matrix as csv file with labels as row and column headers
    labels = [convert_result_name(r) for r in spe_case.subgroups["all"]]
    print(labels)
    np.savetxt(
        save_folder / f"{spe_case.variant}_distance_matrix.csv",
        subgroups_global_distance_matrix["all"],
        delimiter=",",
        header=",".join(labels),
    )
    # Store the distance matrix as csv file with indices as row and column headers using pandas
    pd.DataFrame(
        subgroups_global_distance_matrix["all"],
        index=labels,
        columns=labels,
    ).to_csv(
        save_folder / f"{spe_case.variant}_distance_matrix.csv",
        sep=",",
        header=True,
        index=True,
    )

    # TODO store median values as csv file etc

    # Cache results for further analysis, in form of csv file
    for key, matrix in subgroups_global_distance_matrix.items():
        np.savetxt(
            cache_folder / f"{spe_case.variant}_distance_matrix_{key}.csv",
            matrix,
            delimiter=",",
        )

    # ! ---- SENSIITIVITY ANALYSIS ----

    # Correlation analysis
    for key in criteria:
        plot_correlation_test(
            median_scaled_distance_matrix[key],
            subgroups_global_distance_matrix["all"],
            key,
            "SPE11 distance",
        )

    pearson_correlation_result = pearson_correlation_analysis(
        {key: median_scaled_distance_matrix[key] for key in criteria},
        subgroups_global_distance_matrix["all"],
    )

    # Combine median values and pearson correlation analysis in a table
    # TODO

    # Monitor median values
    ic(median_values_distance_matrix)

    print("Pearson correlation analysis:")
    ic(pearson_correlation_result)
