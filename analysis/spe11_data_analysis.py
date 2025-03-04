"""Standardized workflow for data analysis - under development."""

import argparse
import logging
from pathlib import Path

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
    reduce_distance_matrix_to_participating_groups,
    reduce_to_increase_over_time,
    replace_M_values,
    replace_sparse_values,
    rescale_distance_matrix,
    set_zero_boundaryCO2_values,
    split_result_name,
)
from safety import (
    check_boundaryCO2_values,
    check_M_values,
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

logging.basicConfig(level=logging.INFO)

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
        default="../shared_data/data",
        help="path to folder containing group subfolders",
    )

    parser.add_argument(
        "-t",
        "--tablefolder",
        default="../shared_data/evaluation",
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

    args = parser.parse_args()

    # Define SPECase
    variant = args.variant
    spe_case = SPECase(variant)

    # Data management
    data_folder = Path(args.folder)
    evaluation_folder = Path(args.tablefolder)
    sparse_evaluation_folder = evaluation_folder / spe_case.variant / "sparse"
    dense_evaluation_folder = evaluation_folder / spe_case.variant / "dense"
    save_folder = Path(args.output)
    cache_folder = save_folder / "cache"
    save_folder.mkdir(parents=True, exist_ok=True)
    cache_folder.mkdir(parents=True, exist_ok=True)

    # ! ---- READ SPARSE DATA ----

    # Set up team structure - identify teams and submitted results based on the available sparse data
    participants = {}
    participants = identify_sparse_data(data_folder, spe_case, participants)

    # Read extra data - compute convection from spatial maps using the official SPE11-CSP git repo
    immA_values_path = (
        sparse_evaluation_folder / f"{spe_case.variant}_immA_from_spatial_maps.csv"
    )
    immB_values_path = (
        sparse_evaluation_folder / f"{spe_case.variant}_immB_from_spatial_maps.csv"
    )
    m_values_path = (
        sparse_evaluation_folder / f"{spe_case.variant}_mC_from_spatial_maps.csv"
    )
    sealA_values_path = (
        sparse_evaluation_folder / f"{spe_case.variant}_sealA_from_spatial_maps.csv"
    )
    sealB_values_path = (
        sparse_evaluation_folder / f"{spe_case.variant}_sealB_from_spatial_maps.csv"
    )

    immA_values = read_extra_data(immA_values_path, participants)
    immB_values = read_extra_data(immB_values_path, participants)
    m_values = read_extra_data(m_values_path, participants)
    sealA_values = read_extra_data(sealA_values_path, participants)
    sealB_values = read_extra_data(sealB_values_path, participants)

    sealTot_values = {}
    sealTot_values["t"] = sealA_values["t"].copy()
    sealTot_values.update(
        {
            key: sealA_values[key] + sealB_values[key]
            for key in participants.keys()
            if key != "t"
        }
    )

    # ! ---- READ DENSE DATA ----

    # Read data for each participant and clean the data
    errors = []
    for key, participant in participants.items():
        print("Processing submission", key)

        team, resultID = split_result_name(key)

        # Read raw data
        participant["scalar_data"] = read_time_series_data(
            participant["sparse"], spe_case
        )

        # Replace inconsistent data
        if key in spe_case.inconsistent_immobile_saturations:
            participant["scalar_data"] = replace_sparse_values(
                participant["scalar_data"],
                immA_values["t"],
                immA_values[key],
                "immA",
                spe_case,
                key,
            )
            participant["scalar_data"] = replace_sparse_values(
                participant["scalar_data"],
                immB_values["t"],
                immB_values[key],
                "immB",
                spe_case,
                key,
            )
        if key in spe_case.inconsistent_convection:
            participant["scalar_data"] = replace_sparse_values(
                participant["scalar_data"],
                m_values["t"],
                m_values[key],
                "M_C",
                spe_case,
                key,
            )
        if key in spe_case.inconsistent_seal:
            participant["scalar_data"] = replace_sparse_values(
                participant["scalar_data"],
                sealTot_values["t"],
                sealTot_values[key],
                "sealTot",
                spe_case,
                key,
            )
        if key in spe_case.inconsistent_boundary_co2:
            participant["scalar_data"] = set_zero_boundaryCO2_values(
                participant["scalar_data"], spe_case
            )

        # Clean data
        participant["scalar_data"] = clean_data(
            participant["scalar_data"], replace_nan=False
        )

        # Check availbiity of M_values - replace with precomputed values if necessary
        status_m_values = check_M_values(participant["scalar_data"], spe_case, key)
        if not status_m_values:
            raise ValueError("M-values are missing - stop the analysis.")
            participant["scalar_data"] = replace_M_values(
                participant["scalar_data"],
                m_values["t"],
                m_values[key],
                spe_case=spe_case,
                key=key,
            )
            assert check_M_values(participant["scalar_data"], spe_case, key), (
                "Replacement of M-values failed - recomputation required."
            )

        # Check availability of boundary CO2 data - data not required for the analysis at the moment
        status_boundary_co2_values = check_boundaryCO2_values(
            participant["scalar_data"], spe_case, key
        )
        if not status_boundary_co2_values:
            raise ValueError("Boundary CO2 values are missing - stop the analysis.")
            logger.info(f"Set zero boundary CO2 values for submission {key}.")
            participant["scalar_data"] = set_zero_boundaryCO2_values(
                participant["scalar_data"], spe_case
            )

        # Replace nan values with zeros
        participant["scalar_data"] = clean_data(
            participant["scalar_data"], replace_nan=True
        )

        status_sanity = check_sanity(participant["scalar_data"], key)
        if not status_sanity:
            raise ValueError("Error in the data - stop the analysis.")
            errors.append(key)
            logging.info(f"Error in data - exclude submission {key}.")
            continue

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
        dense_evaluation_folder, participants, spe_case
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
    if spe_case.variant in ["spe11b", "spe11c"]:
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
    criteria = {
        "all": [
            # Sparse data
            "mobA",
            "immA",
            "dissA",
            "mobB",
            "immB",
            "dissB",
            "M_C",
            "sealTot",
            # Dense data
            "late_pressure_l2",
            "early_pressure_l2s",
            "early_mass_w1",
            "late_mass_w1",
        ]
        + (
            ["boundaryCO2", "early_temperature_l2s", "late_temperature_l2s"]
            if spe_case.variant in ["spe11b", "spe11c"]
            else []
        ),
        # Sparse vs. dense
        "sparse": [
            "mobA",
            "immA",
            "dissA",
            "mobB",
            "immB",
            "dissB",
            "M_C",
            "sealTot",
        ]
        + (["boundaryCO2"] if spe_case.variant in ["spe11b", "spe11c"] else []),
        "dense": [
            "early_pressure_l2",
            "late_pressure_l2",
            "early_pressure_l2s",
            "late_pressure_l2s",
            "early_mass_w1",
            "late_mass_w1",
        ]
        + (
            ["early_temperature_l2s", "late_temperature_l2s"]
            if spe_case.variant in ["spe11b", "spe11c"]
            else []
        ),
    }

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

    global_distance_matrix = {
        key: mean_matrix(
            median_scaled_distance_matrix,
            keys=criteria[key],
            mean_type="ag",
        )
        for key in criteria.keys()
    }

    # Restrict data to subgroups
    subgroups_global_distance_matrix = {
        (subgroup, key): reduce_distance_matrix_to_participating_groups(
            matrix,
            spe_case.subgroups[subgroup],
            participant_index,
        )
        for key, matrix in global_distance_matrix.items()
        for subgroup in spe_case.subgroups.keys()
    }

    # ! ---- DISTANCE MATRIX VALUES -----
    plot_heatmap_distance_matrix(
        subgroups_global_distance_matrix[("all", "all")],
        [convert_result_name(r) for r in spe_case.subgroups["all"]],
        spe_case,
        add_mean_to_diagonal=True,
        mean_type="ag",
        path=save_folder / f"{spe_case.variant}_heatmap_distance_matrix_all.png",
    )

    # ! ---- SENSIITIVITY ANALYSIS ----

    # Monitor median values
    ic(median_values_distance_matrix)

    # Correlation analysis
    for key in criteria["all"]:
        plot_correlation_test(
            median_scaled_distance_matrix[key],
            subgroups_global_distance_matrix[("all", "all")],
            key,
            "SPE11 distance",
        )

    pearson_correlation_result = pearson_correlation_analysis(
        {key: median_scaled_distance_matrix[key] for key in criteria["all"]},
        subgroups_global_distance_matrix[("all", "all")],
    )

    print("Pearson correlation analysis:")
    ic(pearson_correlation_result)

    # ! ---- MINIMUM AG MEAN DISTANCE ON ALL SUBMISSIONS ----

    print()
    print("-" * 80)
    print("MINIMUM AG MEAN DISTANCE ON ALL SUBMISSIONS:")
    print("-" * 80)
    print()

    mean_distance_all = {
        key: mean_distance_to_group(
            subgroups_global_distance_matrix[("all", "all")],
            key,
            spe_case.subgroups["all"],
            mean_type="ag",
        )
        for key in participants.keys()
    }
    std_distance_all = {
        key: std_distance_to_group(
            subgroups_global_distance_matrix[("all", "all")],
            key,
            spe_case.subgroups["all"],
            mean_type="ag",
        )
        for key in participants.keys()
    }
    print(
        f"""The group with the smallest mean distance to all other participants is """
        f"""{argmin_distance(mean_distance_all)} with distance """
        f"""{mean_distance_all[argmin_distance(mean_distance_all)]} and std """
        f"""{std_distance_all[argmin_distance(mean_distance_all)]}."""
    )

    # ! ---- MINIMUM MEAN DISTANCE ON SRG ----

    mean_distance_srg = {
        key: mean_distance_to_group(
            subgroups_global_distance_matrix[("all", "all")],
            key,
            spe_case.subgroups["all"],
            mean_type="ag",
        )
        for key in spe_case.subgroups["base case"]
    }
    print(
        f"The SRG group with the smallest mean distance to all (!) other participants is {argmin_distance(mean_distance_srg)}"
    )

    # ! ---- MEDIAN ANALYSIS BASED ON ALL SUBMISSIONS ----

    # Visual inspection of the "median" data
    plot_linkage_clustering_with_colors(
        spe_case,
        subgroups_global_distance_matrix[("all", "all")],
        [convert_result_name(r) for r in spe_case.subgroups["all"]],
        path=save_folder / f"{spe_case.variant}_dendrogram",
    )

    # Quantitative inspection of the "median" data
    median_cluster_all = determine_median_cluster(
        subgroups_global_distance_matrix[("all", "all")],
        spe_case.subgroups["all"],
        mean_type="ag",
    )

    # Determine the mean distances of the two representatives of the median cluster to the other participants
    mean_distance_cluster_all = {
        key: mean_distance_to_group(
            subgroups_global_distance_matrix[("all", "all")],
            key,
            spe_case.subgroups["all"],
            mean_type="ag",
        )
        for key in median_cluster_all
    }
    print(
        f"The median cluster within all submissions consists of {median_cluster_all}, and {argmin_distance(mean_distance_cluster_all)} is the closest to all other participants."
    )

    centroid_analysis(
        subgroups_global_distance_matrix[("all", "all")],
        spe_case.subgroups["all"],
        mean_type="ag",
    )

    # ! ---- CREATE A DISCRETE PROBABILIY DISTRIBUTION FROM ALL SUBMISSIONS ----

    probability_distribution_all = np.sort(
        squareform(subgroups_global_distance_matrix[("all", "all")])
    )
    variability_ag_all = variability_analysis(
        subgroups_global_distance_matrix[("all", "all")],
        mean_type="ag",
    )

    probability_distribution_srg = np.sort(
        squareform(subgroups_global_distance_matrix[("base case", "all")])
    )
    variability_ag_srg = variability_analysis(
        subgroups_global_distance_matrix[("base case", "all")], mean_type="ag"
    )

    # ! ---- GROUP WITH MINIMUM MEAN DISTANCE TO ALL PARTICIPANTS ----

    mean_distance_all = {
        key: mean_distance_to_group(
            subgroups_global_distance_matrix[("all", "all")],
            key,
            spe_case.subgroups["all"],
            mean_type="ag",
        )
        for key in participants.keys()
    }
    print(
        f"The group with the smallest mean distance to all other participants is {argmin_distance(mean_distance_all)}"
    )

    # ! ---- VARIABILITY ANALYSIS ALL ----

    variability_ag_all = variability_analysis(
        subgroups_global_distance_matrix[("all", "all")], mean_type="ag"
    )

    print()
    print("-" * 80)
    print("VARIABILITY ANALYSIS ALL SUBMISSIONS:")
    print("-" * 80)
    print()
    print(
        f"Mean (ag) distance within all submissions is {variability_ag_all['mean']} +- {variability_ag_all['margin_of_error']} with a standard deviation of {variability_ag_all['std']}."
    )

    # ! ---- VARIABILITY ANALYSIS SRG ----

    variability_ag_srg = variability_analysis(
        subgroups_global_distance_matrix[("base case", "all")], mean_type="ag"
    )

    variability_ag_srg_commercial = variability_analysis(
        subgroups_global_distance_matrix[("base case commercial", "all")],
        mean_type="ag",
    )

    variability_ag_srg_non_commercial = variability_analysis(
        subgroups_global_distance_matrix[("base case non-commercial", "all")],
        mean_type="ag",
    )

    print()
    print("-" * 80)
    print("VARIABILITY ANALYSIS SRG:")
    print("-" * 80)
    print()
    print(
        f"Mean (ag) distance within SRG is {variability_ag_srg['mean']} +- {variability_ag_srg['margin_of_error']} with a standard deviation of {variability_ag_srg['std']}."
    )
    print()
    print()
    print(
        f"Mean (ag) distance within SRG (commercial) is {variability_ag_srg_commercial['mean']} +- {variability_ag_srg_commercial['margin_of_error']} with a standard deviation of {variability_ag_srg['std']}."
    )
    print()
    print()
    print(
        f"Mean (ag) distance within SRG (non-commercial) is {variability_ag_srg_non_commercial['mean']} +- {variability_ag_srg_non_commercial['margin_of_error']} with a standard deviation of {variability_ag_srg['std']}."
    )
    print()

    # ! ---- TEST NULL HYPOTHESIS THAT BASELINE VARIABILITY < MEAN ----
    p_srg_smaller_all = left_tailed_test_for_alternative_hypothesis(
        statistics_smaller=variability_ag_srg,
        statistics_greater=variability_ag_all,
        transformed=True,
    )
    p_srg_commericial_smaller_srg = left_tailed_test_for_alternative_hypothesis(
        statistics_smaller=variability_ag_srg_commercial,
        statistics_greater=variability_ag_srg,
        transformed=True,
    )
    p_srg_non_commercial_smaller_srg = left_tailed_test_for_alternative_hypothesis(
        statistics_smaller=variability_ag_srg_non_commercial,
        statistics_greater=variability_ag_srg,
        transformed=True,
    )
    print()
    print("-" * 80)
    print("HYPOTHESIS TESTING - SRG < ALL:")
    print("-" * 80)
    print()
    print(
        f"Hypothesis that the variability within the SRG is smaller than the variability within all submissions has a p-value of {p_srg_smaller_all}.\n"
    )
    print()
    print(
        f"Hypothesis that the variability within the SRG (commercial) is smaller than the SRG variability has a p-value of {p_srg_commericial_smaller_srg}.\n"
    )
    print()
    print(
        f"Hypothesis that the variability within the SRG (non-commercial) is smaller than the SRG variability has a p-value of {p_srg_non_commercial_smaller_srg}.\n"
    )

    # ! ---- DIFFERENT SIMULATORS ----

    if spe_case.variant in ["spe11a"]:
        print("No different simulators available for variant A.")
    elif spe_case.variant in ["spe11b"]:
        # ! ---- VARIABILITY ANALYSIS TETRATECH (DIFFERENT SIMULATORS) ----

        print()
        print("-" * 80)
        print("ANALYSIS SAME GROUP - SAME GROUP / DIFFERENT SIMULATOR:")
        print("-" * 80)
        print()
        distance_values_tetratech = squareform(
            subgroups_global_distance_matrix[("tetratech", "all")]
        )
        print(
            f"The distance between the two tetratech submissions is {distance_values_tetratech}."
        )

        # ! ---- VARIABILITY ANALYSIS UT-CSEE (DIFFERENT SIMULATORS) ----

        distance_values_ut_csee = squareform(
            subgroups_global_distance_matrix[("ut-csee", "all")]
        )
        print(
            f"The distance between the two ut-csee submissions is {distance_values_ut_csee}."
        )

        variability_analysis_different_simulators = variability_analysis(
            np.array([distance_values_tetratech, distance_values_ut_csee]).flatten(),
            mean_type="ag",
        )

        p_value_different_simulator_smaller_srg = (
            left_tailed_test_for_alternative_hypothesis(
                statistics_smaller=variability_analysis_different_simulators,
                statistics_greater=variability_ag_srg,
                transformed=True,
            )
        )
        probability_different_simulator_smaller_srg = cummulative_distribution_function(
            value=max(
                np.max(distance_values_tetratech), np.max(distance_values_ut_csee)
            ),
            distribution=np.sort(
                squareform(subgroups_global_distance_matrix[("all", "all")])
            ),
        )

        print()
        print("-" * 80)
        print("HYPOHESIS TESTING - DIFFERENT SIMULATORS < SRG:")
        print("-" * 80)
        print()
        print(
            f"Hypothesis that the variability within the tetratech and ut-csee submissions is smaller than the one within the SRG has a p-value of {p_value_different_simulator_smaller_srg}.\n"
        )
        print(
            f"The probability that a variability is smaller then the variability in the tetratech and ut-csee submissions is {probability_different_simulator_smaller_srg} after discretization.\n"
        )

    # ! ---- DIFFERENT GROUPS ----

    if spe_case.variant in ["spe11a"]:
        # ! ---- VARIABILITY ANALYSIS OPM (DIFFERENT GROUPS) ----

        print()
        print("-" * 80)
        print("ANALYSIS SAME GROUP - DIFFERENT GROUPS / SAME SIMULATOR:")
        print("-" * 80)
        print()

        print(f"Distances for all participants using Dumux:")
        for key1, key2 in [
            (k1, k2)
            for k1 in spe_case.subgroups["opm"]
            for k2 in spe_case.subgroups["opm"]
            if k1 != k2
        ]:
            print(
                f"Distance between {key1} and {key2}: {find_distance(subgroups_global_distance_matrix[('all', 'all')], (key1, key2), spe_case.subgroups['all'])}"
            )

        # ! ---- VARIABILITY ANALYSIS SLB IX (DIFFERENT GROUPS) ----

        print(f"Distance matrix for all participants using SLB IX:")
        for key1, key2 in [
            (k1, k2)
            for k1 in spe_case.subgroups["slb-IX"]
            for k2 in spe_case.subgroups["slb-IX"]
            if k1 != k2
        ]:
            print(
                f"Distance between {key1} and {key2}: {find_distance(subgroups_global_distance_matrix[('all', 'all')], (key1, key2), spe_case.subgroups['all'])}"
            )

        variability_analysis_different_groups = variability_analysis(
            np.array(
                squareform(subgroups_global_distance_matrix[("opm", "all")]).tolist()
                + squareform(
                    subgroups_global_distance_matrix[("slb-IX", "all")]
                ).tolist()
            ).flatten(),
            mean_type="ag",
        )

        # ! ---- HYPOTHESIS TESTING - DIFFERENT GROUPS < SRG ----

        print()
        print("-" * 80)
        print("HYPOTHESIS TESTING - DIFFERENT GROUPS < SRG:")
        print("-" * 80)
        print()
        p_value_different_group_smaller_srg = (
            left_tailed_test_for_alternative_hypothesis(
                statistics_smaller=variability_analysis_different_groups,
                statistics_greater=variability_ag_srg,
                transformed=True,
            )
        )
        p_value_different_group_different_from_srg = (
            two_tailed_test_for_alternative_hypothesis(
                statistics_smaller=variability_analysis_different_groups,
                statistics_greater=variability_ag_srg,
                transformed=True,
            )
        )
        print(
            f"Hypothesis that the variability within the opm and slb-IX submissions is smaller than the one within the SRG has a p-value of {p_value_different_group_smaller_srg}.\n"
        )
        print(
            f"Hypothesis that the variability within the opm and slb-IX submissions is different from the one within the SRG has a p-value of {p_value_different_group_different_from_srg}.\n"
        )

    elif spe_case.variant in ["spe11b"]:
        # ! ---- VARIABILITY ANALYSIS DUMUX (DIFFERENT GROUPS) ----

        print()
        print("-" * 80)
        print("ANALYSIS SAME GROUP - DIFFERENT GROUPS / SAME SIMULATOR:")
        print("-" * 80)
        print()

        print(f"Distances for all participants using Dumux:")
        for key1, key2 in [
            (k1, k2)
            for k1 in spe_case.subgroups["dumux"]
            for k2 in spe_case.subgroups["dumux"]
            if k1 != k2
        ]:
            print(
                f"Distance between {key1} and {key2}: {find_distance(subgroups_global_distance_matrix[('all', 'all')], (key1, key2), spe_case.subgroups['all'])}"
            )

        # ! ---- VARIABILITY ANALYSIS SLB IX (DIFFERENT GROUPS) ----

        print(f"Distance matrix for all participants using SLB IX:")
        for key1, key2 in [
            (k1, k2)
            for k1 in spe_case.subgroups["slb-IX"]
            for k2 in spe_case.subgroups["slb-IX"]
            if k1 != k2
        ]:
            print(
                f"Distance between {key1} and {key2}: {find_distance(subgroups_global_distance_matrix[('all', 'all')], (key1, key2), spe_case.subgroups['all'])}"
            )

        variability_analysis_different_groups = variability_analysis(
            np.array(
                squareform(subgroups_global_distance_matrix[("dumux", "all")]).tolist()
                + squareform(
                    subgroups_global_distance_matrix[("slb-IX", "all")]
                ).tolist()
            ).flatten(),
            mean_type="ag",
        )

        # ! ---- HYPOTHESIS TESTING - DIFFERENT GROUPS < SRG ----

        print()
        print("-" * 80)
        print("HYPOTHESIS TESTING - DIFFERENT GROUPS < SRG:")
        print("-" * 80)
        print()
        p_value_different_group_smaller_srg = (
            left_tailed_test_for_alternative_hypothesis(
                statistics_smaller=variability_analysis_different_groups,
                statistics_greater=variability_ag_srg,
                transformed=True,
            )
        )
        p_value_different_group_different_from_srg = (
            two_tailed_test_for_alternative_hypothesis(
                statistics_smaller=variability_analysis_different_groups,
                statistics_greater=variability_ag_srg,
                transformed=True,
            )
        )
        print(
            f"Hypothesis that the variability within the dumux and slb-IX submissions is smaller than the one within the SRG has a p-value of {p_value_different_group_smaller_srg}.\n"
        )
        print(
            f"Hypothesis that the variability within the dumux and slb-IX submissions is different from the one within the SRG has a p-value of {p_value_different_group_different_from_srg}.\n"
        )

    else:
        raise NotImplementedError(f"Variant {spe_case.variant} not implemented.")

    # ! ---- MESH REFINEMENT ANALYSIS ----

    # plot_linkage_clustering_with_colors(
    #     subgroups_global_distance_matrix[("refined-vs-base-case", "all")],
    #     spe_case.subgroups["refined-vs-base-case"],
    #     linkage_type,
    #     spe_case,
    #     # f"{spe_case.variant}_dendrogram-refined-vs-base-case",
    # )

    if spe_case.variant in ["spe11a"]:
        print(
            "OPM1 vs OPM4:",
            find_distance(
                subgroups_global_distance_matrix[("all", "all")],
                ("opm1", "opm4"),
                spe_case.subgroups["all"],
            ),
        )
        print(
            "GEOS1 vs GEOS2:",
            find_distance(
                subgroups_global_distance_matrix[("all", "all")],
                ("geos1", "geos2"),
                spe_case.subgroups["all"],
            ),
        )

        print(
            "GEOS2 vs OPM4:",
            find_distance(
                subgroups_global_distance_matrix[("all", "all")],
                ("geos2", "opm4"),
                spe_case.subgroups["all"],
            ),
        )

    if spe_case.variant in ["spe11b"]:
        variability_mesh_opengosim = variability_analysis(
            subgroups_global_distance_matrix[
                ("mesh_refinement_study_opengosim", "all")
            ],
            mean_type="ag",
        )

        print()
        print("-" * 80)
        print("VARIABILITY ANALYSIS OPENGEOSIM SUBMISSIONS:")
        print("-" * 80)
        print()
        print(
            f"Mean (ag) distance within opengeosim submissions is {variability_mesh_opengosim['mean']} +- {variability_mesh_opengosim['margin_of_error']}."
        )

        # def mesh_refinement_study...
        mesh_refinement = {
            "opengosim": [
                (
                    spe_case.groups_and_mesh_size[key_j]
                    / spe_case.groups_and_mesh_size[key_i]
                )
                ** 0.5
                for key_i, key_j in [
                    ("opengosim1", "opengosim2"),
                    ("opengosim2", "opengosim3"),
                    ("opengosim1", "opengosim3"),
                ]
            ],
            "sintef": [
                (
                    spe_case.groups_and_mesh_size[key_j]
                    / spe_case.groups_and_mesh_size[key_i]
                )
                ** 0.5
                for key_i, key_j in [
                    ("sintef1", "sintef3"),
                    ("sintef2", "sintef4"),
                ]
            ],
            "ifpen": [
                (
                    spe_case.groups_and_mesh_size["ifpen2"]
                    / spe_case.groups_and_mesh_size["ifpen1"]
                )
                ** 0.5,
            ],
            "geos": [
                (
                    spe_case.groups_and_mesh_size["geos2"]
                    / spe_case.groups_and_mesh_size["geos1"]
                )
                ** 0.5,
            ],
            "opm": [
                (
                    spe_case.groups_and_mesh_size["opm4"]
                    / spe_case.groups_and_mesh_size["opm1"]
                )
                ** 0.5,
            ],
            # "rice": [
            #    (
            #        spe_case.groups_and_mesh_size["rice2"]
            #        / spe_case.groups_and_mesh_size["rice1"]
            #    )
            #    ** 0.5,
            # ],
        }
        mesh_refinement_info = {
            "opengosim": [
                "100k vs. 200k",  # "OpenGoSim1->2",
                "200k vs. 1.6M",  # "OpenGoSim2->3",
                "100k vs. 1.6M",  # "OpenGoSim1->3",
            ],
            "sintef": ["25k vs. 100k", "25k vs. 100k"],  # "SINTEF1->3", "SINTEF2->4"],
            "ifpen": ["100k vs. 1.6M"],  # ["IFPEN1->2"],
            "geos": ["100k vs. 1.6M"],
            "opm": ["100k vs. 10M"],
        }
        distance_values_refinement = {
            "opengosim": [
                find_distance(
                    subgroups_global_distance_matrix[("all", "all")],
                    (key_i, key_j),
                    spe_case.subgroups["all"],
                )
                for key_i, key_j in [
                    ("opengosim1", "opengosim2"),
                    ("opengosim2", "opengosim3"),
                    ("opengosim1", "opengosim3"),
                ]
            ],
            "sintef": [
                find_distance(
                    subgroups_global_distance_matrix[("all", "all")],
                    (key_i, key_j),
                    spe_case.subgroups["all"],
                )
                for key_i, key_j in [("sintef1", "sintef3"), ("sintef2", "sintef4")]
            ],
            "ifpen": find_distance(
                subgroups_global_distance_matrix[("all", "all")],
                ("ifpen1", "ifpen2"),
                spe_case.subgroups["all"],
            ),
            "geos": find_distance(
                subgroups_global_distance_matrix[("all", "all")],
                ("geos1", "geos2"),
                spe_case.subgroups["all"],
            ),
            "opm": find_distance(
                subgroups_global_distance_matrix[("all", "all")],
                ("opm1", "opm4"),
                spe_case.subgroups["all"],
            ),
            # "rice": find_distance(
            #    subgroups_global_distance_matrix[("all", "all")],
            #    ("rice1", "rice2"),
            #    spe_case.subgroups["all"],
            # ),
        }
        distance_values_moderate_static_refinement = {
            "opengosim": distance_values_refinement["opengosim"],
            "sintef": distance_values_refinement["sintef"],
            "ifpen": distance_values_refinement["ifpen"],
            "geos": distance_values_refinement["geos"],
        }
        variability_mesh_refinement_team_inter = variability_analysis(
            np.hstack(
                [distance_values_refinement[key] for key in distance_values_refinement]
            ).flatten(),
            mean_type="ag",
        )
        p_value_variability_mesh_refinement_team_inter_smaller_than_srg = (
            left_tailed_test_for_alternative_hypothesis(
                statistics_smaller=variability_mesh_refinement_team_inter,
                statistics_greater=variability_ag_srg,
                transformed=True,
            )
        )

        variability_moderate_static_refinement = variability_analysis(
            np.hstack(
                [
                    distance_values_moderate_static_refinement[key]
                    for key in distance_values_moderate_static_refinement
                ]
            ).flatten(),
            mean_type="ag",
        )

        plt.figure("Mesh refinement study")
        for key in sorted(mesh_refinement.keys()):
            sct = plt.scatter(
                mesh_refinement[key],
                np.array(distance_values_refinement[key]),
                marker="o",
                s=100,
                label=update_team_name(key),
                color=spe_case.groups_and_colors[key],
                edgecolors="gray",
                zorder=20,
            )
            for i, txt in enumerate(mesh_refinement_info[key]):
                try:
                    plt.text(
                        mesh_refinement[key][i] + 0.15,
                        distance_values_refinement[key][i] - 0.025,
                        txt,
                        fontsize=8,
                        color=[0.3] * 3,
                        zorder=30,
                    )
                except IndexError:
                    plt.text(
                        mesh_refinement[key][i] + 0.15,
                        distance_values_refinement[key] - 0.025,
                        txt,
                        fontsize=8,
                        color=[0.3] * 3,
                        zorder=30,
                    )

        # Add mean value line
        plt.axhline(
            y=variability_moderate_static_refinement["mean"],
            color="black",
            label="Mean variability|moderate refinement",
        )
        # Add margin of error
        plt.axhline(
            y=variability_moderate_static_refinement["lower_limit_margin_of_error"],
            color="gray",
            # linestyle="--",
        )
        plt.axhline(
            y=variability_moderate_static_refinement["upper_limit_margin_of_error"],
            color="gray",
            # linestyle="--",
        )
        plt.fill_between(
            [0, 10.5],
            variability_moderate_static_refinement["lower_limit_margin_of_error"],
            variability_moderate_static_refinement["upper_limit_margin_of_error"],
            color=[0.9] * 3,
        )
        # Add the same for the SRG
        plt.axhline(
            y=variability_ag_srg["mean"],
            color="black",
            linestyle="--",
            label="Baseline variability",
        )
        plt.axhline(
            y=variability_ag_srg["lower_limit_margin_of_error"],
            color="gray",
            linestyle="--",
        )
        plt.axhline(
            y=variability_ag_srg["upper_limit_margin_of_error"],
            color="gray",
            linestyle="--",
        )
        plt.fill_between(
            [0, 10.5],
            variability_ag_srg["lower_limit_margin_of_error"],
            variability_ag_srg["upper_limit_margin_of_error"],
            color=[0.9] * 3,
        )

        plt.legend(loc="upper left")
        plt.xlabel("Mesh refinement factor")
        plt.xticks([1, 2, 4, 8, 10])
        plt.yticks([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1])
        # Add grid behind the scatter plot
        plt.grid(zorder=10, color="lightgray")
        plt.xlim(0.5, 10.5)
        plt.ylim(0.25, 1.15)
        plt.ylabel("SPE11 distance")
        plt.tight_layout()
        plt.savefig(
            save_folder / f"{spe_case.variant}_mesh_refinement_study.png", dpi=1000
        )
        plt.show()

        p_value_variability_moderate_static_refinement_smaller_than_srg = (
            left_tailed_test_for_alternative_hypothesis(
                statistics_smaller=variability_moderate_static_refinement,
                statistics_greater=variability_ag_srg,
                transformed=True,
            )
        )

        print()
        print("-" * 80)
        print("VARIABILITY ANALYSIS MESH REFINEMENT STUDY:")
        print("-" * 80)
        print()
        print(
            f"Mean (ag) distance within mesh refinement study submissions is {variability_mesh_refinement_team_inter['mean']} +- {variability_mesh_refinement_team_inter['margin_of_error']} with a standard deviation of {variability_mesh_refinement_team_inter['std']}.\n"
        )
        print()
        print(
            f" Mean (ag) distance within moderate static refinement submissions is {variability_moderate_static_refinement['mean']} +- {variability_moderate_static_refinement['margin_of_error']} with a standard deviation of {variability_moderate_static_refinement['std']}.\n"
        )

        print()
        print("-" * 80)
        print("HYPOTHESIS TESTING - MESH REFINEMENT STUDY < SRG:")
        print("-" * 80)
        print()
        print(
            f"Hypothesis that the variability within the mesh refinement study submissions is smaller than the one within the SRG has a p-value of {p_value_variability_mesh_refinement_team_inter_smaller_than_srg}.\n"
        )
        print()

        print(
            f"Hypothesis that the variability within the moderate static refinement submissions is smaller than the one within the SRG has a p-value of {p_value_variability_moderate_static_refinement_smaller_than_srg}.\n"
        )

        # Variability within the two groups

        variability_100k = variability_analysis(
            subgroups_global_distance_matrix[("100k-cartesian-mesh", "all")],
            mean_type="ag",
        )
        variability_1M = variability_analysis(
            subgroups_global_distance_matrix[("1.6m-cartesian-mesh", "all")],
            mean_type="ag",
        )
        p_value_variability_100k_1M_different = (
            two_tailed_test_for_alternative_hypothesis(
                statistics_smaller=variability_100k,
                statistics_greater=variability_1M,
                transformed=True,
            )
        )

        print()
        print("-" * 80)
        print("VARIABILITY ANALYSIS 100k SUBMISSIONS:")
        print("-" * 80)
        print()
        print(
            f"Mean (ag) distance within 100k submissions is {variability_100k['mean']} +- {variability_100k['margin_of_error']}."
        )

        print()
        print("-" * 80)
        print("VARIABILITY ANALYSIS 1.6M SUBMISSIONS:")
        print("-" * 80)
        print()
        print(
            f"Mean (ag) distance within 1.6M submissions is {variability_1M['mean']} +- {variability_1M['margin_of_error']}."
        )

        print()
        print("-" * 80)
        print("VARIABILITY ANALYSIS 100k vs. 1.6M:")
        print("-" * 80)
        print()
        print(
            f"Hypothesis that the variability within the 100k submissions is different from the one within the 1.6M submissions has a p-value of {p_value_variability_100k_1M_different}.\n"
        )

        # ! ---- FACIES-ADAPTED MESHES ----

        # plot_linkage_clustering_with_colors(
        #     subgroups_global_distance_matrix[("facies-adapted-vs-base-case", "all")],
        #     spe_case.subgroups["facies-adapted-vs-base-case"],
        #     linkage_type,
        #     spe_case,
        #     # f"{spe_case.variant}_dendrogram-refined-vs-base-case",
        # )
    if spe_case.variant in ["spe11c"]:
        print(
            "OPM1 vs CAU-KIEL:",
            find_distance(
                subgroups_global_distance_matrix[("all", "all")],
                ("opm1", "cau-kiel"),
                spe_case.subgroups["all"],
            ),
        )
        print(
            "SLB1 vs CAU-KIEL:",
            find_distance(
                subgroups_global_distance_matrix[("all", "all")],
                ("opm1", "cau-kiel"),
                spe_case.subgroups["all"],
            ),
        )
