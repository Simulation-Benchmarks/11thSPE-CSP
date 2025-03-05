"""Standardized workflow for data analysis - under development."""

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from icecream import ic
from scipy.spatial.distance import squareform
import pandas as pd

from cluster_analysis import (
    argmin_distance,
    centroid_analysis,
    determine_median_cluster,
    mean_distance_to_group,
    std_distance_to_group,
)
from datastructure import SPECase, convert_result_name, update_team_name
from spe11_io import (
    find_distance,
    identify_sparse_data,
    reduce_distance_matrix_to_subset,
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

    # Data management
    parser.add_argument(
        "-option",
        "--option",
        help="the type of the SPE11 analysis.",
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
        "-verbose",
        "--verbose",
        action="store_true",
        help="increase output verbosity",
    )

    args = parser.parse_args()

    # Define SPECase
    variant = args.variant
    spe_case = SPECase(variant)

    # Fetch type of analysis
    option = args.option

    # Data management
    data_folder = Path(args.folder)
    save_folder = Path(args.output)
    cache_folder = save_folder / "cache"

    # Group management
    groups = [g.lower() for g in args.groups]

    # Verbosity
    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    # ! ---- READ SPARSE DATA ----

    # Set up team structure - identify teams and submitted results based on the available sparse data
    participants = {}
    participants = identify_sparse_data(data_folder, spe_case, participants)

    # Read distance matrix using pandas - note it contains a header and index
    distance_matrix_with_labels = pd.read_csv(
        save_folder / f"{spe_case.variant}_distance_matrix.csv",
        delimiter=",",
        index_col=0,
    )
    distance_matrix = distance_matrix_with_labels.to_numpy()
    all_groups = distance_matrix_with_labels.columns.to_numpy().tolist()
    all_groups_lower = [g.lower() for g in all_groups]

    # ! ---- MINIMUM AG MEAN DISTANCE ON ALL SUBMISSIONS ----

    if "median-analysis" in args.option:
        print()
        print("-" * 80)
        print("MINIMUM AG MEAN DISTANCE ON ALL SUBMISSIONS:")
        print("-" * 80)
        print()

        mean_distance_all = {
            key: mean_distance_to_group(
                distance_matrix,
                key.lower(),
                all_groups_lower,
                mean_type="ag",
            )
            for key in all_groups
        }
        std_distance_all = {
            key: std_distance_to_group(
                distance_matrix,
                key.lower(),
                all_groups_lower,
                mean_type="ag",
            )
            for key in all_groups
        }
        print(
            f"""The group with the smallest mean distance to all other participants is """
            f"""{argmin_distance(mean_distance_all)}.\n"""
            f""" With distance {mean_distance_all[argmin_distance(mean_distance_all)]} and """
            f"""std {std_distance_all[argmin_distance(mean_distance_all)]}.\n"""
        )

    # ! ---- MEDIAN ANALYSIS BASED ON ALL SUBMISSIONS ----

    if args.option == "median-centroid-analysis":
        # ! ---- GROUP WITH MINIMUM MEAN DISTANCE TO ALL PARTICIPANTS ----

        # Quantitative inspection of the "median" data
        median_cluster_all = determine_median_cluster(
            distance_matrix,
            all_groups,
            mean_type="ag",
        )

        # Determine the mean distances of the two representatives of the median cluster to the other participants
        mean_distance_cluster_all = {
            key: mean_distance_to_group(
                distance_matrix,
                key,
                all_groups,
                mean_type="ag",
            )
            for key in median_cluster_all
        }
        print(
            f"The median cluster within all submissions consists of {median_cluster_all}, and {argmin_distance(mean_distance_cluster_all)} is the closest to all other participants."
        )

        centroid_analysis(
            distance_matrix,
            all_groups,
            mean_type="ag",
        )

    # ! ---- VARIABILITY ANALYSIS BASE CASE ----

    assert False

    if args.option == "variability-analysis-base-case":
        ...

    # Extract sub-matrix for the base case
    base_case_indices = [
        spe_case.subgroups["all"].index(key) for key in spe_case.subgroups["base case"]
    ]
    subgroups_global_distance_matrix = reduce_distance_matrix_to_subset(
        distance_matrix, spe_case.subgroups["base case"], all_groups
    )

    # ! ---- CREATE A DISCRETE PROBABILIY DISTRIBUTION FROM ALL SUBMISSIONS ----

    variability_ag_srg = variability_analysis(
        subgroups_global_distance_matrix["base case"], mean_type="ag"
    )

    # ! ---- VARIABILITY ANALYSIS ALL ----

    variability_ag_all = variability_analysis(distance_matrix, mean_type="ag")

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
        subgroups_global_distance_matrix["base case"], mean_type="ag"
    )

    variability_ag_srg_commercial = variability_analysis(
        subgroups_global_distance_matrix["base case commercial"],
        mean_type="ag",
    )

    variability_ag_srg_non_commercial = variability_analysis(
        subgroups_global_distance_matrix["base case non-commercial"],
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
            subgroups_global_distance_matrix["tetratech"]
        )
        print(
            f"The distance between the two tetratech submissions is {distance_values_tetratech}."
        )

        # ! ---- VARIABILITY ANALYSIS UT-CSEE (DIFFERENT SIMULATORS) ----

        distance_values_ut_csee = squareform(
            subgroups_global_distance_matrix["ut-csee"]
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
            distribution=np.sort(squareform(distance_matrix)),
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
                squareform(subgroups_global_distance_matrix["opm"]).tolist()
                + squareform(subgroups_global_distance_matrix["slb-IX"]).tolist()
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
                squareform(subgroups_global_distance_matrix["dumux"]).tolist()
                + squareform(subgroups_global_distance_matrix["slb-IX"]).tolist()
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
    #     subgroups_global_distance_matrix["refined-vs-base-case"],
    #     spe_case.subgroups["refined-vs-base-case"],
    #     linkage_type,
    #     spe_case,
    #     # f"{spe_case.variant}_dendrogram-refined-vs-base-case",
    # )

    if spe_case.variant in ["spe11a"]:
        print(
            "OPM1 vs OPM4:",
            find_distance(
                distance_matrix,
                ("opm1", "opm4"),
                spe_case.subgroups["all"],
            ),
        )
        print(
            "GEOS1 vs GEOS2:",
            find_distance(
                distance_matrix,
                ("geos1", "geos2"),
                spe_case.subgroups["all"],
            ),
        )

        print(
            "GEOS2 vs OPM4:",
            find_distance(
                distance_matrix,
                ("geos2", "opm4"),
                spe_case.subgroups["all"],
            ),
        )

    if spe_case.variant in ["spe11b"]:
        variability_mesh_opengosim = variability_analysis(
            subgroups_global_distance_matrix["mesh_refinement_study_opengosim"],
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
                    distance_matrix,
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
                    distance_matrix,
                    (key_i, key_j),
                    spe_case.subgroups["all"],
                )
                for key_i, key_j in [("sintef1", "sintef3"), ("sintef2", "sintef4")]
            ],
            "ifpen": find_distance(
                distance_matrix,
                ("ifpen1", "ifpen2"),
                spe_case.subgroups["all"],
            ),
            "geos": find_distance(
                distance_matrix,
                ("geos1", "geos2"),
                spe_case.subgroups["all"],
            ),
            "opm": find_distance(
                distance_matrix,
                ("opm1", "opm4"),
                spe_case.subgroups["all"],
            ),
            # "rice": find_distance(
            #    distance_matrix,
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
            subgroups_global_distance_matrix["100k-cartesian-mesh"],
            mean_type="ag",
        )
        variability_1M = variability_analysis(
            subgroups_global_distance_matrix["1.6m-cartesian-mesh"],
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
        #     subgroups_global_distance_matrix["facies-adapted-vs-base-case"],
        #     spe_case.subgroups["facies-adapted-vs-base-case"],
        #     linkage_type,
        #     spe_case,
        #     # f"{spe_case.variant}_dendrogram-refined-vs-base-case",
        # )
    if spe_case.variant in ["spe11c"]:
        print(
            "OPM1 vs CAU-KIEL:",
            find_distance(
                distance_matrix,
                ("opm1", "cau-kiel"),
                spe_case.subgroups["all"],
            ),
        )
        print(
            "SLB1 vs CAU-KIEL:",
            find_distance(
                distance_matrix,
                ("opm1", "cau-kiel"),
                spe_case.subgroups["all"],
            ),
        )
