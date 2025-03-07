"""Standardized workflow for data analysis - under development."""

import argparse
import logging
from pathlib import Path

import numpy as np
from scipy.spatial.distance import squareform
import pandas as pd

from datastructure import SPECase
from utilities import (
    reduce_distance_matrix_to_subset,
)
from variability_analysis import (
    left_tailed_test_for_alternative_hypothesis,
    variability_analysis,
)

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
        help="names of groups to be analyzed",
        action="append",
        required=False,
    )

    parser.add_argument(
        "-g-smaller",
        "--smaller-groups",
        nargs="+",
        help="names of groups with smaller distance, taking reported numbers",
        action="append",
        required=False,
    )

    parser.add_argument(
        "-g-greater",
        "--greater-groups",
        nargs="+",
        help="names of groups with greater distance, taking reported numbers",
        action="append",
        required=False,
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
    save_folder = Path(args.output)
    cache_folder = save_folder / "cache"

    # Group management
    if args.groups:
        groups = args.groups
        groups_lower = [[g.lower() for g in group] for group in groups]
    else:
        groups = []
        groups_lower = []

    if args.smaller_groups:
        groups_smaller = args.smaller_groups
        groups_smaller_lower = [[g.lower() for g in group] for group in groups_smaller]
    else:
        groups_smaller = []
        groups_smaller_lower = []

    if args.greater_groups:
        groups_greater = args.greater_groups
        groups_greater_lower = [[g.lower() for g in group] for group in groups_greater]
    else:
        groups_greater = []
        groups_greater_lower = []

    # Safety check
    assert len(groups) > 0 or (len(groups_smaller) > 0 and len(groups_greater) > 0), (
        "No groups specified."
    )
    assert len(groups) == 0 or (
        len(groups_smaller) == 0 and len(groups_greater) == 0
    ), "Groups and smaller/greater groups cannot be specified at the same time."
    single_group_analysis = len(groups) > 0
    single_group_single_category_analysis = len(groups) == 1

    # Verbosity
    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    # ! ---- READ SPARSE DATA ----

    # Read distance matrix using pandas - note it contains a header and index
    distance_matrix_with_labels = pd.read_csv(
        save_folder / f"{spe_case.variant}_distance_matrix.csv",
        delimiter=",",
        index_col=0,
    )
    distance_matrix = distance_matrix_with_labels.to_numpy()
    all_groups = distance_matrix_with_labels.columns.to_numpy().tolist()
    all_groups_lower = [g.lower() for g in all_groups]

    if single_group_analysis:
        # ! ---- VARIABILITY ANALYSIS FOR SINGLE GROUP ----
        distances = []
        num_groups = 0
        for group in groups:
            reduced_distance_matrix = reduce_distance_matrix_to_subset(
                distance_matrix, group, all_groups
            )
            np.fill_diagonal(reduced_distance_matrix, 0)
            distances.append(squareform(reduced_distance_matrix))
            num_groups += len(group)
        distances = np.array(distances).flatten()

        variability_ag = variability_analysis(distances, mean_type="ag")

        print()
        print("-" * 80)
        print("VARIABILITY ANALYSIS:")
        print("-" * 80)
        print()
        print(
            f"The mean distance among the {num_groups} groups\n{groups}\n"
            ""
            f"""is {variability_ag["mean"]} +- {variability_ag["margin_of_error"]} (margin of error) based on \n"""
            f"""with a standard deviation of {variability_ag["std"]}."""
        )

    else:
        # ! ---- VARIABILITY FOR SMALLER AND GREATER GROUPS ----

        smaller_distances = []
        num_smaller_groups = 0
        for group in groups_smaller:
            reduced_distance_matrix = reduce_distance_matrix_to_subset(
                distance_matrix, group, all_groups
            )
            np.fill_diagonal(reduced_distance_matrix, 0)
            smaller_distances.append(squareform(reduced_distance_matrix))
            num_smaller_groups += len(group)
        smaller_distances = np.hstack(smaller_distances)

        greater_distances = []
        num_greater_groups = 0
        for group in groups_greater:
            reduced_distance_matrix = reduce_distance_matrix_to_subset(
                distance_matrix, group, all_groups
            )
            np.fill_diagonal(reduced_distance_matrix, 0)
            greater_distances.append(squareform(reduced_distance_matrix))
            num_greater_groups += len(group)
        greater_distances = np.hstack(greater_distances)

        variability_smaller_ag = variability_analysis(smaller_distances, mean_type="ag")
        variability_greater_ag = variability_analysis(greater_distances, mean_type="ag")

        p_value_smaller_greater = left_tailed_test_for_alternative_hypothesis(
            statistics_smaller=variability_smaller_ag,
            statistics_greater=variability_greater_ag,
            transformed=True,
        )
        print()
        print("-" * 80)
        print("COMPARATIVE VARIABILITY ANALYSIS:")
        print("-" * 80)
        print()
        print(
            f"""The mean distance among the {len(groups_smaller)} groups\n{groups_smaller}\n"""
            f"""is {variability_smaller_ag["mean"]} +- {variability_smaller_ag["margin_of_error"]} (margin of error) based on \n"""
            f"""with a standard deviation of {variability_smaller_ag["std"]}.\n"""
        )
        print(
            f"""The mean distance among the {len(groups_greater)} groups\n{groups_greater}\n"""
            f"""is {variability_greater_ag["mean"]} +- {variability_greater_ag["margin_of_error"]} (margin of error) based on \n"""
            f"""with a standard deviation of {variability_greater_ag["std"]}.\n"""
        )
        print(
            f"The hypothesis that the variability within\n{groups_smaller}\nis smaller than the variability within\n{groups_greater}\nhas a p-value of {p_value_smaller_greater}.\n"
        )
