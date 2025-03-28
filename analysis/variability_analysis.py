# SPDX-FileCopyrightText: 2025 Jakub W. Both <jakub.both@uib.no>
#
# SPDX-License-Identifier: MIT
from typing import Literal
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import logging
from pathlib import Path
from scipy.spatial.distance import squareform
from scipy.stats import t
from datastructure import SPECase
from utilities import reduce_distance_matrix_to_subset, nonlinear_transform


def variability_analysis(
    distance_matrix: np.ndarray,
    mean_type: Literal["arithmetic", "geometric", "ag"] = "ag",
):
    if isinstance(distance_matrix, list):
        distance_matrix = np.array(distance_matrix)

    if len(distance_matrix.shape) == 1:
        distance_values = distance_matrix.copy()
    elif len(distance_matrix.shape) == 2:
        distance_values = squareform(distance_matrix)

    # Transform values
    if mean_type == "arithmetic":
        transform_type = "linear"
    elif mean_type == "geometric":
        transform_type = "log"
    elif mean_type == "ag":
        transform_type = "symlog"
    else:
        raise ValueError("Mean type must be arithmetic, geometric, or ag.")

    transformed_values = nonlinear_transform(
        distance_values, transform_type=transform_type
    )

    # Plot distribution of distance_values_all
    if False:
        plt.figure("transformed")
        plt.title("Distribution of transformed distances for all submissions")
        plt.hist(transformed_values, bins=20)
        plt.figure("original")
        plt.title("Distribution of original distances for all submissions")
        plt.hist(distance_values, bins=20)
        plt.show()

    # Statistics on transformed values
    transformed_mean = np.mean(transformed_values)
    transformed_std = np.std(transformed_values)
    transformed_std_error = transformed_std / np.sqrt(len(transformed_values))
    z_value = 1.96  # for confidence_level = 0.95
    transformed_margin_of_error = transformed_std_error * z_value
    transformed_upper_limit_margin_of_error = (
        transformed_mean + transformed_margin_of_error
    )
    transformed_lower_limit_margin_of_error = (
        transformed_mean - transformed_margin_of_error
    )
    transformed_upper_limit_std = transformed_mean + transformed_std
    transformed_lower_limit_std = transformed_mean - transformed_std
    transformed_upper_limit_std_error = transformed_mean + transformed_std_error
    transformed_lower_limit_std_error = transformed_mean - transformed_std_error

    # Transform back
    mean = nonlinear_transform(
        transformed_mean, transform_type=transform_type, inverse=True
    )
    upper_limit_std_error = nonlinear_transform(
        transformed_upper_limit_std_error, transform_type=transform_type, inverse=True
    )
    lower_limit_std_error = nonlinear_transform(
        transformed_lower_limit_std_error, transform_type=transform_type, inverse=True
    )
    upper_limit_margin_of_error = nonlinear_transform(
        transformed_upper_limit_margin_of_error,
        transform_type=transform_type,
        inverse=True,
    )
    lower_limit_margin_of_error = nonlinear_transform(
        transformed_lower_limit_margin_of_error,
        transform_type=transform_type,
        inverse=True,
    )
    upper_limit_std = nonlinear_transform(
        transformed_upper_limit_std, transform_type=transform_type, inverse=True
    )
    lower_limit_std = nonlinear_transform(
        transformed_lower_limit_std, transform_type=transform_type, inverse=True
    )
    margin_of_error = max(
        upper_limit_margin_of_error - mean, mean - lower_limit_margin_of_error
    )
    std = max(upper_limit_std - mean, mean - lower_limit_std)
    std_error = max(upper_limit_std_error - mean, mean - lower_limit_std_error)

    return {
        # In original space
        "sample_size": len(distance_values),
        "mean": mean,
        "std": std,
        "std_error": std_error,
        "margin_of_error": margin_of_error,
        "upper_limit_margin_of_error": upper_limit_margin_of_error,
        "lower_limit_margin_of_error": lower_limit_margin_of_error,
        "upper_limit_std": upper_limit_std,
        "lower_limit_std": lower_limit_std,
        # In transformed space
        "transformed_mean": transformed_mean,
        "transformed_std": transformed_std,
        "transformed_std_error": transformed_std_error,
    }


def variability_analysis_discrete(
    distribution: np.ndarray,
    distance_matrix: np.ndarray,
    mean_type: Literal["arithmetic", "geometric", "ag"] = "ag",
):
    # Statistics on disctribution
    N = len(distribution)
    lower_5 = np.percentile(distribution, 5)
    upper_5 = np.percentile(distribution, 95)
    confidence_interval = upper_5 - lower_5
    print("Confidence interval: ", confidence_interval)
    print("Lower 5th percentile: ", lower_5)
    print("Upper 5th percentile: ", upper_5)

    if len(distance_matrix.shape) == 1:
        distance_values = distance_matrix.copy()
    elif len(distance_matrix.shape) == 2:
        distance_values = squareform(distance_matrix)

    median = np.median(distance_values)
    print("Median: ", median)

    # Transform values
    if mean_type == "arithmetic":
        transform_type = "linear"
    elif mean_type == "geometric":
        transform_type = "log"
    elif mean_type == "ag":
        transform_type = "symlog"
    else:
        raise ValueError("Mean type must be arithmetic, geometric, or ag.")

    transformed_values = nonlinear_transform(
        distance_values, transform_type=transform_type
    )

    # Plot distribution of distance_values_all
    if True:
        plt.figure("transformed")
        plt.title("Distribution of transformed distances for all submissions")
        plt.hist(transformed_values, bins=20)
        plt.figure("original")
        plt.title("Distribution of original distances for all submissions")
        plt.hist(distance_values, bins=20)
        plt.show()

    # Statistics on transformed values
    transformed_mean = np.mean(transformed_values)
    transformed_std = np.std(transformed_values)
    transformed_std_error = transformed_std / np.sqrt(len(transformed_values))
    z_value = 1.96  # for confidence_level = 0.95
    transformed_margin_of_error = transformed_std_error * z_value
    transformed_upper_limit_95 = transformed_mean + transformed_margin_of_error
    transformed_lower_limit_95 = transformed_mean - transformed_margin_of_error
    transformed_upper_limit_std = transformed_mean + transformed_std
    transformed_lower_limit_std = transformed_mean - transformed_std

    # Transform back
    mean = nonlinear_transform(
        transformed_mean, transform_type=transform_type, inverse=True
    )
    upper_limit_95 = nonlinear_transform(
        transformed_upper_limit_95, transform_type=transform_type, inverse=True
    )
    lower_limit_95 = nonlinear_transform(
        transformed_lower_limit_95, transform_type=transform_type, inverse=True
    )
    upper_limit_std = nonlinear_transform(
        transformed_upper_limit_std, transform_type=transform_type, inverse=True
    )
    lower_limit_std = nonlinear_transform(
        transformed_lower_limit_std, transform_type=transform_type, inverse=True
    )
    margin_of_error = max(upper_limit_95 - mean, mean - lower_limit_95)
    std = max(upper_limit_std - mean, mean - lower_limit_std)


def _cdf(distribution: np.ndarray, value: float, df: float) -> float:
    """Return the cumulative distribution function of value in a given distribution for certain degree of freedom."""
    return np.sum(distribution < value) / len(distribution)


def left_tailed_test_for_alternative_hypothesis(
    statistics_smaller: dict,
    statistics_greater: dict,
    transformed: bool = False,
    distribution: np.ndarray = None,
):
    """Return p-value for that statistics_smaller is smaller than statistics_greater.

    In other words, the null hypothesis is that the statistics_smaller is equal to the
    statistics_greater, and the alternative hypothesis is that the statistics_smaller is
    smaller than the statistics_greater.

    """

    # Fetch statistics
    if transformed:
        mean_smaller = statistics_smaller["transformed_mean"]
        std_error_smaller = statistics_smaller["transformed_std_error"]
        sample_size_smaller = statistics_smaller["sample_size"]
        mean_greater = statistics_greater["transformed_mean"]
        std_error_greater = statistics_greater["transformed_std_error"]
        sample_size_greater = statistics_greater["sample_size"]
    else:
        mean_smaller = statistics_smaller["mean"]
        std_error_smaller = statistics_smaller["std_error"]
        sample_size_smaller = statistics_smaller["sample_size"]
        mean_greater = statistics_greater["mean"]
        std_error_greater = statistics_greater["std_error"]
        sample_size_greater = statistics_greater["sample_size"]

    # T-score and degrees of freedom using Welch's t-test formula for whether
    # mean_smaller < mean_greater
    t_score = (mean_smaller - mean_greater) / np.sqrt(
        std_error_smaller**2 + std_error_greater**2
    )
    df = (std_error_smaller**2 + std_error_greater**2) ** 2 / (
        (std_error_smaller**2) ** 2 / (sample_size_smaller - 1)
        + (std_error_greater**2) ** 2 / (sample_size_greater - 1)
    )

    # Left-tailed test for p-value
    if distribution is None:
        p_value = t.cdf(t_score, df)
    else:
        assert False, "Not implemented yet."
        p_value = _cdf(distribution, t_score, df)
    return p_value


def two_tailed_test_for_alternative_hypothesis(
    statistics_smaller: dict,
    statistics_greater: dict,
    transformed: bool = False,
    distribution: np.ndarray = None,
):
    """Return p-value for that statistics_smaller is different from statistics_greater.

    In other words, the null hypothesis is that the statistics_smaller is equal to the
    statistics_greater, and the alternative hypothesis is that the statistics_smaller is
    different from the statistics_greater.

    """

    # Fetch statistics
    if transformed:
        mean_smaller = statistics_smaller["transformed_mean"]
        std_error_smaller = statistics_smaller["transformed_std_error"]
        sample_size_smaller = statistics_smaller["sample_size"]
        mean_greater = statistics_greater["transformed_mean"]
        std_error_greater = statistics_greater["transformed_std_error"]
        sample_size_greater = statistics_greater["sample_size"]
    else:
        mean_smaller = statistics_smaller["mean"]
        std_error_smaller = statistics_smaller["margin_of_error"]
        sample_size_smaller = statistics_smaller["sample_size"]
        mean_greater = statistics_greater["mean"]
        std_error_greater = statistics_greater["margin_of_error"]
        sample_size_greater = statistics_greater["sample_size"]

    # T-score and degrees of freedom using Welch's t-test formula for whether
    # mean_smaller < mean_greater
    t_score = (mean_smaller - mean_greater) / np.sqrt(
        std_error_smaller**2 + std_error_greater**2
    )
    df = (std_error_smaller**2 + std_error_greater**2) ** 2 / (
        (std_error_smaller**2) ** 2 / (sample_size_smaller - 1)
        + (std_error_greater**2) ** 2 / (sample_size_greater - 1)
    )

    # Two-tailed test for p-value
    if distribution is None:
        p_value = 2 * (1 - t.cdf(abs(t_score), df))
    else:
        raise ValueError("Not implemented yet.")
    return p_value


def cummulative_distribution_function(
    value: float,
    distribution: np.ndarray,
):
    """Return probability that a value is less or equal value, in a given distribution."""

    # Fetch statistics
    distribution_sorted = np.sort(distribution)

    # One-tailed test for p-value
    cdf = np.sum(distribution_sorted < value) / len(distribution_sorted)
    return cdf


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
            f"""The mean distance among the {num_smaller_groups} groups\n{groups_smaller}\n"""
            f"""is {variability_smaller_ag["mean"]} +- {variability_smaller_ag["margin_of_error"]} (margin of error) based on \n"""
            f"""with a standard deviation of {variability_smaller_ag["std"]}.\n"""
        )
        print(
            f"""The mean distance among the {num_greater_groups} groups\n{groups_greater}\n"""
            f"""is {variability_greater_ag["mean"]} +- {variability_greater_ag["margin_of_error"]} (margin of error) based on \n"""
            f"""with a standard deviation of {variability_greater_ag["std"]}.\n"""
        )
        print(
            f"The hypothesis that the variability within\n{groups_smaller}\nis smaller than the variability within\n{groups_greater}\nhas a p-value of {p_value_smaller_greater}.\n"
        )
