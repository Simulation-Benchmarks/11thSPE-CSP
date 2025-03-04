from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import squareform
from scipy.stats import t

from spe11_io import nonlinear_transform


def variability_analysis(
    distance_matrix: np.ndarray,
    mean_type: Literal["arithmetic", "geometric", "ag"] = "ag",
):
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
