# SPDX-FileCopyrightText: 2025 Jakub W. Both <jakub.both@uib.no>
#
# SPDX-License-Identifier: MIT
"""Compute statistical significance through different correlation tests."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import kendalltau


def plot_correlation_test(data1, data2, key1, key2):

    # Pearson Correlation
    p_corr = pearson_correlation(data1, data2)

    # Linear regression
    slope, intercept = np.polyfit(np.ravel(data1), np.ravel(data2), 1)

    # Residual
    residual = np.ravel(data2) - (slope * np.ravel(data1) + intercept)
    print(
        p_corr,
        "res",
        np.sqrt(np.mean(residual**2)),
        np.sqrt(np.mean(residual**2)) / (np.max(data2) - np.min(data2)),
    )

    plt.plot(np.ravel(data1), np.ravel(data2), "o", color="blue", label="Data")
    plt.plot(
        np.ravel(data1),
        slope * np.ravel(data1) + intercept,
        color="red",
        label="y = %.2fx + %.2f" % (slope, intercept),
    )
    plt.text(0.5, 1, f"Pearson correlation: {p_corr:.2f}", fontsize=12)
    plt.xlabel(f"Metric ({key1})")
    plt.ylabel(f"Metric ({key2})")
    plt.title("Correlation Test")
    plt.legend()
    Path("sensitivity_analysis").mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(f"sensitivity_analysis/correlation_{key1}_{key2}.png")
    plt.close()
    # plt.show()


def pearson_correlation(data1, data2):
    """Compute the Pearson correlation coefficient between two datasets."""
    # print(np.corrcoef(np.ravel(data1), np.ravel(data2)).shape)
    # print(np.corrcoef(np.ravel(data1), np.ravel(data2))[0, 1])
    # print(np.corrcoef(np.vstack((np.ravel(data1), np.ravel(data2)))))
    # cov = np.cov(data1, data2)
    # std1 = np.std(data1)
    # std2 = np.std(data2)
    # print(cov / (std1 * std2))
    # assert False
    return np.corrcoef(np.ravel(data1), np.ravel(data2))[0, 1]


def pearson_correlation_analysis(data_dict, ref_data) -> dict:
    """Compute the correlation between the reference data and the other datasets."""
    correlation_dict = {}
    for key, data in data_dict.items():
        correlation_dict[key] = pearson_correlation(data, ref_data)
    return correlation_dict


def pearson_cross_correlation_analysis(data):

    correlation_dict = {}
    correlation_matrix = np.zeros((len(data), len(data)))
    for counter1, (key1, data1) in enumerate(data.items()):
        for counter2, (key2, data2) in enumerate(data.items()):
            correlation_dict[f"{key1}_{key2}"] = pearson_correlation(data1, data2)
            correlation_matrix[counter1, counter2] = pearson_correlation(data1, data2)

    # Put on the diag the mean value of the correlation in the respective row
    for i in range(len(data)):
        correlation_matrix[i, i] = np.mean(correlation_matrix[i, :])

    def crop_cmap(cmap, min_val=0.0, max_val=1.0):
        cmap = plt.get_cmap(cmap)
        colors = cmap(np.linspace(min_val, max_val, cmap.N))
        return LinearSegmentedColormap.from_list("cropped_cmap", colors)

    cropped_cmap = crop_cmap("binary", 0.0, 0.4)  # Using 20% of the 'viridis' colormap

    plt.figure(figsize=(14, 14))
    ax = plt.gca()
    cax = ax.matshow(correlation_matrix, cmap=cropped_cmap)
    # Add colorbar with label 'Pearson Correlation Coefficient'
    cbar = plt.colorbar(cax)
    cbar.set_label("Pearson Correlation Coefficient")

    plt.xticks(range(len(data)), data.keys(), rotation=90)
    plt.yticks(range(len(data)), data.keys())
    for i in range(len(data)):
        for j in range(len(data)):
            plt.text(
                j,
                i,
                f"{correlation_matrix[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=8,
            )
    plt.show()
    return correlation_dict


def kendall_tau_correlation_analysis(data_dict, ref_data) -> dict:
    """Compute the correlation between the reference data and the other datasets."""
    correlation_dict = {}
    for key, data in data_dict.items():
        correlation_dict[key] = kendall_tau_correlation(data, ref_data)
    return correlation_dict


def kendall_tau_correlation(data1, data2):
    """Compute the Kendall Tau correlation coefficient between two datasets."""
    return kendalltau(np.ravel(data1), np.ravel(data2))


def kendall_tau_distance(matrix_data1, matrix_data2):
    data1 = np.ravel(matrix_data1)
    data2 = np.ravel(matrix_data2)
    assert len(data1) == len(data2), "Data1 and Data2 must have the same length"
    n = len(data1)
    num_pairs = n * (n - 1) // 2
    kendall_indices = [
        (i, j)
        for i in range(n)
        for j in range(i + 1, n)
        if (data1[i] - data1[j]) * (data2[i] - data2[j]) < 0
    ]
    return len(kendall_indices) / num_pairs


def kendall_tau_cross_correlation_analysis(data):

    correlation_dict = {}
    correlation_matrix = np.zeros((len(data), len(data)))
    for counter1, (key1, data1) in enumerate(data.items()):
        for counter2, (key2, data2) in enumerate(data.items()):
            distance = kendall_tau_distance(data1, data2)
            correlation_dict[f"{key1}_{key2}"] = distance
            correlation_matrix[counter1, counter2] = distance

    # Put on the diag the mean value of the correlation in the respective row
    for i in range(len(data)):
        correlation_matrix[i, i] = np.mean(correlation_matrix[i, :])

    def crop_cmap(cmap, min_val=0.0, max_val=1.0):
        cmap = plt.get_cmap(cmap)
        colors = cmap(np.linspace(min_val, max_val, cmap.N))
        return LinearSegmentedColormap.from_list("cropped_cmap", colors)

    cropped_cmap = crop_cmap("binary", 0.0, 0.4)  # Using 20% of the 'viridis' colormap

    plt.figure(figsize=(14, 14))
    ax = plt.gca()
    cax = ax.matshow(correlation_matrix, cmap=cropped_cmap)
    for i in range(len(data)):
        for j in range(len(data)):
            plt.text(
                j,
                i,
                f"{correlation_matrix[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=8,
            )
    cbar = plt.colorbar(cax)
    cbar.set_label("Kendall tau distance")
    plt.xticks(range(len(data)), data.keys(), rotation=90)
    plt.yticks(range(len(data)), data.keys())
    plt.show()
    return correlation_dict


def kendall_tau_cross_distance_analysis(data):

    distance_dict = {}
    distance_matrix = np.zeros((len(data), len(data)))
    for counter1, (key1, data1) in enumerate(data.items()):
        for counter2, (key2, data2) in enumerate(data.items()):
            distance = kendall_tau_distance(data1, data2)
            distance_dict[f"{key1}_{key2}"] = distance
            distance_matrix[counter1, counter2] = distance

    # Put on the diag the mean value of the distance in the respective row
    for i in range(len(data)):
        distance_matrix[i, i] = np.mean(distance_matrix[i, :])

    def crop_cmap(cmap, min_val=0.0, max_val=1.0):
        cmap = plt.get_cmap(cmap)
        colors = cmap(np.linspace(min_val, max_val, cmap.N))
        return LinearSegmentedColormap.from_list("cropped_cmap", colors)

    cropped_cmap = crop_cmap("binary", 0.0, 0.4)  # Using 20% of the 'viridis' colormap

    plt.figure(figsize=(14, 14))
    ax = plt.gca()
    cax = ax.matshow(distance_matrix, cmap=cropped_cmap)
    for i in range(len(data)):
        for j in range(len(data)):
            plt.text(
                j,
                i,
                f"{distance_matrix[i, j]:.2f}",
                ha="center",
                va="center",
                fontsize=8,
            )
    cbar = plt.colorbar(cax)
    cbar.set_label("Kendall tau distance")
    plt.xticks(range(len(data)), data.keys(), rotation=90)
    plt.yticks(range(len(data)), data.keys())
    plt.show()
    return distance_dict
