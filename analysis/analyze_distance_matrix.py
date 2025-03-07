"""Standardized workflow for data analysis - under development."""

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from cluster_analysis import (
    argmin_distance,
    centroid_analysis,
    determine_median_cluster,
    mean_distance_to_group,
    std_distance_to_group,
    plot_linkage_clustering_with_colors,
)
from datastructure import SPECase
from utilities import (
    reduce_distance_matrix_to_subset,
    plot_heatmap_distance_matrix,
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
        action="append",
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

    # Verbosity
    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    # ! ---- READ SPARSE DATA ----

    # Read distance matrix using pandas - note it contains a header and index
    all_distance_matrix_with_labels = pd.read_csv(
        save_folder / f"{spe_case.variant}_distance_matrix.csv",
        delimiter=",",
        index_col=0,
    )
    all_distance_matrix = all_distance_matrix_with_labels.to_numpy()
    all_groups = all_distance_matrix_with_labels.columns.to_numpy().tolist()
    all_groups_lower = [g.lower() for g in all_groups]

    # Group management
    if args.groups:
        groups = args.groups
    else:
        groups = [all_groups]
    groups_lower = [[g.lower() for g in group] for group in groups]
    groups_flat = [item for sublist in groups for item in sublist]
    groups_flat_lower = [g.lower() for g in groups_flat]

    match args.option:
        case "print-distances":
            # Collect all pairwise combinations (without self-comparisons and duplicates)
            for group in groups:
                for i, key0 in enumerate(group):
                    for j, key1 in enumerate(group):
                        if i < j:
                            ind0 = all_groups_lower.index(key0.lower())
                            ind1 = all_groups_lower.index(key1.lower())
                            print(
                                f"Distance between {key0} and {key1}: {all_distance_matrix[ind0, ind1]}"
                            )

        case "show-distance-matrix":
            assert len(groups) == 1, "Only one group can be selected for this analysis."

            # Restrict to the groups of interest
            distance_matrix = reduce_distance_matrix_to_subset(
                all_distance_matrix, groups_flat, all_groups
            )
            plot_heatmap_distance_matrix(
                distance_matrix,
                groups_flat,
                spe_case,
                path=save_folder / f"{spe_case.variant}_heatmap_distance_matrix.png",
            )

        case "show-clustering":
            # Restrict to the groups of interest
            distance_matrix = reduce_distance_matrix_to_subset(
                all_distance_matrix, groups_flat, all_groups
            )

            # Need to make sure the diagonal is zero as preparation for clustering
            np.fill_diagonal(distance_matrix, 0)

            # Visual inspection in form of a dendrogram
            plot_linkage_clustering_with_colors(
                spe_case,
                distance_matrix,
                groups_flat,
                path=save_folder / f"{spe_case.variant}_dendrogram",
            )

        case "find-min-mean-distance":
            # ! ---- MINIMUM AG MEAN DISTANCE  ----

            assert len(groups) == 1, "Only one group can be selected for this analysis."

            # Restrict to the groups of interest
            distance_matrix = reduce_distance_matrix_to_subset(
                all_distance_matrix, groups[0], all_groups
            )

            mean_distance_all = {
                key: mean_distance_to_group(
                    distance_matrix,
                    key.lower(),
                    groups_lower[0],
                    mean_type="ag",
                )
                for key in all_groups
            }
            std_distance_all = {
                key: std_distance_to_group(
                    distance_matrix,
                    key.lower(),
                    groups_lower[0],
                    mean_type="ag",
                )
                for key in all_groups
            }

            print()
            print("-" * 80)
            print("MINIMUM AG MEAN DISTANCE:")
            print("-" * 80)
            print()
            print(
                f"""Within the group:\n{groups}\n"""
                f"""With distance {mean_distance_all[argmin_distance(mean_distance_all)]} and """
                f"""std {std_distance_all[argmin_distance(mean_distance_all)]}.\n"""
                """The group with the smallest mean distance to all other participants is:\n"""
                f"""{argmin_distance(mean_distance_all)}\n"""
            )

        case "median-centroid-analysis":
            assert len(groups) == 1, "Only one group can be selected for this analysis."

            # Restrict to the groups of interest
            distance_matrix = reduce_distance_matrix_to_subset(
                all_distance_matrix, groups[0], all_groups
            )

            # Quantitative inspection of the "median" data
            median_cluster_all = determine_median_cluster(
                distance_matrix,
                groups_lower[0],
                mean_type="ag",
            )

            # Determine the mean distances of the two representatives of the median cluster to the other participants
            mean_distance_cluster_all = {
                key: mean_distance_to_group(
                    distance_matrix,
                    key,
                    groups[0],
                    mean_type="ag",
                )
                for key in median_cluster_all
            }
            print(
                f"The median cluster within all submissions consists of {median_cluster_all}, and {argmin_distance(mean_distance_cluster_all)} is the closest to all other participants."
            )

            centroid_analysis(
                distance_matrix,
                groups[0],
                mean_type="ag",
            )

        case "find-min-distance":
            # Extract the full distance matrix for all groups
            distance_matrix = reduce_distance_matrix_to_subset(
                all_distance_matrix, groups_flat, all_groups
            )

            # If multiple groups are provided, remove group-intern distances,
            # otherwise only exclude self-comparisons.
            if len(groups) > 1:
                for group in groups:
                    for key0 in group:
                        for key1 in group:
                            ind0 = groups_flat.index(key0)
                            ind1 = groups_flat.index(key1)
                            distance_matrix[ind0, ind1] = np.inf
            else:
                np.fill_diagonal(distance_matrix, np.inf)

            # Find the minimum distance and the two participants with the smallest distance
            min_distance = np.min(distance_matrix)
            min_indices = np.unravel_index(
                np.argmin(distance_matrix), distance_matrix.shape
            )
            min_participants = (
                groups_flat[min_indices[0]],
                groups_flat[min_indices[1]],
            )

            print()
            print("-" * 80)
            print("MINIMUM DISTANCE:")
            print("-" * 80)
            print()
            print(
                f"The minimum distance is {min_distance} between {min_participants[0]} and {min_participants[1]}."
            )
