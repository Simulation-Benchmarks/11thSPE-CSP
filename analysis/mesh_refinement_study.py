"""Standardized workflow for data analysis - under development."""

import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from datastructure import SPECase
from utilities import (
    find_distance,
    reduce_distance_matrix_to_subset,
)
from variability_analysis import (
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
        required=True,
    )

    parser.add_argument(
        "-m",
        "--mesh-size",
        nargs="+",
        help="mesh sizes to be analyzed",
        action="append",
        required=True,
    )

    parser.add_argument(
        "-l",
        "--label",
        help="label for plotting",
        action="append",
    )

    parser.add_argument(
        "-variability-analysis",
        "--variability-analysis",
        help="perform variability analysis",
        nargs="+",
    )

    parser.add_argument(
        "-reference-group",
        "--reference-group",
        help="reference group for the analysis",
        nargs="+",
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
    groups = args.groups
    groups_lower = [[g.lower() for g in group] for group in groups]

    # Mesh size management
    mesh_sizes = args.mesh_size

    # Label management
    labels = args.label
    if labels is None:
        labels = np.arange(len(groups), dtype=int).tolist()

    # Further input for tuning the analysis
    _variability_analysis = args.variability_analysis
    reference_group = args.reference_group

    # Safety check
    assert len(groups) == len(mesh_sizes), (
        "Groups and mesh sizes must have the same length."
    )
    assert len(groups) == len(labels), "Groups and labels must have the same length."
    assert all(
        [len(group) == len(mesh_size) for group, mesh_size in zip(groups, mesh_sizes)]
    ), "Each group must have the same number of mesh sizes."
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

    def reduce_mesh_size(mesh_size):
        """Convert numbers to strings.

        Examples:
            21001 -> "21k"
            1200010 -> "1.2M"

        """
        mesh_size_str = []
        if not isinstance(mesh_size, list):
            mesh_size = [mesh_size]
        for size in mesh_size:
            if not isinstance(size, int):
                size = int(size)
            if size > 1_000_000:
                mesh_size_str.append(f"{size / 1_000_000:.1f}M")
            elif size > 1_000:
                mesh_size_str.append(f"{size / 1_000:.0f}k")
            else:
                mesh_size_str.append(f"{size}")
        return mesh_size_str

    # Store mesh refinement factors for consecutive mesh sizes within the same group
    mesh_refinement_factor = {}
    mesh_refinement_info = {}
    mesh_refinement_distance_values = {}
    counter = 0
    for group, mesh_size in zip(groups, mesh_sizes):
        mrf = []
        mri = []
        dv = []
        # Consider all increasing pairs of mesh sizes
        for i in range(len(mesh_size) - 1):
            for j in range(i + 1, len(mesh_size)):
                mrf.append(np.sqrt(float(mesh_size[j]) / float(mesh_size[i])))
                mri.append(
                    f"{reduce_mesh_size(mesh_size)[i]} vs. {reduce_mesh_size(mesh_size)[j]}"
                )
                dv.append(
                    find_distance(distance_matrix, (group[i], group[j]), all_groups)
                )
        mesh_refinement_factor[counter] = mrf
        mesh_refinement_info[counter] = mri
        mesh_refinement_distance_values[counter] = dv
        counter += 1

    plt.figure("Mesh refinement study")
    fig = plt.gcf()
    fig.set_size_inches(14, 14)
    ax = plt.gca()
    for key in range(counter):
        try:
            color = spe_case.groups_and_colors[labels[key]]
        except KeyError:
            # Choose a color from the default color cycle
            color = plt.gca()._get_lines.get_next_color()
        sct = plt.scatter(
            mesh_refinement_factor[key],
            np.array(mesh_refinement_distance_values[key]),
            marker="o",
            s=100,
            label=labels[key],
            color=color,
            edgecolors="gray",
            zorder=20,
        )
        for i, txt in enumerate(mesh_refinement_info[key]):
            try:
                plt.text(
                    mesh_refinement_factor[key][i] + 0.15,
                    mesh_refinement_distance_values[key][i] - 0.025,
                    txt,
                    fontsize=14,
                    color=[0.3] * 3,
                    zorder=30,
                )
            except IndexError:
                plt.text(
                    mesh_refinement_factor[key][i] + 0.15,
                    mesh_refinement_distance_values[key] - 0.025,
                    txt,
                    fontsize=14,
                    color=[0.3] * 3,
                    zorder=30,
                )

    # Include variability analysis
    if _variability_analysis:
        distance_values = []
        for key in _variability_analysis:
            counter = labels.index(key)
            distance_values.append(mesh_refinement_distance_values[counter])
        distance_values = np.hstack(distance_values).flatten()

        variability_refinement = variability_analysis(
            distance_values,
            mean_type="ag",
        )

        # Add mean value line
        plt.axhline(
            y=variability_refinement["mean"],
            color="black",
            label="Mean variability",
        )
        # Add margin of error
        plt.axhline(
            y=variability_refinement["lower_limit_margin_of_error"],
            color="gray",
            # linestyle="--",
        )
        plt.axhline(
            y=variability_refinement["upper_limit_margin_of_error"],
            color="gray",
            # linestyle="--",
        )
        plt.fill_between(
            [0, 10.5],
            variability_refinement["lower_limit_margin_of_error"],
            variability_refinement["upper_limit_margin_of_error"],
            color=[0.9] * 3,
        )

    if reference_group:
        reference_distance_matrix = reduce_distance_matrix_to_subset(
            distance_matrix, reference_group, all_groups
        )
        np.fill_diagonal(reference_distance_matrix, 0)
        reference_variability = variability_analysis(
            reference_distance_matrix,
            mean_type="ag",
        )
        # Add the same for the SRG
        plt.axhline(
            y=reference_variability["mean"],
            color="black",
            linestyle="--",
            label="Reference variability",
        )
        plt.axhline(
            y=reference_variability["lower_limit_margin_of_error"],
            color="gray",
            linestyle="--",
        )
        plt.axhline(
            y=reference_variability["upper_limit_margin_of_error"],
            color="gray",
            linestyle="--",
        )
        plt.fill_between(
            [0, 10.5],
            reference_variability["lower_limit_margin_of_error"],
            reference_variability["upper_limit_margin_of_error"],
            color=[0.9] * 3,
        )

    plt.legend(loc="upper left", fontsize=14)
    plt.xlabel("Mesh refinement factor", fontsize=14)
    # Add grid behind the scatter plot
    plt.grid(zorder=10, color="lightgray")
    plt.xlim(0.5, 10.5)
    plt.ylabel("SPE11 distance", fontsize=14)
    # Set fontsize large
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)

    plt.tight_layout()
    plt.savefig(save_folder / f"{spe_case.variant}_mesh_refinement_study.png", dpi=1000)
    plt.show()
    print(
        f"Mesh refinement study plot saved under {save_folder / spe_case.variant}_mesh_refinement_study.png."
    )
