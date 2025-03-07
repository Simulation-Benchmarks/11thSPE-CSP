"""Data classes for sparse and dense data."""

import numpy as np


class SPECase:
    def __init__(self, variant):
        """Initialize the SPE case."""

        self.variant = variant

        # Define the test case
        if self.variant == "spe11a":
            # Columns in the sparse reporting data which are relevant for the analysis
            self.data_format = {
                "t": 0,  # Time [s]
                "p1": 1,  # Pressure at observation point 1 [Pa]
                "p2": 2,  # Pressure at observation point 2 [Pa]
                "mobA": 3,  # Mobile free phase Box A [kg]
                "immA": 4,  # Immobile free phase Box A [kg]
                "dissA": 5,  # Dissolved CO2 in water Box A [kg]
                "sealA": 6,  # CO2 in seal Box A [kg]
                "mobB": 7,  # Mobile free phase Box B [kg]
                "immB": 8,  # Immobile free phase Box B [kg]
                "dissB": 9,  # Dissolved CO2 in water Box B [kg]
                "sealB": 10,  # CO2 in seal Box B [kg]
                "mC": 11,  # Convection measure [m]
                "sealTot": 12,  # Total CO2 in seal [kg]
            }

            self.spe11_distance_median = {
                "dissA": 0.07176848993149854,
                "dissB": 0.039408382947875735,
                "early_mass_w1": 0.0031770509108494235,
                "early_pressure_l2": 2962.2983042516043,
                "early_pressure_l2s": 366.1756663833846,
                "full_mass_w1": 0.01151651695441281,
                "full_pressure_l2": 6357.29909636252,
                "full_pressure_l2s": 2291.5109135015327,
                "immA": 0.0003004940803581391,
                "immB": 0.0008559326532330515,
                "late_mass_w1": 0.010015734349925,
                "late_pressure_l2": 5538.2541384088,
                "late_pressure_l2s": 2007.7063387277124,
                "mC": 117.09094179404104,
                "mobA": 0.07597671198450695,
                "mobB": 0.005546844025910249,
                "p1": 9971.718918358361,
                "p2": 607.5545021007509,
                "sealA": 0.0083075259191785,
                "sealB": 0.0035292158797797244,
                "sealTot": 0.028765405258409812,
                "t": 1.0,
            }

            # Isothermal
            self.non_isothermal = False

            # Reporting times - sparse data analysis
            self.time_unit = 3600  # seconds in an hour
            self.injection_stop = 5  # in hours
            self.reporting_time_unit = "h"
            self.reporting_times_sparse = (
                np.arange(0, 120 + 0.1, 0.1) * self.time_unit
            )  # in seconds
            self.reporting_times_dense = np.arange(0, 120 + 1, 1)[
                1:
            ]  # in hours - exclude 0 for integration

            # Special subgroups for the data
            self.subgroups = {
                # All submissions
                "all": sorted(
                    [
                        "calgary",
                        "cau-kiel",
                        "csiro1",
                        "csiro2",
                        "ctc-cne",
                        "geos1",
                        "geos2",
                        "ifpen",
                        "opengosim",
                        "opm1",
                        "opm2",
                        "opm3",
                        "opm4",
                        "pau-inria",
                        "pflotran",
                        "slb1",
                        "slb2",
                        "tetratech",
                        "ut-csee1",
                        "ut-csee2",
                        "ut-csee3",
                    ]
                ),
                # Standard reporting group - see Figure in the paper
                "base case": sorted(
                    [
                        "calgary",
                        "cau-kiel",
                        "csiro2",
                        "ctc-cne",
                        "geos1",
                        "ifpen",
                        "opengosim",
                        "opm1",
                        "pflotran",
                        "slb1",
                        "slb2",
                        "tetratech",
                        "ut-csee1",
                        "ut-csee2",
                        "ut-csee3",
                    ]
                ),
                # Standard reporting group - see Figure in the paper
                "base case non-commercial": sorted(
                    [
                        "calgary",
                        "cau-kiel",
                        "csiro2",
                        "geos1",
                        "ifpen",
                        "opengosim",
                        "opm1",
                        "pflotran",
                    ]
                ),
                # Standard reporting group - see Figure in the paper
                "base case commercial": sorted(
                    [
                        "ctc-cne",
                        "slb1",
                        "slb2",
                        "tetratech",
                        "ut-csee1",
                    ]
                ),
                # Different groups using the same simulators
                "opm": sorted(["cau-kiel", "opm1"]),
                "slb-IX": sorted(["slb2", "ctc-cne", "ut-csee1"]),
                # Meshes (all)
                "refined": sorted(["csiro1", "geos2", "opm4"]),
                "mesh_refinement_study_opm": sorted(["opm1", "opm4"]),
                "mesh_refinement_study_geos": sorted(["geos1", "geos2"]),
                # Facies adapted grids
                "facies adapted": sorted(
                    [
                        "pau-inria",
                        "opm3",
                    ]
                ),
            }

            self.results_and_categories = {
                "calgary": "base case",
                "cau-kiel": "base case",
                "csiro1": "refined",
                "csiro2": "base case",
                "ctc-cne": "base case",
                "geos1": "base case",
                "geos2": "refined",
                "ifpen": "base case",
                "opengosim": "base case",
                "opm1": "base case",
                "opm2": "special cases",
                "opm3": "facies adapted",
                "opm4": "refined",
                "pau-inria": "facies adapted",
                "pflotran": "base case",
                "slb1": "base case",
                "slb2": "base case",
                "tetratech": "base case",
                "ut-csee1": "base case",
                "ut-csee2": "base case",
                "ut-csee3": "base case",
            }

        elif variant == "spe11b":
            # Columns in the sparse reporting data which are relevant for the analysis
            self.data_format = {
                "t": 0,  # Time [s]
                "p1": 1,  # Pressure at observation point 1 [Pa]
                "p2": 2,  # Pressure at observation point 2 [Pa]
                "mobA": 3,  # Mobile free phase Box A [kg]
                "immA": 4,  # Immobile free phase Box A [kg]
                "dissA": 5,  # Dissolved CO2 in water Box A [kg]
                "sealA": 6,  # CO2 in seal Box A [kg]
                "mobB": 7,  # Mobile free phase Box B [kg]
                "immB": 8,  # Immobile free phase Box B [kg]
                "dissB": 9,  # Dissolved CO2 in water Box B [kg]
                "sealB": 10,  # CO2 in seal Box B [kg]
                "mC": 11,  # Convection measure [m]
                "sealTot": 12,  # Total CO2 in seal [kg]
                "boundaryCO2": 13,  # CO2 at the boundary [kg]
            }

            # Fix SPE11 distance median values based on official submissions
            self.spe11_distance_median = {
                "boundaryCO2": 604684682.3895249,
                "dissA": 53866193808.94402,
                "dissB": 29798204007.810722,
                "early_mass_w1": 4479583379725.279,
                "early_pressure_l2": 57833574123.78914,
                "early_pressure_l2s": 3790073212.2605157,
                "early_temperature_l2s": 40686.35307077792,
                "full_mass_w1": 20736931126443.348,
                "full_pressure_l2": 108603916995.56699,
                "full_pressure_l2s": 8011312689.846781,
                "full_temperature_l2s": 66258.07925225252,
                "immA": 2574286403.793894,
                "immB": 14216936478.372185,
                "late_mass_w1": 20250177682468.125,
                "late_pressure_l2": 92062835303.04178,
                "late_pressure_l2s": 5843305321.074707,
                "late_temperature_l2s": 50317.42351695232,
                "mC": 69710213.04146881,
                "mobA": 160975516363.48978,
                "mobB": 53803314384.14842,
                "p1": 167862794448.8785,
                "p2": 162638953414.91946,
                "sealA": 4298701516.103369,
                "sealB": 1854934823.843052,
                "sealTot": 17903125016.94225,
                "t": 1.0,
            }

            # Non-isothermal
            self.non_isothermal = True

            # Reporting times - sparse data analysis
            self.time_unit = 31536000
            self.injection_stop = 50
            self.reporting_time_unit = "y"
            self.reporting_times_sparse = (
                np.arange(0, 1000 + 0.1, 0.1) * self.time_unit
            )  # in seconds
            self.reporting_times_dense = np.arange(0, 1000 + 5, 5)[
                1:
            ]  # in years - exclude 0 - for integration

            # Special subgroups for the data
            self.subgroups = {
                # All submissions
                "all": sorted(
                    [
                        "calgary",
                        "cau-kiel",
                        "csiro",
                        "ctc-cne",
                        "darts",
                        "geos1",
                        "geos2",
                        "ifpen1",
                        "ifpen2",
                        "kfupm",
                        "opengosim1",
                        "opengosim2",
                        "opengosim3",
                        "opm1",
                        "opm2",
                        "opm3",
                        "opm4",
                        "pau-inria",
                        "pflotran",
                        "rice1",
                        "rice2",
                        "sintef1",
                        "sintef2",
                        "sintef3",
                        "sintef4",
                        "slb",
                        "stuttgart1",
                        "stuttgart2",
                        "stuttgart3",
                        "stuttgart4",
                        "tetratech1",
                        "tetratech2",
                        "ut-csee1",
                        "ut-csee2",
                    ]
                ),
                # Standard reporting group - see Figure in the paper
                "base case": sorted(
                    [
                        "pau-inria",
                        "stuttgart1",
                        "calgary",
                        "csiro",
                        "geos1",
                        "ifpen1",
                        "opengosim1",
                        "opm1",
                        "pflotran",
                        "tetratech1",
                        "ctc-cne",
                        "slb",
                        "ut-csee1",
                        "kfupm",
                        "tetratech2",
                        "ut-csee2",
                    ]
                ),
                # Standard reporting group - see Figure in the paper
                "base case non-commercial": sorted(
                    [
                        "pau-inria",
                        "stuttgart1",
                        "calgary",
                        "csiro",
                        "geos1",
                        "ifpen1",
                        "opengosim1",
                        "opm1",
                        "pflotran",
                    ]
                ),
                # Standard reporting group - see Figure in the paper
                "base case commercial": sorted(
                    [
                        "ctc-cne",
                        "kfupm",
                        "slb",
                        "tetratech1",
                        "tetratech2",
                        "ut-csee1",
                        "ut-csee2",
                    ]
                ),
                # Different simulators
                "tetratech": sorted(["tetratech1", "tetratech2"]),
                "ut-csee": sorted(["ut-csee1", "ut-csee2"]),
                # Different groups using the same simulators
                "dumux": sorted(["pau-inria", "stuttgart1"]),
                "slb-IX": sorted(["slb", "ctc-cne", "ut-csee1"]),
                # Meshes (all)
                "refined": sorted(
                    [
                        "cau-kiel",
                        "darts",
                        "geos2",
                        "ifpen2",
                        "opengosim2",
                        "opengosim3",
                        "opm4",
                    ]
                ),
                "mesh_refinement_study_opengosim": sorted(
                    ["opengosim1", "opengosim2", "opengosim3"]
                ),
                "mesh_refinement_study_sintef": sorted(["sintef1", "sintef3"]),
                "mesh_refinement_study_opm": sorted(["opm1", "opm4"]),
                "mesh_refinement_study_geos": sorted(["geos1", "geos2"]),
                "mesh_refinement_study_ifpen": sorted(["ifpen1", "ifpen2"]),
                # Facies adapted grids
                "facies adapted": sorted(
                    [
                        "opm2",
                        "sintef1",
                        "sintef2",
                        "sintef3",
                        "sintef4",
                    ]
                ),
            }
            self.subgroups["refined-vs-base case"] = sorted(
                self.subgroups["refined"] + self.subgroups["base case"]
            )
            self.subgroups["facies adapted-vs-base case"] = sorted(
                self.subgroups["facies adapted"] + self.subgroups["base case"]
            )
            self.subgroups["mesh_refinement_study"] = sorted(
                self.subgroups["mesh_refinement_study_opengosim"]
                + self.subgroups["mesh_refinement_study_opm"]
                + self.subgroups["mesh_refinement_study_geos"]
                + self.subgroups["mesh_refinement_study_ifpen"]
                + self.subgroups["mesh_refinement_study_sintef"]
            )
            self.subgroups["100k-cartesian-mesh"] = self.subgroups["base case"]
            self.subgroups["1.6m-cartesian-mesh"] = sorted(
                [
                    "ifpen2",
                    "geos2",
                    "opengosim3",
                    "opm4",
                ]
            )

            # Path to gmsh file for structured grid provided by the official SPE11 repository
            # It has been created using 'python .\make_structured_mesh.py --variant B -nx 840 -ny 120'
            # msh = Path("data/msh/spe11b_structured.msh")

            self.results_and_categories = {
                "calgary": "base case",
                "cau-kiel": "refined",
                "csiro": "base case",
                "ctc-cne": "base case",
                "darts": "refined",
                "geos1": "base case",
                "geos2": "refined",
                "ifpen1": "base case",
                "ifpen2": "refined",
                "kfupm": "base case",
                "opengosim1": "base case",
                "opengosim2": "refined",
                "opengosim3": "refined",
                "opm1": "base case",
                "opm2": "facies adapted",
                "opm3": "special cases",
                "opm4": "refined",
                "pau-inria": "base case",
                "pflotran": "base case",
                "rice1": "special cases",
                "rice2": "special cases",
                "sintef1": "facies adapted",
                "sintef2": "facies adapted",
                "sintef3": "facies adapted",
                "sintef4": "facies adapted",
                "slb": "base case",
                "stuttgart1": "base case",
                "stuttgart2": "special cases",
                "stuttgart3": "special cases",
                "stuttgart4": "special cases",
                "tetratech1": "base case",
                "tetratech2": "base case",
                "ut-csee1": "base case",
                "ut-csee2": "base case",
            }

            self.groups_and_mesh_size = {
                "sintef1": 25000,
                "sintef2": 25000,
                "sintef3": 100000,
                "sintef4": 100000,
                "ifpen1": 100000,
                "ifpen2": 1600000,
                "geos1": 100000,
                "geos2": 1600000,
                "opengosim1": 100000,
                "opengosim2": 200000,
                "opengosim3": 1600000,
                "opm1": 100000,
                "opm4": 10000000,
                "rice1": 10000,
                "rice2": 30000,
            }

        elif variant == "spe11c":
            # Columns in the sparse reporting data which are relevant for the analysis
            self.data_format = {
                "t": 0,  # Time [s]
                "p1": 1,  # Pressure at observation point 1 [Pa]
                "p2": 2,  # Pressure at observation point 2 [Pa]
                "mobA": 3,  # Mobile free phase Box A [kg]
                "immA": 4,  # Immobile free phase Box A [kg]
                "dissA": 5,  # Dissolved CO2 in water Box A [kg]
                "sealA": 6,  # CO2 in seal Box A [kg]
                "mobB": 7,  # Mobile free phase Box B [kg]
                "immB": 8,  # Immobile free phase Box B [kg]
                "dissB": 9,  # Dissolved CO2 in water Box B [kg]
                "sealB": 10,  # CO2 in seal Box B [kg]
                "mC": 11,  # Convection measure [m]
                "sealTot": 12,  # Total CO2 in seal [kg]
                "boundaryCO2": 13,  # CO2 at the boundary [kg]
            }

            self.spe11_distance_median = {
                "boundaryCO2": 47050470125889.16,
                "dissA": 130579678392805.52,
                "dissB": 13844297362.569794,
                "early_mass_w1": 3.466805485936395e17,
                "early_pressure_l2": 414141484905.823,
                "early_pressure_l2s": 126930592687.4999,
                "early_temperature_l2s": 2208546.5209046425,
                "full_mass_w1": 1.3269133295327048e18,
                "full_pressure_l2": 1044940693045.5431,
                "full_pressure_l2s": 223816392797.7573,
                "full_temperature_l2s": 3376825.498822822,
                "immA": 8341245970240.63,
                "immB": 78844943893.82198,
                "late_mass_w1": 1.2764058004583936e18,
                "late_pressure_l2": 943893526296.7959,
                "late_pressure_l2s": 185145149416.61313,
                "late_temperature_l2s": 2703567.0430007833,
                "mC": 124346821500.29715,
                "mobA": 226634096329277.88,
                "mobB": 391440466795.95404,
                "p1": 22577384735.74473,
                "p2": 21530076342.5297,
                "sealA": 13497337504810.178,
                "sealB": 7176169.471371533,
                "sealTot": 19572365609395.59,
                "t": 1.0,
            }

            # Non-isothermal
            self.non_isothermal = True

            # Reporting times
            self.time_unit = 31536000
            self.injection_stop = 50  # in years
            self.reporting_time_unit = "y"
            self.reporting_times_sparse = (
                np.arange(0, 1000 + 0.1, 0.1) * self.time_unit
            )  # in seconds
            self.reporting_times_dense = np.array(
                [
                    5,
                    10,
                    15,
                    20,
                    25,
                    30,
                    35,
                    40,
                    45,
                    50,
                    75,
                    100,
                    150,
                    200,
                    250,
                    300,
                    350,
                    400,
                    450,
                    500,
                    600,
                    700,
                    800,
                    900,
                    1000,
                ]
            )

            self.subgroups = {
                # All submissions
                "all": sorted(
                    [
                        "calgary",
                        "cau-kiel",
                        "csiro",
                        "ctc-cne",
                        "geos1",
                        "geos2",
                        "ifpen",
                        "opengosim1",
                        "opengosim2",
                        "opm1",
                        "opm2",
                        "opm3",
                        "opm4",
                        "pau-inria",
                        "pflotran",
                        "sintef1",
                        "sintef2",
                        "sintef3",
                        "slb",
                        "tetratech1",
                        "tetratech2",
                        "ut-csee",
                    ]
                ),
                # Standard reporting group - see Figure in the paper
                "base case": sorted(
                    [
                        "calgary",
                        "cau-kiel",
                        "csiro",
                        "ctc-cne",
                        "geos1",
                        "ifpen",
                        "opengosim1",
                        "opm1",
                        "pau-inria",
                        "slb",
                        "tetratech1",
                        "tetratech2",
                    ]
                ),
                # Standard reporting group - see Figure in the paper
                "base case commercial": sorted(
                    [
                        "cau-kiel",
                        "calgary",
                        "csiro",
                        "geos1",
                        "ifpen",
                        "opengosim1",
                        "opm1",
                        "pau-inria",
                    ]
                ),
                # Standard reporting group - see Figure in the paper
                "base case non-commercial": sorted(
                    [
                        "ctc-cne",
                        "slb",
                        "tetratech1",
                        "tetratech2",
                    ]
                ),
            }
            self.results_and_categories = {
                "calgary": "base case",
                "cau-kiel": "base case",
                "csiro": "base case",
                "ctc-cne": "base case",
                "geos1": "base case",
                "geos2": "refined",
                "ifpen": "base case",
                "opengosim1": "base case",
                "opengosim2": "refined",
                "opm1": "base case",
                "opm2": "facies adapted",
                "opm3": "special cases",
                "opm4": "refined",
                "pau-inria": "base case",
                "pflotran": "refined",
                "sintef1": "facies adapted",
                "sintef2": "facies adapted",
                "sintef3": "facies adapted",
                "slb": "base case",
                "tetratech1": "base case",
                "tetratech2": "base case",
                "ut-csee": "refined",
            }

            self.groups_and_mesh_size = {
                #    "sintef1": 25000,
                #    "sintef2": 25000,
                #    "sintef3": 100000,
                #    "sintef4": 100000,
                #    "ifpen1": 100000,
                #    "ifpen2": 1600000,
                #    "geos1": 100000,
                #    "geos2": 1600000,
                #    "opengosim1": 100000,
                #    "opengosim2": 200000,
                #    "opengosim3": 1600000,
                #    "opm1": 100000,
                #    "opm4": 10000000,
            }

        else:
            raise ValueError(f"Variant {variant} not known.")

        # ! ---- Common for all variants ----

        # Data formats
        self.data_format_from_index = {v: k for k, v in self.data_format.items()}

        # Timings - sparse (in seconds) and dense (in hours/years) data
        self.reporting_times_sparse_delta = np.insert(
            np.diff(self.reporting_times_sparse), 0, 0
        )
        self.time_weight_sparse = np.ones(len(self.reporting_times_sparse))
        self.time_weight_sparse[
            self.reporting_times_sparse > self.injection_stop * self.time_unit
        ] = 0.1
        self.reporting_times_dense_delta = np.insert(
            np.diff(self.reporting_times_dense), 0, 0
        )
        self.time_weight_dense = np.ones(len(self.reporting_times_dense))
        self.time_weight_dense[self.reporting_times_dense > self.injection_stop] = 0.1

        # Split the dense data time interval into early and late reporting times
        time_cut = np.searchsorted(self.reporting_times_dense, self.injection_stop)
        self.early_reporting_times_dense = self.reporting_times_dense[0 : time_cut + 1]
        self.late_reporting_times_dense = self.reporting_times_dense[time_cut:]
        self.early_reporting_duration = (
            self.reporting_times_dense[time_cut] - self.reporting_times_dense[0]
        )
        self.late_reporting_duration = (
            self.reporting_times_dense[-1] - self.reporting_times_dense[time_cut]
        )
        self.early_reporting_times_dense_delta = np.insert(
            np.diff(self.reporting_times_dense[0 : time_cut + 1]), 0, 0
        )
        self.late_reporting_times_dense_delta = np.diff(
            self.reporting_times_dense[time_cut:]
        )
        self.early_time_weight_dense = 1 * np.ones_like(
            self.early_reporting_times_dense_delta
        )
        self.late_time_weight_dense = 0.1 * np.ones_like(
            self.late_reporting_times_dense_delta
        )

        # Colors
        self.categories_and_colors = {
            "base case": "#91afdcff",
            "refined": "#77cf65ff",
            "facies adapted": "#ff9e54ff",
            "special cases": "#ffffffff",
        }

        # Colors from ../evaluation/groups_and_colors.py
        self.groups_and_colors = {
            "calgary": "#90B0DD",
            "cau-kiel": "#C40F10",
            "csiro": "#0B5197",
            "ctc-cne": "#FFA052",
            "darts": "#FF5A03",
            "geos": "#B0B0B0",
            "ifpen": "#000000",
            "kfupm": "#FF7573",
            "opengosim": "#AD92C3",
            "opm": "#7141A3",
            "pau-inria": "#AC7A71",
            "pflotran": "#683229",
            "rice": "#F39ABF",
            "sintef": "#D651A9",
            "slb": "#75D166",
            "stuttgart": "#127F12",
            "tetratech": "#07A4BB",
            "ut-csee": "#A1A30C",
        }

        # Find appropriate colors for text to be writtin on group_and_colors
        self.groups_and_text_colors = {
            "calgary": "black",
            "cau-kiel": "white",
            "csiro": "white",
            "ctc-cne": "black",
            "darts": "black",
            "geos": "black",
            "ifpen": "white",
            "kfupm": "black",
            "opengosim": "black",
            "opm": "white",
            "pau-inria": "black",
            "pflotran": "white",
            "rice": "black",
            "sintef": "black",
            "slb": "black",
            "stuttgart": "black",
            "tetratech": "black",
            "ut-csee": "black",
        }
