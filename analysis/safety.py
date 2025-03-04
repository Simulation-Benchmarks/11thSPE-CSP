"""Sanity checks."""

import logging
from warnings import warn

import numpy as np

from datastructure import SPECase


def standard_check(data, key: str):
    if np.any(np.isclose(data, -1)):
        logging.info(f"Missing values detected for submission {key}.")
        return False
    if np.any(np.isclose(data, -999)):
        logging.info(f"Missing values detected for submission {key}.")
        return False
    if np.isnan(data).any():
        logging.info(f"NaN values detected for submission {key}.")
        return False
    if np.any(data < 0):
        logging.info(
            f"Negative values (min value {np.min(data)}) detected for submission {key}."
        )
    if np.any(data == np.inf):
        logging.info(f"Positive infinity values detected for submission {key}.")
        return False
    if np.any(data == -np.inf):
        logging.info(f"Negative infinity values detected for submission {key}.")
        return False
    return True


def check_M_values(data, speCase: SPECase, key: str):
    m_values = data[:, speCase.data_format["M_C"]]
    return standard_check(m_values, key)


def check_boundaryCO2_values(data, speCase: SPECase, key: str):
    if speCase.variant == "spe11a":
        # No boundary CO2 values for variant A
        return True
    boundaryCO2_values = data[:, speCase.data_format["boundaryCO2"]]
    return standard_check(boundaryCO2_values, key)


def check_sanity(data, key: str):
    return standard_check(data, key)


def clean_data(data, replace_nan=False):
    """Clean data using hardcoded values based on usage by SPE11 participants."""

    clean_data = data.copy()

    if np.any(np.isclose(clean_data, -999)):
        clean_data[np.isclose(clean_data, -999)] = 0

    if replace_nan and np.any(np.isnan(clean_data)):
        clean_data[np.isnan(clean_data)] = 0

    return clean_data


def clean_names(names):
    warn("This function should be removed after the data is updated.")
    clean_names = []

    for name in names:
        # Clean name if something like 'b_NAME2' which follows a preliminary naming convention
        if "_" in name:
            count = name.count("_")
            assert count == 1, f"Name {names} has {count} underscores."
            assert len(name.split("_")[0]) == 1, (
                f"Name {name} has more than one character before the underscore."
            )
            name = name.split("_")[1]

        # Check if name ends with digit
        if name[-1].isdigit():
            digit = name[-1]
            name = name[:-1]
        else:
            digit = None

        # Clean name if among previous group names
        conversion = {
            "calgary": "calgary",
            "csiro": "csiro",
            "ctc-cne": "ctc-cne",
            "delft-darts": "darts",
            "geos": "geos",
            "ifpen": "ifpen",
            "kfupm": "kfupm",
            "kiel": "cau-kiel",
            "opengosim": "opengosim",
            "opm": "opm",
            "pau-inria": "pau-inria",
            "pnnl": "pflotran",
            "rice": "rice",
            "sintef": "sintef",
            "slb": "slb",
            "stuttgart": "stuttgart",
            "tetratech-rps": "tetratech",
            "ut-csee-pge": "ut-csee",
        }

        if name in conversion:
            name = conversion[name]

        # Add digit back
        if digit:
            name = f"{name}{digit}"

        clean_names.append(name)

    return clean_names
