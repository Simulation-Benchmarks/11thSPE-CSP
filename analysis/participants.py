"""Fixed metadata for the official participants of SPE11."""


def get_result_name_from_multiple(folder):

    # folder has the structur .../team/speID/result<N> and is a Path object. I need to split these three parts
    team = Path(folder).parts[-3]
    speID = Path(folder).parts[-2]
    result = Path(folder).parts[-1]

    # Reduce to identifiers
    result = result.split("result")[-1]
    speID = speID[-1]

    return f"{speID}_{team}{result}"


def get_result_name_from_unique(folder):

    # folder has the structur .../team/speID/result<N> and is a Path object. I need to split these three parts
    team = Path(folder).parts[-2]
    speID = Path(folder).parts[-1]

    # Reduce to identifiers
    speID = speID[-1]

    return f"{speID}_{team}"


def get_std_name(name):
    """Return the standard name of the participant.

    Parameters
    ----------
    name : str
        The name of the participant.

    Returns
    -------
    str
        The standard name of the participant.

    """
    conversion = {
        "calgary": "calgary",
        "csiro": "csiro",
        "ctc-cne": "ctc-cne",
        "delft-darts": "darts",
        "geos": "geos",
        "ifpen": "ifpen",
        "kfupm": "kfupm",
        "kiel": "kiel",
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
    if name not in conversion:
        raise ValueError(f"Unknown participant '{name}'.")
    return conversion[name]
