"""
pyrad.io.config
===============

Functions for reading pyrad config files

.. autosummary::
    :toctree: generated/

    read_config
    get_num_elements
    string_to_datatype
    get_array
    get_struct
    get_array_type
    init_array

"""
import re
import os
import yaml

# Allows parsing of environment variables in YAML safe loader
env_var_pattern = re.compile(r"\$\{([^}^{]+)\}")


def env_var_constructor(loader, node):
    value = loader.construct_scalar(node)
    return os.path.expandvars(value)


yaml.SafeLoader.add_implicit_resolver("!env_var", env_var_pattern, None)
yaml.SafeLoader.add_constructor("!env_var", env_var_constructor)

# Define type mappings
SCALAR_TYPES = {
    "BYTE": int,
    "BOOL": lambda x: str(x).lower() in ("true", "1"),
    "INT": int,
    "LONG": int,
    "HEX": lambda x: int(x, 16)
    if isinstance(x, str) and x.startswith("0x")
    else int(x),
    "EXP": float,
    "FLOAT": float,
    "DOUBLE": float,
    "STRING": str,
}

ARRAY_TYPE_NAMES = {
    "BYTARR": "BYTE",
    "INTARR": "INT",
    "LONARR": "LONG",
    "HEXARR": "HEX",
    "EXPARR": "EXP",
    "FLTARR": "FLOAT",
    "DBLARR": "DOUBLE",
    "STRARR": "STRING",
}


def read_config(fname, cfg=None, defaults=None):
    """
    Read a pyrad config file.
    It can use either the classical pyrad config syntax or yaml files

    Parameters
    ----------
    fname : str
        Name of the configuration file to read.

    cfg : dict of dicts, optional
        dictionary of dictionaries containing configuration parameters where
        the new parameters will be placed

    defaults: dict of dicts, optional
        dictionary of dictionaries containing default values. If a key is
        contained in the defaults dict but not in the configuration file,
        the value from the defaults dict will be assigned

    Returns
    -------
    cfg : dict of dicts
        dictionary of dictionaries containing the configuration parameters

    """
    # if config dictionary does not exist yet create it
    if cfg is None:
        cfg = dict()

    # check if the file can be read
    try:
        with open(fname, "r", encoding="utf-8", errors="ignore") as cfgfile:
            # Figure out if it is yaml or not
            if fname.endswith(".yaml") or fname.endswith(".yml"):
                cfg_new = yaml.load(cfgfile, Loader=yaml.SafeLoader)
            else:
                cfg_new = _read_config_pyrad(cfgfile)
    except Exception:
        raise
        raise Exception("ERROR: Could not find|open config file '" + fname + "'")

    # Merge with existing dict
    cfg = merge_dicts(cfg, cfg_new)

    # if default does not exist, create it
    if defaults is None:
        defaults = dict()

    # Verify that all keys in default are in newly created config
    cfg = merge_dicts(cfg, defaults)
    return cfg


def _read_config_pyrad(cfgfile):
    lines = [
        os.path.expandvars(line).split("#")[0].strip()
        for line in cfgfile.readlines()
        if line.strip() and not line.strip().startswith("#")
    ]

    def parse_block(line, index, expected_count):
        result = {}
        count = 0
        while index < len(lines):
            line = lines[index]
            parts = line.split(maxsplit=2)

            if len(parts) < 2:
                raise ValueError(f"Invalid line at {index + 1}: '{line}'")

            key = parts[0]
            type_str = parts[1].upper()

            if type_str == "STRUCT":
                if len(parts) != 3:
                    raise ValueError(
                        f"Missing struct count at line {index + 1}: '{line}'"
                    )
                struct_count = int(parts[2])
                index += 1
                struct_data, index = parse_block(line, index, struct_count)
                result[key] = struct_data
                count += 1

            elif type_str in ARRAY_TYPE_NAMES:
                if len(parts) != 3:
                    raise ValueError(
                        f"Missing array count at line {index + 1}: '{line}'"
                    )
                arr_count = int(parts[2])
                index += 1
                arr = []
                for _ in range(arr_count):
                    if index >= len(lines):
                        raise ValueError(
                            f"Expected {arr_count} array items for '{key}'"
                        )
                    arr.append(lines[index])
                    index += 1
                result[key] = arr
                count += 1

            elif type_str in SCALAR_TYPES:
                if len(parts) != 3:
                    raise ValueError(
                        f"Expected value for key '{key}' at line {index + 1}"
                    )
                raw_val = parts[2]
                parser = SCALAR_TYPES[type_str]
                value = parser(raw_val)
                result[key] = value
                index += 1
                count += 1

            else:
                break  # End of current block

            if expected_count is not None and count == expected_count:
                break

        if expected_count is not None and count != expected_count:
            raise ValueError(
                f"At line {line}, expected {expected_count} elements, but parsed {count}"
            )
        return result, index

    config = {}
    index = 0
    while index < len(lines):
        line = lines[index]
        parts = line.split(maxsplit=2)
        if len(parts) < 2:
            raise ValueError(
                f"Invalid top-level declaration at line {index+1}: '{line}'"
            )
        key, type_str = parts[0], parts[1].upper()

        if type_str == "STRUCT":
            if len(parts) != 3:
                raise ValueError(f"Missing struct count at line {index + 1}: '{line}'")
            struct_count = int(parts[2])
            index += 1
            struct_data, index = parse_block(line, index, struct_count)
            config[key] = struct_data

        elif type_str in ARRAY_TYPE_NAMES:
            if len(parts) != 3:
                raise ValueError(f"Missing array count at line {index + 1}: '{line}'")
            arr_count = int(parts[2])
            index += 1
            arr = []
            for _ in range(arr_count):
                if index >= len(lines):
                    raise ValueError(f"Expected {arr_count} array items for '{key}'")
                if lines[index] == "\n" or "STRUCT" in lines[index]:
                    raise ValueError(
                        f"Parsed illegal line while reading array for {line}"
                    )
                arr.append(lines[index])
                index += 1
            config[key] = arr

        elif type_str in SCALAR_TYPES:
            if len(parts) != 3:
                raise ValueError(f"Expected value for key '{key}' at line {index + 1}")
            raw_val = parts[2]
            parser = SCALAR_TYPES[type_str]
            value = parser(raw_val)
            config[key] = value
            index += 1
        else:
            raise ValueError(f"Unknown type '{type_str}' at line {index+1}")

    if "dataSetList" in config:
        if len(config["dataSetList"]):
            # Check if all elements are in config
            for dset in config["dataSetList"]:
                if ":" in dset:
                    dset = dset.split(":")[1]
                if dset not in config:
                    raise ValueError(
                        f"Dataset {dset} specificied in dataSetList was not found in config file"
                    )
    return config


def merge_dicts(ref, defaults):
    """
    Merge two nested dictionaries recursively.

    Parameters:
    -----------
    ref : dict
        The base dictionary to be merged into.
    defaults : dict
        The dictionary whose keys and values are merged into dictionary ref.

    Returns:
    --------
    ref: dict
        The base dictionary after addition of missing keys from defaults
    """
    for key, value in defaults.items():
        if key in ref:
            if isinstance(ref[key], dict) and isinstance(value, dict):
                merge_dicts(ref[key], value)
        else:
            ref[key] = value
    return ref
