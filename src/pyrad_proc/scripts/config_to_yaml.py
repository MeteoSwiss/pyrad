import argparse
import sys
import yaml
from typing import Tuple, Any

import pyrad
from pyrad.io.config import read_config


def convert_to_yaml(input_file):
    config = read_config(input_file)
    return yaml.dump(config, sort_keys=False, allow_unicode=True)


def main():
    parser = argparse.ArgumentParser(
        description="Convert legacy config format to YAML."
    )
    parser.add_argument("input_file", help="Path to the legacy config file")
    parser.add_argument(
        "-o", "--output", help="Path to write the YAML output (default: stdout)"
    )

    args = parser.parse_args()
    try:
        yaml_output = convert_to_yaml(args.input_file)
        if args.output:
            with open(args.output, "w") as out_file:
                out_file.write(yaml_output)
        else:
            print(yaml_output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
