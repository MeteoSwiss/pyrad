#!/usr/bin/env bash
set -euo pipefail
python3 parse_pyrad_name_mappings.py
python3 parse_pyrad_processes.py
python3 parse_pyrad_products.py
