# -*- coding: utf-8 -*-
"""
Scriptcript used to install the pyrad examples in a test environment
configuration file that points to that data.
The test data is downloaded in the `PYRAD_DATA_PATH` environmental variable
"""

import os

from pyrad.test import download_pyrad_examples

tox_test_data_dir = os.environ["PYSTEPS_DATA_PATH"]

download_pyrad_examples(tox_test_data_dir, force=True)
)

