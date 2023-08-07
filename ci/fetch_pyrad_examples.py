# -*- coding: utf-8 -*-
"""
Scriptcript used to install the pyrad examples in a test environment
configuration file that points to that data.
The test data is downloaded in the `PYRAD_DATA_PATH` environmental variable
"""

import gzip
import json
import os
import shutil
import sys
import time
from datetime import datetime, timedelta
from distutils.dir_util import copy_tree
from logging.handlers import RotatingFileHandler
from tempfile import NamedTemporaryFile, TemporaryDirectory
from urllib import request
from urllib.error import HTTPError
from zipfile import ZipFile


def download_pyrad_examples(dir_path, force=True):
    """
    Download pyrad examples from github.
    Parameters
    ----------
    dir_path: str
        Path to directory where the pyrad data will be placed.
    force: bool
        If the destination directory exits and force=False, a DirectoryNotEmpty
        exception if raised.
        If force=True, the data will we downloaded in the destination directory and may
        override existing files.
    """

    # Check if directory exists but is not empty
    if os.path.exists(dir_path) and os.path.isdir(dir_path):
        if os.listdir(dir_path) and not force:
            raise DirectoryNotEmpty(
                dir_path + "is not empty.\n"
                "Set force=True force the extraction of the files."
            )
    else:
        os.makedirs(dir_path)

    # NOTE:
    # The http response from github can either contain Content-Length (size of the file)
    # or use chunked Transfer-Encoding.
    # If Transfer-Encoding is chunked, then the Content-Length is not available since
    # the content is dynamically generated and we can't know the length a
    # priori easily.

    print("Downloading pyrad_examples from github.")
    tmp_file_name, _ = request.urlretrieve(
        "https://github.com/MeteoSwiss/pyrad-examples/archive/gh_actions.zip"
    )

    with ZipFile(tmp_file_name, "r") as zip_obj:
        tmp_dir = TemporaryDirectory()

        # Extract all the contents of zip file in the temp directory
        common_path = os.path.commonprefix(zip_obj.namelist())

        zip_obj.extractall(tmp_dir.name)

        copy_tree(os.path.join(tmp_dir.name, common_path), dir_path)


tox_test_data_dir = os.environ["PYRAD_EXAMPLES_PATH"]

download_pyrad_examples(tox_test_data_dir, force=True)
