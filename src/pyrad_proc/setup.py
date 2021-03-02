#!/usr/bin/env python
"""Pyrad: Python Radar Toolkit

Pyrad is a Python module containing
the utilities that run the MeteoSwiss radar processing framework.
It is designed so that it accepts a growing number of radar data types.
The core of the processing is performed by the module Py-ART.

"""

import os
import shutil
import sys
import re
import subprocess
import glob
import builtins
from datetime import datetime
import getpass
import setuptools  # for 'develop' mode

DOCLINES = __doc__.split("\n")

CLASSIFIERS = """\
Development Status :: 0 - Prototype
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: BSD License
Programming Language :: Python
Programming Language :: Python :: 3
Programming Language :: Python :: 3.6
Programming Language :: Python :: 3.7
Programming Language :: C
Programming Language :: Cython
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Atmospheric Science
Operating System :: POSIX :: Linux
"""


NAME = 'pyrad_mch'
MAINTAINER = "MeteoSwiss Pyrad Developers"
MAINTAINER_EMAIL = "daniel.wolfensberger@meteoswiss.ch"
DESCRIPTION = DOCLINES[0]
LONG_DESCRIPTION = "\n".join(DOCLINES[2:])
URL = "https://github.com/MeteoSwiss/pyrad.git"
DOWNLOAD_URL = "https://github.com/MeteoSwiss/pyrad.git"
LICENSE = 'BSD'
CLASSIFIERS = filter(None, CLASSIFIERS.split('\n'))
PLATFORMS = ["Linux"]
MAJOR = 0
MINOR = 6
MICRO = 3
ISRELEASED = True
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)
SCRIPTS = glob.glob('scripts/*')
COMPILE_DATE_TIME = datetime.utcnow().strftime("%Y-%m-%d %H:%M")
USERNAME = getpass.getuser()


# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError as ee:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

# This is a bit hackish: we are setting a global variable so that the main
# pyrad __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet. While ugly, it's
# a lot more robust than what was previously being used.
builtins.__PYRAD_SETUP__ = True


def write_version_py(filename='pyrad/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM PYRAD_PROC SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
compile_date_time = '%(compile_date_time)s'
username = '%(username)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    # Adding the git rev number needs to be done inside write_version_py(),
    # otherwise the import of pyrad.version messes up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('../../.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('pyrad/version.py'):
        # must be a source distribution, use existing version file
        try:
            from pyrad.version import git_revision as GIT_REVISION
        except ImportError:
            raise ImportError("Unable to import git_revision. Try removing "
                              "pyrad/version.py and the build directory "
                              "before building.")
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev+' + GIT_REVISION[:7]

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'compile_date_time': COMPILE_DATE_TIME,
                       'username': USERNAME,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('pyrad')
    config.add_data_files(('pyrad', '*.txt'))

    return config


def setup_package():

    # rewrite version file
    write_version_py()

    from numpy.distutils.core import setup

    setup(
        name=NAME,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        license=LICENSE,
        classifiers=CLASSIFIERS,
        platforms=PLATFORMS,
        configuration=configuration,
        scripts=SCRIPTS,
    )

if __name__ == '__main__':
    setup_package()
