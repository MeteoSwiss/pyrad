"""
Pyrad: The Python Radar Toolkit
=====================================

"""

# Detect if we're being called as part of Pyrad's setup procedure
try:
    __PYRAD_SETUP__
except NameError:
    __PYRAD_SETUP__ = False

if __PYRAD_SETUP__:
    import sys as _sys
    _sys.stderr.write("Running from Pyrad source directory.\n")
    del _sys
else:

    # Make sure that deprecation warnings get printed by default
    import warnings as _warnings
    _warnings.simplefilter("always", DeprecationWarning)

    # print out helpful message if build fails or importing from source tree
    # fvj built not checked for the moment
    # from . import __check_build

    # versioning
    from .version import git_revision as __git_revision__
    from .version import version as __version__

    # import subpackages
    from . import graph
    from . import io
    from . import proc
    from . import prod
    from . import util
    from . import flow
    # root level functions
    # non at the moment
