"""
Pyrad: The Python Radar Toolkit
=====================================

"""
import os

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
    # print out helpful message if build fails or importing from source tree
    # fvj built not checked for the moment
    # from . import __check_build

    # import subpackages
    from . import graph  # noqa
    from . import io  # noqa
    from . import proc  # noqa
    from . import prod  # noqa
    from . import util  # noqa
    from . import flow  # noqa

    if os.environ:
        from .util import enable_debug_on_error

        enable_debug_on_error()

    # root level functions
    # non at the moment
