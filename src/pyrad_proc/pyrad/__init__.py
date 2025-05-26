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


    ORANGE = "\033[38;5;208m"
    RESET = "\033[0m"

    def orange_warning_format(message, category, filename, lineno, file=None, line=None):
        print(f"{ORANGE}{category.__name__}:{message}\n{RESET}")


    _warnings.showwarning = orange_warning_format

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

    # root level functions
    # non at the moment
