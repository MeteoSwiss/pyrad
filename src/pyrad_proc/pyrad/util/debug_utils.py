# Make sure that deprecation warnings get printed by default
import warnings as _warnings
import pdb
import traceback
import sys
import os
import inspect

_warnings.simplefilter("always", DeprecationWarning)
ORANGE = "\033[38;5;208m"
RESET = "\033[0m"


def orange_warning_format(message, category, filename, lineno, file=None, line=None):
    print(f"{ORANGE}{category.__name__}:{message}\n{RESET}")


_warnings.showwarning = orange_warning_format


def warn(message, use_debug=True):
    _warnings.warn(message)
    sys.stderr.flush()
    if os.environ.get("PYRAD_DEBUG", False) and use_debug:
        trace = traceback.format_stack()
        idx_warn = (["warn" in line for line in trace]).index(True)
        print(f"\033[91mException triggered from {trace[idx_warn]}\033[0m")
        parent_frame = inspect.currentframe().f_back
        pdb.Pdb().set_trace(parent_frame)


def custom_excepthook(exc_type, exc_value, exc_traceback):
    # Print the exception
    traceback.print_exception(exc_type, exc_value, exc_traceback)

    # If the environment variable is set, enter debugger
    if os.environ.get("PYRAD_DEBUG", False):
        print("\nStarting pdb due to unhandled exception...\n")
        pdb.post_mortem(exc_traceback)


def enable_debug_on_error():
    sys.excepthook = custom_excepthook
