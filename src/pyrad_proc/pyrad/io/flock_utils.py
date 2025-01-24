try:
    import fcntl

    FCNTL_AVAIL = True
except ModuleNotFoundError:
    import msvcrt  # For windows

    FCNTL_AVAIL = False

import os
import time
from warnings import warn
import errno


def lock_file(file):
    if FCNTL_AVAIL:
        while True:
            try:
                fcntl.flock(file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except OSError as e:
                if e.errno == errno.EAGAIN:
                    time.sleep(0.1)
                elif e.errno == errno.EBADF:
                    warn(
                        "WARNING: No file locking is possible (NFS mount?), "
                        + "expect strange issues with multiprocessing..."
                    )
                    break
                else:
                    raise
    else:
        while True:
            try:
                # Attempt to acquire an exclusive lock on the file
                msvcrt.locking(
                    os.open(file, os.O_RDWR | os.O_CREAT), msvcrt.LK_NBLCK, 0
                )
                break
            except OSError as e:
                if e.errno == 13:  # Permission denied
                    time.sleep(0.1)
                elif e.errno == 13:  # No such file or directory
                    warn(
                        "WARNING: No file locking is possible (NFS mount?), "
                        + "expect strange issues with multiprocessing..."
                    )
                    break
                else:
                    raise


def unlock_file(file):
    if os.name == "posix":
        fcntl.flock(file, fcntl.LOCK_UN)
    else:
        msvcrt.locking(file.fileno(), msvcrt.LK_UNLCK, 1)
