try:
    import fcntl
    FCNTL_AVAIL = True
except ModuleNotFoundError:
    import msvcrt  # For windows
    FCNTL_AVAIL = False

import os
import time
from warnings import warn

def lock_file(file):
    if os.name == 'posix':
        while True:
            try:
                fcntl.flock(file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except OSError as e:
                    if e.errno == 11:
                        time.sleep(0.1)
                    else:
                        warn("File locking failed, proceeding without locking.")
                        break
    else:
        while True:
            try:
                msvcrt.locking(file.fileno(), msvcrt.LK_NBLCK, 1)
                break
            except OSError:
                time.sleep(0.1)

def unlock_file(file):
    if os.name == 'posix':
        fcntl.flock(file, fcntl.LOCK_UN)
    else:
        msvcrt.locking(file.fileno(), msvcrt.LK_UNLCK, 1)