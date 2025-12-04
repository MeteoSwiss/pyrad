from datetime import datetime, timezone
import numpy as np
import cftime


def cftodatetime(obj):
    """
    Convert output from netCDF4.num2date into timezone-aware UTC datetime(s).

    Parameters
    ----------
    obj : datetime, cftime object, array-like, or masked array
        Result of netCDF4.num2date.

    Returns
    -------
    datetime or numpy.ndarray of datetime
        Timezone-aware UTC datetime(s).
    """

    def convert_single(x):
        # Handle missing values
        if x is None:
            return None

        # If already a Python datetime
        if isinstance(x, datetime):
            # If naive -> add UTC timezone
            return (
                x.replace(tzinfo=timezone.utc)
                if x.tzinfo is None
                else x.astimezone(timezone.utc)
            )

        # Handle cftime objects
        if isinstance(x, cftime.datetime):
            return datetime(
                x.year,
                x.month,
                x.day,
                x.hour,
                x.minute,
                x.second,
                x.microsecond,
                tzinfo=timezone.utc,
            )

        raise TypeError(f"Unsupported type: {type(x)}")

    # If scalar
    if not isinstance(obj, (list, tuple, np.ndarray)):
        return convert_single(obj)

    # If array-like
    obj_arr = np.asarray(obj, dtype=object)

    return np.vectorize(convert_single, otypes=[object])(obj_arr)
