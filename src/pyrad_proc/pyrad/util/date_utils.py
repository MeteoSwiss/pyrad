from datetime import datetime, timezone
import numpy as np
import cftime
import pandas as pd


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


# Strict formats you want to support explicitly
STRICT_DT_FORMATS = [
    "%Y%m%d%H%M%S",  # 20260413000414
    "%Y%m%d%H%M",  # 202604130004
    "%Y%m%d",  # 20260413
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%d %H:%M",
    "%Y-%m-%d",
    "%Y/%m/%d %H:%M:%S",
    "%Y/%m/%d %H:%M",
    "%Y/%m/%d",
]


def _try_strict_formats(s_nonnull: pd.Series, min_valid_ratio: float = 0.95):
    """
    Try exact datetime formats and return parsed series if one matches well.
    """
    s_str = s_nonnull.astype("string").str.strip()

    for fmt in STRICT_DT_FORMATS:
        parsed = pd.to_datetime(s_str, format=fmt, errors="coerce")
        valid_ratio = parsed.notna().mean()
        if valid_ratio >= min_valid_ratio:
            return parsed, fmt

    return None, None


def _looks_like_compact_datetime(s_nonnull: pd.Series) -> bool:
    """
    Detect compact numeric-like datetimes such as YYYYMMDDHHMMSS or YYYYMMDD.
    """
    s_str = s_nonnull.astype("string").str.strip()
    lengths = s_str.str.len()
    if lengths.nunique() != 1:
        return False

    length = lengths.iloc[0]
    if length not in (8, 12, 14):
        return False

    # all digits
    if not s_str.str.fullmatch(r"\d+").all():
        return False

    return True


def detect_datetime_columns(df: pd.DataFrame, min_valid_ratio: float = 0.95):
    datetime_cols = []

    for col in df.columns:
        s = df[col]

        # 1) already datetime
        if pd.api.types.is_datetime64_any_dtype(s):
            datetime_cols.append(col)
            continue

        s_nonnull = s.dropna()
        if s_nonnull.empty:
            continue

        parsed_full = None

        # 2) First try strict known formats on ALL columns, including numeric ones
        #    This is what will catch 20260413000414 correctly.
        if (
            pd.api.types.is_object_dtype(s)
            or pd.api.types.is_string_dtype(s)
            or pd.api.types.is_integer_dtype(s)
            or pd.api.types.is_float_dtype(s)
        ):
            # For float columns, only allow if values are integer-like
            if pd.api.types.is_float_dtype(s):
                float_ok = np.isclose(s_nonnull, np.floor(s_nonnull)).all()
                if not float_ok:
                    continue

            # avoid scientific notation / trailing .0 issues for numeric columns
            if pd.api.types.is_numeric_dtype(s):
                s_for_parse = s_nonnull.astype("Int64").astype("string")
            else:
                s_for_parse = s_nonnull.astype("string").str.strip()

            # Strong hint path for compact integer timestamps
            if _looks_like_compact_datetime(s_for_parse):
                parsed_nonnull, used_fmt = _try_strict_formats(
                    s_for_parse, min_valid_ratio=min_valid_ratio
                )
                if parsed_nonnull is not None:
                    parsed_full = pd.to_datetime(
                        (
                            s.astype("Int64").astype("string")
                            if pd.api.types.is_numeric_dtype(s)
                            else s.astype("string").str.strip()
                        ),
                        format=used_fmt,
                        errors="coerce",
                    )

            # General strict-format path
            if parsed_full is None:
                parsed_nonnull, used_fmt = _try_strict_formats(
                    s_for_parse, min_valid_ratio=min_valid_ratio
                )
                if parsed_nonnull is not None:
                    parsed_full = pd.to_datetime(
                        (
                            s.astype("Int64").astype("string")
                            if pd.api.types.is_numeric_dtype(s)
                            else s.astype("string").str.strip()
                        ),
                        format=used_fmt,
                        errors="coerce",
                    )

        # 3) Conservative fallback: only for string/object columns
        #    Do NOT do this for numeric columns, otherwise IDs may get parsed.
        if parsed_full is None and (
            pd.api.types.is_object_dtype(s) or pd.api.types.is_string_dtype(s)
        ):
            s_str = s_nonnull.astype("string").str.strip()

            # reject mostly numeric free-form strings unless they matched a strict format above
            numeric_ratio = pd.to_numeric(s_str, errors="coerce").notna().mean()
            if numeric_ratio < 0.8:
                parsed_nonnull = pd.to_datetime(s_str, errors="coerce")
                valid_ratio = parsed_nonnull.notna().mean()

                if valid_ratio >= min_valid_ratio:
                    year_min = parsed_nonnull.dt.year.min()
                    year_max = parsed_nonnull.dt.year.max()

                    if 1900 <= year_min <= 2100 and 1900 <= year_max <= 2100:
                        parsed_full = pd.to_datetime(
                            s.astype("string").str.strip(), errors="coerce"
                        )

        # 4) Final acceptance checks
        if parsed_full is not None:
            parsed_nonnull = parsed_full.dropna()
            if not parsed_nonnull.empty:
                year_min = parsed_nonnull.dt.year.min()
                year_max = parsed_nonnull.dt.year.max()

                if 1900 <= year_min <= 2100 and 1900 <= year_max <= 2100:
                    df[col] = parsed_full
                    datetime_cols.append(col)

    return datetime_cols
