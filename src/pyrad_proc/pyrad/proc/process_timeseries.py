"""
pyrad.proc.process_timeseries
=============================

Functions to obtain time series of radar data

.. autosummary::
    :toctree: generated/

    process_point_measurement
    process_multiple_points
    process_qvp
    process_rqvp
    process_evp
    process_svp
    process_time_height
    process_ts_along_coord

"""

import datetime
from ..util import warn
import numpy as np
from netCDF4 import num2date

import pyart

from ..util.date_utils import cftodatetime
from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.read_data_sensor import read_coord_sensors


def process_point_measurement(procstatus, dscfg, radar_list=None):
    """
    Obtains the radar data at a point location. Multiple radars are supported only
    if latlon is True.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement,
            can be any datatype supported by pyrad and available in the data.
            If multiple radars are provided, the datatype must be the same for all
            of them.
        agg_method : string. Dataset keyword
            Which aggregation method to use to combine data that is within the
            tolerance in AziTol, EleTol and RngTol
            'nearest' : will only get nearest point to prescribed lon/lat/alt or
                ele/azi/rng
            'nearest_valid' : will only get the nearest valid point, within
                tolerance AziTol, EleTol and RngTol to prescribed
                lon/lat/alt or ele/azi/rng (ignore missing values)
            'average' : will average (while ignore missing values), all values
                that fall within the tolerance in AziTol, EleTol and RngTol
            'none' : will not perform any averaging and will get all values that
                fall within the tolerance in AziTol, EleTol and RngTol, each
                with its individual timestamp
            'min' : minimum value within the tolerance in AziTol, EleTol and RngTol
            'max' : maximum value within the tolerance in AziTol, EleTol and RngTol
            'maxZH' : value at the location of the maximum in reflectivity within the tolerance in AziTol, EleTol and RngTol.
                for each radar. Accepted reflectivity datatype names are: dBuZ, dBZ, dBZv, dBuZv, dBZc, dBZvc
            Default is 'nearest'
        latlon : boolean. Dataset keyword
            if True position is obtained from latitude, longitude information,
            otherwise position is obtained from antenna coordinates
            (range, azimuth, elevation).
        truealt : boolean. Dataset keyword
            if True the user input altitude is used to determine the point of
            interest.
            if False use the altitude at a given radar elevation ele over the
            point of interest. Default is False.
        lon : float. Dataset keyword
            the longitude [deg]. Use when latlon is True.
        lat : float. Dataset keyword
            the latitude [deg]. Use when latlon is True.
        alt : float. Dataset keyword
            altitude [m MSL]. Use when latlon is True and truealt is True
        ele : float. Dataset keyword
            radar elevation [deg]. Use when latlon is False or when latlon is
            True and truealt is False
        azi : float. Dataset keyword
            radar azimuth [deg]. Use when latlon is False
        rng : float. Dataset keyword
            range from radar [m]. Use when latlon is False
        AziTol : float. Dataset keyword
            azimuthal tolerance to determine which radar azimuth to use [deg]
        EleTol : float. Dataset keyword
            elevation tolerance to determine which radar elevation to use [deg]
        RngTol : float. Dataset keyword
            range tolerance to determine which radar bin to use [m]
        fill_value : float or None
            If not None masked values are going to be filled by this value

    radar_list : list of Radar objects
          Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the data and metadata at the point of interest
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    refl_datatypes = {"dBuZ", "dBZ", "dBZv", "dBuZv", "dBZc", "dBZvc"}
    valid_agg_methods = [
        "nearest",
        "nearest_valid",
        "average",
        "max",
        "min",
        "maxZH",
        "none",
    ]

    # ------------------------------------------------------------------
    # Parse requested radars + datatype(s)
    # ------------------------------------------------------------------
    req_by_radar = {}
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        req_by_radar.setdefault(radarnr, []).append(datatype)

    req = list(req_by_radar.items())

    if not dscfg.get("latlon", False) and len(req) > 1:
        warn("Multi-radar point measurement requires dscfg['latlon']=True.")
        return None, None

    agg_method = dscfg.get("agg_method", "nearest")
    if agg_method not in valid_agg_methods:
        warn("Invalid agg_method. Using 'nearest'.")
        agg_method = "nearest"

    # ------------------------------------------------------------------
    # Finalization
    # ------------------------------------------------------------------
    if procstatus == 2:
        if dscfg["initialized"] == 0:
            return None, None

        out = {}
        ind_out = {}

        for radarnr, datatype_list in req:
            if radarnr not in dscfg["global_data"]:
                out[radarnr] = None
                ind_out[radarnr] = None
                continue

            gd = dscfg["global_data"][radarnr]

            out[radarnr] = {
                "time": gd["time"],
                "value": gd["value"],
                "ref_time": gd["ref_time"],
                "datatype": gd["datatype"],
                "point_coordinates_WGS84_lon_lat_alt": gd[
                    "point_coordinates_WGS84_lon_lat_alt"
                ],
                "antenna_coordinates_az_el_r": gd["antenna_coordinates_az_el_r"],
                "radarnr": radarnr,
            }

            ind_out[radarnr] = gd.get("ind_rad", None)

        out["final"] = True
        return out, ind_out

    # ------------------------------------------------------------------
    # Processing
    # ------------------------------------------------------------------
    if radar_list is None:
        warn("ERROR: radar_list is None")
        return None, None

    lon = dscfg["lon"]
    lat = dscfg["lat"]
    fill_value = dscfg.get("fill_value", None)

    out = {}
    ind_out = {}

    for radarnr, datatype_list in req:
        # --------------------------------------------------------------
        # Determine main field and, if needed, reflectivity field
        # --------------------------------------------------------------
        if agg_method == "maxZH":
            refl_candidates = [dt for dt in datatype_list if dt in refl_datatypes]
            value_candidates = [dt for dt in datatype_list if dt not in refl_datatypes]

            if not refl_candidates:
                warn(
                    f"{radarnr}: agg_method='maxZH' requires one reflectivity "
                    f"datatype in {sorted(refl_datatypes)}."
                )
                out[radarnr] = None
                continue

            refl_datatype = refl_candidates[0]
            refl_field_name = get_fieldname_pyart(refl_datatype)

            # If only reflectivity is provided, return reflectivity at max ZH.
            if value_candidates:
                datatype = value_candidates[0]
            else:
                datatype = refl_datatype
        else:
            datatype = datatype_list[0]
            refl_datatype = None
            refl_field_name = None

        field_name = get_fieldname_pyart(datatype)

        # --------------------------------------------------------------
        # Map radar number to radar_list index
        # --------------------------------------------------------------
        ind_rad = int(radarnr[5:8]) - 1
        ind_out[radarnr] = ind_rad

        if ind_rad < 0 or ind_rad >= len(radar_list) or radar_list[ind_rad] is None:
            warn(f"ERROR: No valid radar for {radarnr}")
            out[radarnr] = None
            continue

        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn(f"{radarnr}: Field '{field_name}' not available. Skipping.")
            out[radarnr] = None
            continue

        if agg_method == "maxZH" and refl_field_name not in radar.fields:
            warn(
                f"{radarnr}: Reflectivity field '{refl_field_name}' not available. "
                "Skipping."
            )
            out[radarnr] = None
            continue

        # --------------------------------------------------------------
        # Projection: lon/lat -> x/y
        # --------------------------------------------------------------
        projparams = {
            "proj": "pyart_aeqd",
            "lon_0": radar.longitude["data"],
            "lat_0": radar.latitude["data"],
        }

        x, y = pyart.core.geographic_to_cartesian(lon, lat, projparams)

        # --------------------------------------------------------------
        # Determine altitude
        # --------------------------------------------------------------
        if not dscfg.get("truealt", False):
            if hasattr(pyart.config, "_KE"):
                ke = pyart.config._KE
            else:
                ke = 4.0 / 3.0

            a = 6378100.0
            re = a * ke

            elrad = dscfg["ele"] * np.pi / 180.0
            r_ground = np.sqrt(x**2.0 + y**2.0)
            r_tmp = r_ground / np.cos(elrad)

            alt = (
                radar.altitude["data"]
                + np.sqrt(r_tmp**2.0 + re**2.0 + 2.0 * r_tmp * re * np.sin(elrad))
                - re
            )
            alt = alt[0]
        else:
            alt = dscfg["alt"]

        # --------------------------------------------------------------
        # Cartesian -> antenna coordinates
        # --------------------------------------------------------------
        r, az, el = pyart.core.cartesian_to_antenna(x, y, alt - radar.altitude["data"])

        r = r[0]
        az = az[0]
        el = el[0]

        # --------------------------------------------------------------
        # Tolerance checks
        # --------------------------------------------------------------
        d_az = np.abs(radar.azimuth["data"] - az)
        if np.min(d_az) > dscfg["AziTol"]:
            warn(
                f"{radarnr}: No bin found (az tol). "
                f"min d_az={np.min(d_az)} > {dscfg['AziTol']}"
            )
            out[radarnr] = None
            continue

        d_el = np.abs(radar.elevation["data"] - el)
        if np.min(d_el) > dscfg["EleTol"]:
            warn(
                f"{radarnr}: No bin found (el tol). "
                f"min d_el={np.min(d_el)} > {dscfg['EleTol']}"
            )
            out[radarnr] = None
            continue

        d_r = np.abs(radar.range["data"] - r)
        if np.min(d_r) > dscfg["RngTol"]:
            warn(
                f"{radarnr}: No bin found (range tol). "
                f"min d_r={np.min(d_r)} > {dscfg['RngTol']}"
            )
            out[radarnr] = None
            continue

        # --------------------------------------------------------------
        # Extract / aggregate value
        # --------------------------------------------------------------
        is_pointing = False
        selected_range_indices = None

        # Determine all rays and ranges falling inside the tolerances first.
        # This is also needed to detect pointing scans for agg_method='nearest'.
        tol_ray_arr = np.where((d_az <= dscfg["AziTol"]) & (d_el <= dscfg["EleTol"]))[0]

        tol_r_arr = np.where(d_r <= dscfg["RngTol"])[0]

        if tol_ray_arr.size == 0 or tol_r_arr.size == 0:
            warn(f"{radarnr}: No gates found within tolerance.")
            out[radarnr] = None
            continue

        # In a pointing scan, every radar ray falls within the angular tolerance.
        is_pointing = tol_ray_arr.size == np.asarray(radar.azimuth["data"]).size

        field_data = radar.fields[field_name]["data"]

        if agg_method == "nearest":
            # Nearest range gate is common to every ray.
            nearest_local_r = np.argmin(d_r[tol_r_arr])
            nearest_ind_r = tol_r_arr[nearest_local_r]

            if is_pointing:
                # Preserve every ray and its timestamp.
                ind_ray_arr = tol_ray_arr
                ind_r_arr = np.atleast_1d(nearest_ind_r)

                val = field_data[ind_ray_arr, nearest_ind_r]

                selected_range_indices = np.full(
                    ind_ray_arr.size,
                    nearest_ind_r,
                    dtype=int,
                )

            else:
                # Normal scan: retain the original nearest-ray behaviour,
                # but restrict the search to rays inside the tolerance.
                angular_distance = np.abs(
                    radar.azimuth["data"][tol_ray_arr] - az
                ) + np.abs(radar.elevation["data"][tol_ray_arr] - el)

                nearest_local_ray = np.argmin(angular_distance)
                nearest_ind_ray = tol_ray_arr[nearest_local_ray]

                ind_ray_arr = np.atleast_1d(nearest_ind_ray)
                ind_r_arr = np.atleast_1d(nearest_ind_r)

                val = field_data[nearest_ind_ray, nearest_ind_r]

                selected_range_indices = ind_r_arr.copy()

        else:
            ind_ray_arr = tol_ray_arr
            ind_r_arr = tol_r_arr

            data_region = field_data[ind_ray_arr][:, ind_r_arr]

            if is_pointing:
                # Preserve the ray/time dimension and aggregate only over range.
                nearest_local_r = np.argmin(d_r[ind_r_arr])
                nearest_ind_r = ind_r_arr[nearest_local_r]

                if agg_method == "average":
                    val = np.ma.mean(data_region, axis=1)

                    selected_range_indices = np.full(
                        ind_ray_arr.size,
                        nearest_ind_r,
                        dtype=int,
                    )

                elif agg_method == "max":
                    val = np.ma.max(data_region, axis=1)

                    # Store the actual range of the maximum for each ray.
                    local_r_indices = np.ma.argmax(data_region, axis=1)
                    selected_range_indices = ind_r_arr[local_r_indices]

                elif agg_method == "min":
                    val = np.ma.min(data_region, axis=1)

                    # Store the actual range of the minimum for each ray.
                    local_r_indices = np.ma.argmin(data_region, axis=1)
                    selected_range_indices = ind_r_arr[local_r_indices]

                elif agg_method == "none":
                    # Shape:
                    #     nrays x nranges_within_tolerance
                    val = data_region

                    # There is no unique range for each ray. Use the nearest
                    # requested range in the coordinate metadata.
                    selected_range_indices = np.full(
                        ind_ray_arr.size,
                        nearest_ind_r,
                        dtype=int,
                    )

                elif agg_method == "maxZH":
                    refl_region = radar.fields[refl_field_name]["data"][ind_ray_arr][
                        :, ind_r_arr
                    ]

                    valid_count = np.ma.count(
                        refl_region,
                        axis=1,
                    )

                    val = np.ma.masked_all(
                        ind_ray_arr.size,
                        dtype=data_region.dtype,
                    )

                    selected_range_indices = np.full(
                        ind_ray_arr.size,
                        nearest_ind_r,
                        dtype=int,
                    )

                    valid_rays = valid_count > 0

                    if not np.any(valid_rays):
                        warn(
                            f"{radarnr}: No valid reflectivity values "
                            "found for maxZH."
                        )
                        out[radarnr] = None
                        continue

                    valid_ray_positions = np.where(valid_rays)[0]

                    max_local_r = np.ma.argmax(
                        refl_region[valid_rays],
                        axis=1,
                    )

                    val[valid_rays] = data_region[
                        valid_ray_positions,
                        max_local_r,
                    ]

                    selected_range_indices[valid_rays] = ind_r_arr[max_local_r]

                elif agg_method == "nearest_valid":
                    range_distance = np.abs(radar.range["data"][ind_r_arr] - r)

                    val = np.ma.masked_all(
                        ind_ray_arr.size,
                        dtype=data_region.dtype,
                    )

                    selected_range_indices = np.full(
                        ind_ray_arr.size,
                        nearest_ind_r,
                        dtype=int,
                    )

                    data_mask = np.ma.getmaskarray(data_region)

                    for iray in range(ind_ray_arr.size):
                        valid_ranges = ~data_mask[iray]

                        if not np.any(valid_ranges):
                            continue

                        candidate_local_indices = np.where(valid_ranges)[0]

                        best_local_r = candidate_local_indices[
                            np.argmin(range_distance[candidate_local_indices])
                        ]

                        val[iray] = data_region[
                            iray,
                            best_local_r,
                        ]

                        selected_range_indices[iray] = ind_r_arr[best_local_r]

                    if np.ma.count(val) == 0:
                        warn(f"{radarnr}: No valid gates found " "within tolerance.")
                        out[radarnr] = None
                        continue

            else:
                # ----------------------------------------------------------
                # Normal scanning case: aggregate over rays and ranges
                # ----------------------------------------------------------
                if agg_method == "average":
                    val = np.ma.mean(data_region)

                elif agg_method == "max":
                    if np.ma.count(data_region) == 0:
                        warn(f"{radarnr}: No valid gates found " "within tolerance.")
                        out[radarnr] = None
                        continue

                    local_index = np.unravel_index(
                        np.ma.argmax(data_region),
                        data_region.shape,
                    )

                    val = data_region[local_index]

                    ind_ray_arr = np.atleast_1d(ind_ray_arr[local_index[0]])
                    ind_r_arr = np.atleast_1d(ind_r_arr[local_index[1]])
                    selected_range_indices = ind_r_arr.copy()

                elif agg_method == "min":
                    if np.ma.count(data_region) == 0:
                        warn(f"{radarnr}: No valid gates found " "within tolerance.")
                        out[radarnr] = None
                        continue

                    local_index = np.unravel_index(
                        np.ma.argmin(data_region),
                        data_region.shape,
                    )

                    val = data_region[local_index]

                    ind_ray_arr = np.atleast_1d(ind_ray_arr[local_index[0]])
                    ind_r_arr = np.atleast_1d(ind_r_arr[local_index[1]])
                    selected_range_indices = ind_r_arr.copy()

                elif agg_method == "none":
                    val = data_region

                elif agg_method == "maxZH":
                    refl_region = radar.fields[refl_field_name]["data"][ind_ray_arr][
                        :, ind_r_arr
                    ]

                    if np.ma.count(refl_region) == 0:
                        warn(
                            f"{radarnr}: No valid reflectivity "
                            "values found for maxZH."
                        )
                        out[radarnr] = None
                        continue

                    zh_flat_idx = np.ma.argmax(refl_region)
                    zh_idx = np.unravel_index(
                        zh_flat_idx,
                        refl_region.shape,
                    )

                    val = data_region[zh_idx]

                    ind_ray_arr = np.atleast_1d(ind_ray_arr[zh_idx[0]])
                    ind_r_arr = np.atleast_1d(ind_r_arr[zh_idx[1]])

                    selected_range_indices = ind_r_arr.copy()

                elif agg_method == "nearest_valid":
                    # Find the closest geometrical gate among the candidates.
                    ind_clo_ray = np.argmin(
                        np.abs(radar.azimuth["data"] - az)
                        + np.abs(radar.elevation["data"] - el)
                    )
                    ind_clo_r = np.argmin(np.abs(radar.range["data"] - r))

                    clo_x = radar.gate_x["data"][
                        ind_clo_ray,
                        ind_clo_r,
                    ]
                    clo_y = radar.gate_y["data"][
                        ind_clo_ray,
                        ind_clo_r,
                    ]
                    clo_z = radar.gate_z["data"][
                        ind_clo_ray,
                        ind_clo_r,
                    ]

                    dist_x = np.abs(
                        radar.gate_x["data"][ind_ray_arr][:, ind_r_arr] - clo_x
                    )
                    dist_y = np.abs(
                        radar.gate_y["data"][ind_ray_arr][:, ind_r_arr] - clo_y
                    )
                    dist_z = np.abs(
                        radar.gate_z["data"][ind_ray_arr][:, ind_r_arr] - clo_z
                    )

                    dist_to_clo = np.sqrt(dist_x**2 + dist_y**2 + dist_z**2)

                    vmask = field_data[ind_ray_arr][:, ind_r_arr]

                    dist_to_clo[np.ma.getmaskarray(vmask)] = np.inf

                    if not np.isfinite(dist_to_clo).any():
                        warn(f"{radarnr}: No valid gates found " "within tolerance.")
                        out[radarnr] = None
                        continue

                    clo_valid_idx = np.unravel_index(
                        np.argmin(dist_to_clo),
                        dist_to_clo.shape,
                    )

                    val = vmask[clo_valid_idx]

                    ind_ray_arr = np.atleast_1d(ind_ray_arr[clo_valid_idx[0]])
                    ind_r_arr = np.atleast_1d(ind_r_arr[clo_valid_idx[1]])

                    selected_range_indices = ind_r_arr.copy()

                if selected_range_indices is None:
                    nearest_local_r = np.argmin(d_r[ind_r_arr])

                    selected_range_indices = np.full(
                        ind_ray_arr.size,
                        ind_r_arr[nearest_local_r],
                        dtype=int,
                    )

        # --------------------------------------------------------------
        # Normalize output shape
        # --------------------------------------------------------------
        if is_pointing:
            # For aggregated methods, val has shape (nrays,).
            # For 'none', val has shape (nrays, nranges_in_tolerance).
            val = np.ma.asarray(val)
        else:
            val = np.ma.asarray(val).ravel()

        # --------------------------------------------------------------
        # Antenna coordinates actually used
        # --------------------------------------------------------------
        ant_coord = np.empty(
            (3, ind_ray_arr.size),
            dtype=np.float32,
        )

        ant_coord[0, :] = radar.azimuth["data"][ind_ray_arr]
        ant_coord[1, :] = radar.elevation["data"][ind_ray_arr]
        ant_coord[2, :] = radar.range["data"][selected_range_indices]

        # --------------------------------------------------------------
        # Fill masked values if requested
        # --------------------------------------------------------------
        if fill_value is not None:
            val = np.ma.asarray(val).filled(fill_value)

        # --------------------------------------------------------------
        # Time for selected ray(s)
        # --------------------------------------------------------------
        time = cftodatetime(
            num2date(
                radar.time["data"][ind_ray_arr],
                radar.time["units"],
                radar.time["calendar"],
            )
        )

        time = np.atleast_1d(time)

        # --------------------------------------------------------------
        # Persistent storage
        # --------------------------------------------------------------
        if "global_data" not in dscfg:
            dscfg["global_data"] = {}

        if radarnr not in dscfg["global_data"]:
            dscfg["global_data"][radarnr] = {
                "point_coordinates_WGS84_lon_lat_alt": [
                    lon,
                    lat,
                    alt,
                ],
                "antenna_coordinates_az_el_r": [
                    az,
                    el,
                    r,
                ],
                "time": [],
                "value": [],
                "ref_time": dscfg["timeinfo"],
                "ind_rad": ind_rad,
                "datatype": datatype,
            }
            dscfg["initialized"] = 1

        gd = dscfg["global_data"][radarnr]

        gd["ref_time"] = dscfg["timeinfo"]
        gd["datatype"] = datatype

        if is_pointing:
            # Store one entry for every ray timestamp.
            #
            # For average/min/max/maxZH/nearest_valid:
            #     each value entry is a scalar.
            #
            # For agg_method='none':
            #     each value entry is a 1-D array containing all
            #     selected range gates for that ray.
            gd["time"].extend(time.tolist())

            if val.ndim == 1:
                gd["value"].extend([value for value in val])
            else:
                gd["value"].extend([val[iray].copy() for iray in range(val.shape[0])])
        else:
            gd["value"].append(val[0])
            gd["time"].append(time[0])

        # --------------------------------------------------------------
        # Per-radar output
        # --------------------------------------------------------------
        new_dataset = {
            "value": dscfg["global_data"][radarnr]["value"],
            "datatype": datatype,
            "time": dscfg["global_data"][radarnr]["time"],
            "ref_time": dscfg["timeinfo"],
            "point_coordinates_WGS84_lon_lat_alt": [
                lon,
                lat,
                alt,
            ],
            "antenna_coordinates_az_el_r": [
                az,
                el,
                r,
            ],
            "used_antenna_coordinates_az_el_r": ant_coord,
            "radarnr": radarnr,
        }

        out[radarnr] = new_dataset
    out["final"] = False
    return out, ind_out


def process_multiple_points(procstatus, dscfg, radar_list=None):
    """
    Obtains the radar data at multiple points. The points are defined in a file

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement,
            can be any datatype supported by pyrad and available in the data
        truealt : boolean. Dataset keyword
            if True the user input altitude is used to determine the point of
            interest.
            if False use the altitude at a given radar elevation ele over the
            point of interest. Default is False.
        coord_fname : string
            File name containing the points coordinates
        alt_points : float. Dataset keyword
            altitude [m MSL]. Use when latlon is True.
        ele_points : float. Dataset keyword
            radar elevation [deg]. Use when latlon is False or when latlon is
            True and truealt is False
        AziTol : float. Dataset keyword
            azimuthal tolerance to determine which radar azimuth to use [deg]
        EleTol : float. Dataset keyword
            elevation tolerance to determine which radar elevation to use [deg]
        RngTol : float. Dataset keyword
            range tolerance to determine which radar bin to use [m]

    radar_list : list of Radar objects
          Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the data and metadata at the point of interest
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        break
    field_name = get_fieldname_pyart(datatype)
    ind_rad = int(radarnr[5:8]) - 1

    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn("ERROR: No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    if field_name not in radar.fields:
        warn(
            "Unable to extract point measurement information. " + "Field not available"
        )
        return None, None

    # default parameters
    truealt = dscfg.get("truealt", False)
    ele_points = dscfg.get("ele_points", 1.0)
    alt_points = dscfg.get("alt_points", 0.0)
    AziTol = dscfg.get("AziTol", 0.25)
    EleTol = dscfg.get("EleTol", 1.0)
    RngTol = dscfg.get("RngTol", 120.0)

    projparams = dict()
    projparams.update({"proj": "pyart_aeqd"})
    projparams.update({"lon_0": radar.longitude["data"]})
    projparams.update({"lat_0": radar.latitude["data"]})

    lat, lon, point_id = read_coord_sensors(dscfg["coord_fname"])
    x, y = pyart.core.geographic_to_cartesian(lon, lat, projparams)
    npoints = lon.size

    if not truealt:
        ke = 4.0 / 3.0  # constant for effective radius
        a = 6378100.0  # earth radius
        re = a * ke  # effective radius

        elrad = ele_points * np.pi / 180.0
        r_ground = np.sqrt(x**2.0 + y**2.0)
        r = r_ground / np.cos(elrad)
        alt = (
            radar.altitude["data"]
            + np.sqrt(r**2.0 + re**2.0 + 2.0 * r * re * np.sin(elrad))
            - re
        )
    else:
        alt = alt_points + np.zeros(npoints)

    r, az, el = pyart.core.cartesian_to_antenna(x, y, alt - radar.altitude["data"])

    val = np.ma.masked_all(npoints)
    time = np.ma.masked_all(npoints, dtype=datetime.datetime)
    used_lon = np.ma.masked_all(npoints)
    used_lat = np.ma.masked_all(npoints)
    used_alt = np.ma.masked_all(npoints)
    for ind in range(npoints):
        d_az = np.abs(radar.azimuth["data"] - az[ind])
        if np.min(d_az) > AziTol:
            warn(
                " No radar bin found for point (az, el, r):("
                + str(az[ind])
                + ", "
                + str(el[ind])
                + ", "
                + str(r[ind])
                + "). Minimum distance to radar azimuth "
                + str(d_az)
                + f" larger than azimuth tolerance {AziTol}"
            )
            continue

        d_el = np.abs(radar.elevation["data"] - el[ind])
        if np.min(d_el) > EleTol:
            warn(
                " No radar bin found for point (az, el, r):("
                + str(az[ind])
                + ", "
                + str(el[ind])
                + ", "
                + str(r[ind])
                + "). Minimum distance to radar elevation "
                + str(d_el)
                + f" larger than elevation tolerance {EleTol}"
            )
            continue

        d_r = np.abs(radar.range["data"] - r[ind])
        if np.min(d_r) > RngTol:
            warn(
                " No radar bin found for point (az, el, r):("
                + str(az[ind])
                + ", "
                + str(el[ind])
                + ", "
                + str(r[ind])
                + "). Minimum distance to radar range bin "
                + str(d_r)
                + f" larger than range tolerance {RngTol}"
            )
            continue

        ind_ray = np.argmin(
            np.abs(radar.azimuth["data"] - az[ind])
            + np.abs(radar.elevation["data"] - el[ind])
        )
        ind_r = np.argmin(np.abs(radar.range["data"] - r[ind]))

        used_lon[ind] = radar.gate_longitude["data"][ind_ray, ind_r]
        used_lat[ind] = radar.gate_latitude["data"][ind_ray, ind_r]
        used_alt[ind] = radar.gate_altitude["data"][ind_ray, ind_r]

        val[ind] = radar.fields[field_name]["data"].data[ind_ray, ind_r]
        time[ind] = cftodatetime(
            num2date(
                radar.time["data"][ind_ray], radar.time["units"], radar.time["calendar"]
            )
        )

    # prepare for exit
    new_dataset = dict()
    new_dataset.update({"value": val})
    new_dataset.update({"datatype": datatype})
    new_dataset.update({"time": time})
    new_dataset.update({"ref_time": dscfg["timeinfo"]})
    new_dataset.update({"point_coordinates_WGS84_lon": lon})
    new_dataset.update({"point_coordinates_WGS84_lat": lat})
    new_dataset.update({"point_coordinates_WGS84_alt": alt})
    new_dataset.update({"antenna_coordinates_az": az})
    new_dataset.update({"antenna_coordinates_el": el})
    new_dataset.update({"antenna_coordinates_r": r})
    new_dataset.update({"used_point_coordinates_WGS84_lon": used_lon})
    new_dataset.update({"used_point_coordinates_WGS84_lat": used_lat})
    new_dataset.update({"used_point_coordinates_WGS84_alt": used_alt})
    new_dataset.update({"point_id": point_id})
    new_dataset.update({"final": False})

    return new_dataset, ind_rad


def process_qvp(procstatus, dscfg, radar_list=None):
    """
    Computes quasi vertical profiles, by averaging over height levels
    PPI data.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement,
            can be any datatype supported by pyrad and available in the data
        angle : int or float
            If the radar object contains a PPI volume, the sweep number to
            use, if it contains an RHI volume the elevation angle.
            Default 0.
        ang_tol : float
            If the radar object contains an RHI volume, the tolerance in the
            elevation angle for the conversion into PPI
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 50
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to accept average. Default 30.
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the QVP and a keyword stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    Reference
    ---------
    Ryzhkov A., Zhang P., Reeves H., Kumjian M., Tschallener T., Trömel S.,
    Simmer C. 2016: Quasi-Vertical Profiles: A New Way to Look at Polarimetric
    Radar Data. JTECH vol. 33 pp 551-562

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        field_names = []
        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8]) - 1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn("ERROR: No valid radar")
            return None, None

        radar = radar_list[ind_rad]

        # default parameters
        angle = dscfg.get("angle", 0)
        ang_tol = dscfg.get("ang_tol", 1.0)
        hmax = dscfg.get("hmax", 10000.0)
        hres = dscfg.get("hres", 50.0)
        avg_type = dscfg.get("avg_type", "mean")
        nvalid_min = dscfg.get("nvalid_min", 30)
        interp_kind = dscfg.get("interp_kind", "none")

        # initialize dataset
        if not dscfg["initialized"]:
            qvp = pyart.retrieve.compute_qvp(
                radar,
                field_names,
                ref_time=dscfg["timeinfo"],
                angle=angle,
                ang_tol=ang_tol,
                hmax=hmax,
                hres=hres,
                avg_type=avg_type,
                nvalid_min=nvalid_min,
                interp_kind=interp_kind,
                qvp=None,
            )

            if qvp is None:
                warn("Unable to compute QVP")
                return None, None

            global_dict = dict()
            global_dict.update({"start_time": dscfg["timeinfo"]})
            global_dict.update({"radar_out": qvp})
            dscfg["global_data"] = global_dict
            dscfg["initialized"] = 1
        else:
            qvp = pyart.retrieve.compute_qvp(
                radar,
                field_names,
                ref_time=dscfg["timeinfo"],
                angle=angle,
                ang_tol=ang_tol,
                hmax=hmax,
                hres=hres,
                avg_type=avg_type,
                nvalid_min=nvalid_min,
                interp_kind=interp_kind,
                qvp=dscfg["global_data"]["radar_out"],
            )

            if qvp is None:
                warn("Unable to compute QVP")
                return None, None

            dscfg["global_data"]["radar_out"] = qvp

        new_dataset = dict()
        new_dataset.update({"radar_out": qvp})
        new_dataset.update({"radar_type": "temporal"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg["initialized"]:
            return None, None

        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8]) - 1

        qvp = dscfg["global_data"]["radar_out"]

        new_dataset = dict()
        new_dataset.update({"radar_out": qvp})
        new_dataset.update({"radar_type": "final"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad


def process_rqvp(procstatus, dscfg, radar_list=None):
    """
    Computes range defined quasi vertical profiles, by averaging over height
    levels PPI data.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement,
            can be any datatype supported by pyrad and available in the data
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 2.
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to accept average. Default 30.
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'nearest'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation
        rmax : float
            ground range up to which the data is intended for use [m].
            Default 50000.
        weight_power : float
            Power p of the weighting function 1/abs(grng-(rmax-1))**p given to
            the data outside the desired range. -1 will set the weight to 0.
            Default 2.

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the QVP and a keyword stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    Reference
    ---------
    Tobin D.M., Kumjian M.R. 2017: Polarimetric Radar and Surface-Based
    Precipitation-Type Observations of ice Pellet to Freezing Rain
    Transitions. Weather and Forecasting vol. 32 pp 2065-2082

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        field_names = []
        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8]) - 1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn("ERROR: No valid radar")
            return None, None

        radar = radar_list[ind_rad]

        # default parameters
        hmax = dscfg.get("hmax", 10000.0)
        hres = dscfg.get("hres", 2.0)
        avg_type = dscfg.get("avg_type", "mean")
        nvalid_min = dscfg.get("nvalid_min", 30)
        interp_kind = dscfg.get("interp_kind", "nearest")
        rmax = dscfg.get("rmax", 50000.0)
        weight_power = dscfg.get("weight_power", 2.0)

        # initialize dataset
        if not dscfg["initialized"]:
            qvp = pyart.retrieve.compute_rqvp(
                radar,
                field_names,
                ref_time=dscfg["timeinfo"],
                hmax=hmax,
                hres=hres,
                avg_type=avg_type,
                nvalid_min=nvalid_min,
                interp_kind=interp_kind,
                rmax=rmax,
                weight_power=weight_power,
                qvp=None,
            )

            if qvp is None:
                warn("Unable to compute QVP")
                return None, None

            global_dict = dict()
            global_dict.update({"start_time": dscfg["timeinfo"]})
            global_dict.update({"radar_out": qvp})
            dscfg["global_data"] = global_dict
            dscfg["initialized"] = 1
        else:
            qvp = pyart.retrieve.compute_rqvp(
                radar,
                field_names,
                ref_time=dscfg["timeinfo"],
                hmax=hmax,
                hres=hres,
                avg_type=avg_type,
                nvalid_min=nvalid_min,
                interp_kind=interp_kind,
                rmax=rmax,
                weight_power=weight_power,
                qvp=dscfg["global_data"]["radar_out"],
            )

            if qvp is None:
                warn("Unable to compute QVP")
                return None, None

            dscfg["global_data"]["radar_out"] = qvp

        new_dataset = dict()
        new_dataset.update({"radar_out": qvp})
        new_dataset.update({"radar_type": "temporal"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg["initialized"]:
            return None, None

        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8]) - 1

        qvp = dscfg["global_data"]["radar_out"]

        new_dataset = dict()
        new_dataset.update({"radar_out": qvp})
        new_dataset.update({"radar_type": "final"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad


def process_evp(procstatus, dscfg, radar_list=None):
    """
    Computes enhanced vertical profiles, by averaging over height levels
    PPI data.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement,
            can be any datatype supported by pyrad and available in the data
        lat, lon : float
            latitude and longitude of the point of interest [deg]
        latlon_tol : float
            tolerance in latitude and longitude in deg. Default 0.0005
        delta_rng, delta_azi : float
            maximum range distance [m] and azimuth distance [degree] from the
            central point of the evp containing data to average. Default 5000.
            and 10.
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 250.
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to consider the data valid when
            performing the averaging. Default 1
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'.
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the EVP and a keyword stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    Reference
    ---------
    Kaltenboeck R., Ryzhkov A. 2016: A freezing rain storm explored with a
    C-band polarimetric weather radar using the QVP methodology.
    Meteorologische Zeitschrift vol. 26 pp 207-222

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        field_names = []
        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8]) - 1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn("ERROR: No valid radar")
            return None, None

        radar = radar_list[ind_rad]

        # default parameters
        lon = dscfg["lon"]
        lat = dscfg["lat"]
        latlon_tol = dscfg.get("latlon_tol", 0.0005)
        delta_rng = dscfg.get("delta_rng", 15000.0)
        delta_azi = dscfg.get("delta_azi", 10.0)
        hmax = dscfg.get("hmax", 10000.0)
        hres = dscfg.get("hres", 250.0)
        avg_type = dscfg.get("avg_type", "mean")
        nvalid_min = dscfg.get("nvalid_min", 1)
        interp_kind = dscfg.get("interp_kind", "none")

        # initialize dataset
        if not dscfg["initialized"]:
            qvp = pyart.retrieve.compute_evp(
                radar,
                field_names,
                lon,
                lat,
                ref_time=dscfg["timeinfo"],
                latlon_tol=latlon_tol,
                delta_rng=delta_rng,
                delta_azi=delta_azi,
                hmax=hmax,
                hres=hres,
                avg_type=avg_type,
                nvalid_min=nvalid_min,
                interp_kind=interp_kind,
                qvp=None,
            )

            if qvp is None:
                warn("Unable to compute QVP")
                return None, None

            global_dict = dict()
            global_dict.update({"start_time": dscfg["timeinfo"]})
            global_dict.update({"radar_out": qvp})
            dscfg["global_data"] = global_dict
            dscfg["initialized"] = 1
        else:
            qvp = pyart.retrieve.compute_evp(
                radar,
                field_names,
                lon,
                lat,
                ref_time=dscfg["timeinfo"],
                latlon_tol=latlon_tol,
                delta_rng=delta_rng,
                delta_azi=delta_azi,
                hmax=hmax,
                hres=hres,
                avg_type=avg_type,
                nvalid_min=nvalid_min,
                interp_kind=interp_kind,
                qvp=dscfg["global_data"]["radar_out"],
            )

            if qvp is None:
                warn("Unable to compute QVP")
                return None, None

        dscfg["global_data"]["radar_out"] = qvp

        new_dataset = dict()
        new_dataset.update({"radar_out": qvp})
        new_dataset.update({"radar_type": "temporal"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg["initialized"]:
            return None, None

        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8]) - 1

        qvp = dscfg["global_data"]["radar_out"]

        new_dataset = dict()
        new_dataset.update({"radar_out": qvp})
        new_dataset.update({"radar_type": "final"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad


def process_svp(procstatus, dscfg, radar_list=None):
    """
    Computes slanted vertical profiles, by averaging over height levels
    PPI data.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement,
            can be any datatype supported by pyrad and available in the data
        angle : int or float
            If the radar object contains a PPI volume, the sweep number to
            use, if it contains an RHI volume the elevation angle.
            Default 0.
        ang_tol : float
            If the radar object contains an RHI volume, the tolerance in the
            elevation angle for the conversion into PPI. Default 1.
        lat, lon : float
            latitude and longitude of the point of interest [deg]
        latlon_tol : float
            tolerance in latitude and longitude in deg. Default 0.0005
        delta_rng, delta_azi : float
            maximum range distance [m] and azimuth distance [degree] from the
            central point of the svp containing data to average. Default 5000.
            and 10.
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 250.
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to consider the data valid when
            performing the averaging. Default 1
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the svp and a keyword stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    Reference
    ---------
    Bukovcic P., Zrnic D., Zhang G. 2017: Winter Precipitation Liquid-Ice
    Phase Transitions Revealed with Polarimetric Radar and 2DVD Observations
    in Central Oklahoma. JTECH vol. 56 pp 1345-1363

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        field_names = []
        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8]) - 1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn("ERROR: No valid radar")
            return None, None

        radar = radar_list[ind_rad]

        # default parameters
        angle = dscfg.get("angle", 0)
        ang_tol = dscfg.get("ang_tol", 1.0)
        lon = dscfg["lon"]
        lat = dscfg["lat"]
        latlon_tol = dscfg.get("latlon_tol", 0.0005)
        delta_rng = dscfg.get("delta_rng", 30000.0)
        delta_azi = dscfg.get("delta_azi", 10.0)
        hmax = dscfg.get("hmax", 10000.0)
        hres = dscfg.get("hres", 250.0)
        avg_type = dscfg.get("avg_type", "mean")
        nvalid_min = dscfg.get("nvalid_min", 1)
        interp_kind = dscfg.get("interp_kind", "none")

        # initialize dataset
        if not dscfg["initialized"]:
            qvp = pyart.retrieve.compute_svp(
                radar,
                field_names,
                lon,
                lat,
                angle,
                ref_time=dscfg["timeinfo"],
                ang_tol=ang_tol,
                latlon_tol=latlon_tol,
                delta_rng=delta_rng,
                delta_azi=delta_azi,
                hmax=hmax,
                hres=hres,
                avg_type=avg_type,
                nvalid_min=nvalid_min,
                interp_kind=interp_kind,
                qvp=None,
            )

            if qvp is None:
                warn("Unable to compute QVP")
                return None, None

            global_dict = dict()
            global_dict.update({"start_time": dscfg["timeinfo"]})
            global_dict.update({"radar_out": qvp})
            dscfg["global_data"] = global_dict
            dscfg["initialized"] = 1
        else:
            qvp = pyart.retrieve.compute_svp(
                radar,
                field_names,
                lon,
                lat,
                angle,
                ref_time=dscfg["timeinfo"],
                ang_tol=ang_tol,
                latlon_tol=latlon_tol,
                delta_rng=delta_rng,
                delta_azi=delta_azi,
                hmax=hmax,
                hres=hres,
                avg_type=avg_type,
                nvalid_min=nvalid_min,
                interp_kind=interp_kind,
                qvp=dscfg["global_data"]["radar_out"],
            )

            if qvp is None:
                warn("Unable to compute QVP")
                return None, None

        dscfg["global_data"]["radar_out"] = qvp

        new_dataset = dict()
        new_dataset.update({"radar_out": qvp})
        new_dataset.update({"radar_type": "temporal"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg["initialized"]:
            return None, None

        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8]) - 1

        svp = dscfg["global_data"]["radar_out"]

        new_dataset = dict()
        new_dataset.update({"radar_out": svp})
        new_dataset.update({"radar_type": "final"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad


def process_time_height(procstatus, dscfg, radar_list=None):
    """
    Produces time height radar objects at a point of interest defined by
    latitude and longitude. A time-height contains the evolution
    of the vertical structure of radar measurements above the location
    of interest.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement,
            can be any datatype supported by pyrad and available in the data
        lat, lon : float
            latitude and longitude of the point of interest [deg]
        latlon_tol : float
            tolerance in latitude and longitude in deg. Default 0.0005
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 50
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the QVP and a keyword stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        field_names = []
        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8]) - 1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn("ERROR: No valid radar")
            return None, None

        radar = radar_list[ind_rad]

        # default parameters
        lon = dscfg["lon"]
        lat = dscfg["lat"]
        latlon_tol = dscfg.get("latlon_tol", 0.0005)
        hmax = dscfg.get("hmax", 10000.0)
        hres = dscfg.get("hres", 50.0)
        interp_kind = dscfg.get("interp_kind", "none")

        # initialize dataset
        if not dscfg["initialized"]:
            qvp = pyart.retrieve.compute_vp(
                radar,
                field_names,
                lon,
                lat,
                ref_time=dscfg["timeinfo"],
                latlon_tol=latlon_tol,
                hmax=hmax,
                hres=hres,
                interp_kind=interp_kind,
                qvp=None,
            )

            if qvp is None:
                warn("Unable to compute QVP")
                return None, None

            global_dict = dict()
            global_dict.update({"start_time": dscfg["timeinfo"]})
            global_dict.update({"radar_out": qvp})
            dscfg["global_data"] = global_dict
            dscfg["initialized"] = 1
        else:
            qvp = pyart.retrieve.compute_vp(
                radar,
                field_names,
                lon,
                lat,
                ref_time=dscfg["timeinfo"],
                latlon_tol=latlon_tol,
                hmax=hmax,
                hres=hres,
                interp_kind=interp_kind,
                qvp=dscfg["global_data"]["radar_out"],
            )

            if qvp is None:
                warn("Unable to compute QVP")
                return None, None

        dscfg["global_data"]["radar_out"] = qvp

        new_dataset = dict()
        new_dataset.update({"radar_out": qvp})
        new_dataset.update({"radar_type": "temporal"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg["initialized"]:
            return None, None

        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8]) - 1

        qvp = dscfg["global_data"]["radar_out"]

        new_dataset = dict()
        new_dataset.update({"radar_out": qvp})
        new_dataset.update({"radar_type": "final"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad


def process_ts_along_coord(procstatus, dscfg, radar_list=None):
    """
    Produces time series along a particular antenna coordinate

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the time series,
            can be any datatype supported by pyrad and available in the data
        mode : str
            coordinate to extract data along. Can be ALONG_AZI, ALONG_ELE or
            ALONG_RNG
        fixed_range, fixed_azimuth, fixed_elevation : float
            The fixed range [m], azimuth [deg] or elevation [deg] to extract.
            In each mode two of these parameters have to be defined. If they
            are not defined they default to 0.
        ang_tol, rng_tol : float
            The angle tolerance [deg] and range tolerance [m] around the fixed
            range or azimuth/elevation
        value_start, value_stop : float
            The minimum and maximum value at which the data along a coordinate
            start and stop

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the data and a keyword stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
        field_names = []
        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            field_names.append(get_fieldname_pyart(datatype))

        ind_rad = int(radarnr[5:8]) - 1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn("ERROR: No valid radar")
            return None, None

        radar = radar_list[ind_rad]

        mode = dscfg.get("mode", "ALONG_AZI")
        if mode == "ALONG_RNG":
            value_start = dscfg.get("value_start", 0.0)
            value_stop = dscfg.get("value_stop", radar.range["data"][-1])
            ang_tol = dscfg.get("AngTol", 1.0)
            rng_tol = dscfg.get("RngTol", 50.0)
            fixed_range = dscfg.get("fixed_range", None)
            fixed_azimuth = dscfg.get("fixed_azimuth", 0.0)
            fixed_elevation = dscfg.get("fixed_elevation", 0.0)
        elif mode == "ALONG_AZI":
            value_start = dscfg.get("value_start", np.min(radar.azimuth["data"]))
            value_stop = dscfg.get("value_stop", np.max(radar.azimuth["data"]))
            ang_tol = dscfg.get("AngTol", 1.0)
            rng_tol = dscfg.get("RngTol", 50.0)
            fixed_range = dscfg.get("fixed_range", 0.0)
            fixed_azimuth = dscfg.get("fixed_azimuth", None)
            fixed_elevation = dscfg.get("fixed_elevation", 0.0)
        elif mode == "ALONG_ELE":
            value_start = dscfg.get("value_start", np.min(radar.elevation["data"]))
            value_stop = dscfg.get("value_stop", np.max(radar.elevation["data"]))
            ang_tol = dscfg.get("AngTol", 1.0)
            rng_tol = dscfg.get("RngTol", 50.0)
            fixed_range = dscfg.get("fixed_range", 0.0)
            fixed_azimuth = dscfg.get("fixed_azimuth", 0.0)
            fixed_elevation = dscfg.get("fixed_elevation", None)
        else:
            warn("Unknown plotting mode " + dscfg["mode"])
            return None, None

        # initialize dataset
        if not dscfg["initialized"]:
            acoord = pyart.retrieve.compute_ts_along_coord(
                radar,
                field_names,
                mode=mode,
                fixed_range=fixed_range,
                fixed_azimuth=fixed_azimuth,
                fixed_elevation=fixed_elevation,
                ang_tol=ang_tol,
                rng_tol=rng_tol,
                value_start=value_start,
                value_stop=value_stop,
                ref_time=dscfg["timeinfo"],
                acoord=None,
            )

            if acoord is None:
                warn("Unable to compute time series along coordinate")
                return None, None

            global_dict = dict()
            global_dict.update({"start_time": dscfg["timeinfo"]})
            global_dict.update({"radar_out": acoord})
            dscfg["global_data"] = global_dict
            dscfg["initialized"] = 1
        else:
            acoord = pyart.retrieve.compute_ts_along_coord(
                radar,
                field_names,
                mode=mode,
                fixed_range=fixed_range,
                fixed_azimuth=fixed_azimuth,
                fixed_elevation=fixed_elevation,
                ang_tol=ang_tol,
                rng_tol=rng_tol,
                value_start=value_start,
                value_stop=value_stop,
                ref_time=dscfg["timeinfo"],
                acoord=dscfg["global_data"]["radar_out"],
            )

            if acoord is None:
                warn("Unable to compute time series along coordinate")
                return None, None

        dscfg["global_data"]["radar_out"] = acoord

        new_dataset = dict()
        new_dataset.update({"radar_out": acoord})
        new_dataset.update({"radar_type": "temporal"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad

    if procstatus == 2:
        if not dscfg["initialized"]:
            return None, None

        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8]) - 1

        acoord = dscfg["global_data"]["radar_out"]

        new_dataset = dict()
        new_dataset.update({"radar_out": acoord})
        new_dataset.update({"radar_type": "final"})
        new_dataset.update({"start_time": dscfg["global_data"]["start_time"]})

        return new_dataset, ind_rad
