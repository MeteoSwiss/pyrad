"""
pyrad.proc.process_icon
===========================

Functions to manage icon data

.. autosummary::
    :toctree: generated/

    process_icon
    process_hzt
    process_iso0_mf
    process_iso0_grib
    process_icon_lookup_table
    process_hzt_lookup_table
    process_icon_to_radar
    process_icon_coord
    process_hzt_coord

"""

from copy import deepcopy
from warnings import warn
import glob
import os

import numpy as np
from netCDF4 import num2date

import pyart

from ..io.io_aux import get_datatype_fields, find_raw_icon_file
from ..io.io_aux import find_hzt_file, find_iso0_file, find_iso0_grib_file
from ..io.io_aux import get_fieldname_pyart
from ..io.read_data_icon import read_icon_data, read_icon_coord
from ..io.read_data_icon import icon2radar_data, icon2radar_coord
from ..io.read_data_icon import get_icon_fields
from ..io.read_data_radar import interpol_field
from ..io.read_data_hzt import read_hzt_data, hzt2radar_data, hzt2radar_coord
from ..io.read_data_hzt import get_iso0_field
from ..io.read_data_iso0_mf import read_iso0_mf_data, read_iso0_grib_data
from ..io.read_data_iso0_mf import iso2radar_data, grib2radar_data

# from memory_profiler import profile


def process_icon(procstatus, dscfg, radar_list=None):
    """
    Gets icon data and put it in radar coordinates

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            arbitrary data type supported by pyrad
        keep_in_memory : int. Dataset keyword
            if set keeps the icon data dict, the icon coordinates dict and
            the icon field in radar coordinates in memory
        regular_grid : int. Dataset keyword
            if set it is assume that the radar has a grid constant in time and
            there is no need to compute a new icon field if the icon
            data has not changed
        icon_type : str. Dataset keyword
            name of the icon field to process. Default TEMP
        icon_variables : list of strings. Dataset keyword
            Py-art name of the icon fields. Default temperature
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field corresponding to
        icon_variables
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # debugging
    # start_time = time.time()

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    keep_in_memory = dscfg.get("keep_in_memory", 0)
    regular_grid = dscfg.get("regular_grid", 0)

    icon_type = dscfg.get("icon_type", "TEMP")
    if icon_type == "TEMP":
        field_names = ["temperature"]
        if "icon_variables" in dscfg:
            field_names = []
            for var in dscfg["icon_variables"]:
                field_names.append(get_fieldname_pyart(var))
        zmin = None
    elif icon_type == "WIND":
        field_names = ["wind_speed", "wind_direction", "vertical_wind_shear"]
        if "icon_variables" in dscfg:
            field_names = []
            for var in dscfg["icon_variables"]:
                field_names.append(get_fieldname_pyart(var))
        zmin = None
    else:
        warn("Unknown icon data type " + icon_type)
        return None, None

    fname = find_raw_icon_file(dscfg["timeinfo"], icon_type, dscfg, ind_rad=ind_rad)

    if fname is None:
        return None, None

    model = os.path.basename(fname)[14:26]
    if model not in ("icon-ch1-eps", "icon-ch2-eps"):
        warn("Unknown NWP model " + model)
        return None, None

    if keep_in_memory:
        if dscfg["initialized"] == 0:
            icon_coord = read_icon_coord(
                dscfg["iconpath"][ind_rad] + "rad2icon/" + model + "_MDR_3D_const.nc",
                zmin=zmin,
            )
            dscfg["global_data"] = {
                "icon_fname": None,
                "icon_data": None,
                "icon_fields": None,
                "time_index": None,
                "icon_coord": icon_coord,
            }
            dscfg["initialized"] = 1

        icon_coord = dscfg["global_data"]["icon_coord"]
        if fname != dscfg["global_data"]["icon_fname"]:
            # debugging
            # start_time2 = time.time()
            icon_data = read_icon_data(fname, field_names=field_names, celsius=True)
            # print(" reading icon takes %s seconds " %
            #      (time.time() - start_time2))
            if icon_data is None:
                warn("icon data not found")
                return None, None

            dscfg["global_data"]["icon_data"] = icon_data
            dscfg["global_data"]["icon_fname"] = fname
        else:
            print("raw icon data already in memory")
            icon_data = dscfg["global_data"]["icon_data"]
    else:
        icon_coord = read_icon_coord(
            dscfg["iconpath"][ind_rad] + "rad2icon/" + model + "_MDR_3D_const.nc",
            zmin=zmin,
        )

        # debugging
        # start_time2 = time.time()
        icon_data = read_icon_data(fname, field_names=field_names, celsius=True)
        # print(" reading icon takes %s seconds " %
        #      (time.time() - start_time2))
        if icon_data is None:
            warn("icon data not found")
            return None, None

    dticon = num2date(icon_data["time"]["data"][:], icon_data["time"]["units"])
    time_index = np.argmin(abs(dticon - dscfg["timeinfo"]))

    if keep_in_memory and regular_grid:
        if time_index != dscfg["global_data"]["time_index"]:
            icon_fields = icon2radar_data(
                radar,
                icon_coord,
                icon_data,
                time_index=time_index,
                field_names=field_names,
            )
            if icon_fields is None:
                warn("Unable to obtain icon fields")
                return None, None

            dscfg["global_data"]["time_index"] = time_index
            dscfg["global_data"]["icon_fields"] = icon_fields
        else:
            print("icon field already in memory")
            icon_fields = dscfg["global_data"]["icon_fields"]
    else:
        icon_fields = icon2radar_data(
            radar, icon_coord, icon_data, time_index=time_index, field_names=field_names
        )
        if icon_fields is None:
            warn("Unable to obtain icon fields")
            return None, None

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()

    for field in icon_fields:
        for field_name in field:
            new_dataset["radar_out"].add_field(field_name, field[field_name])

    return new_dataset, ind_rad


def process_hzt(procstatus, dscfg, radar_list=None):
    """
    Gets iso0 degree data in HZT format and put it in radar coordinates

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        metranet_read_lib : str. Global keyword
            Type of METRANET reader library used to read the data.
            Can be 'C' or 'python'
        datatype : string. Dataset keyword
            arbitrary data type supported by pyrad
        keep_in_memory : int. Dataset keyword
            if set keeps the icon data dict, the icon coordinates dict and
            the icon field in radar coordinates in memory
        regular_grid : int. Dataset keyword
            if set it is assume that the radar has a grid constant in time and
            there is no need to compute a new icon field if the icon
            data has not changed
        icon_type : str. Dataset keyword
            name of the icon field to process. Default TEMP
        icon_variables : list of strings. Dataset keyword
            Py-art name of the icon fields. Default temperature
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields corresponding to
        icon_variables
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # debugging
    # start_time = time.time()

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    keep_in_memory = dscfg.get("keep_in_memory", 0)
    regular_grid = dscfg.get("regular_grid", 0)

    fname = find_hzt_file(dscfg["timeinfo"], dscfg, ind_rad=ind_rad)

    if fname is None:
        return None, None

    if keep_in_memory:
        if dscfg["initialized"] == 0:
            hzt_data = read_hzt_data(fname, read_lib=dscfg["metranet_read_lib"])
            if hzt_data is None:
                warn("HZT data not found")
                return None, None
            hzt_coord = {"x": hzt_data["x"], "y": hzt_data["y"]}
            dscfg["global_data"] = {
                "hzt_fname": None,
                "hzt_data": None,
                "iso0_field": None,
                "time_index": None,
                "hzt_coord": hzt_coord,
            }
            dscfg["initialized"] = 1

        hzt_coord = dscfg["global_data"]["hzt_coord"]
        if fname != dscfg["global_data"]["hzt_fname"]:
            hzt_data = read_hzt_data(fname, read_lib=dscfg["metranet_read_lib"])
            if hzt_data is None:
                warn("HZT data not found")
                return None, None

            dscfg["global_data"]["hzt_data"] = hzt_data
            dscfg["global_data"]["hzt_fname"] = fname
        else:
            print("raw HZT data already in memory")
            hzt_data = dscfg["global_data"]["hzt_data"]
    else:
        hzt_data = read_hzt_data(fname, read_lib=dscfg["metranet_read_lib"])
        if hzt_data is None:
            warn("HZT data not found")
            return None, None
        hzt_coord = {"x": hzt_data["x"], "y": hzt_data["y"]}

    dthzt = num2date(hzt_data["time"]["data"][:], hzt_data["time"]["units"])
    time_index = np.argmin(abs(dthzt - dscfg["timeinfo"]))

    if keep_in_memory and regular_grid:
        if time_index != dscfg["global_data"]["time_index"]:
            iso0_field = hzt2radar_data(radar, hzt_coord, hzt_data)
            if iso0_field is None:
                warn("Unable to obtain heigth over iso0 field")
                return None, None

            dscfg["global_data"]["time_index"] = time_index
            dscfg["global_data"]["iso0_field"] = iso0_field
        else:
            print("HZT field already in memory")
            iso0_field = dscfg["global_data"]["iso0_field"]
    else:
        iso0_field = hzt2radar_data(radar, hzt_coord, hzt_data)
        if iso0_field is None:
            warn("Unable to obtain HZT fields")
            return None, None

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field("height_over_iso0", iso0_field)

    return new_dataset, ind_rad


def process_iso0_mf(procstatus, dscfg, radar_list=None):
    """
    Gets iso0 degree data in text format and put it in radar coordinates.
    This function is meant to process data received from the MeteoFrance NWP
    model. The model provides a maximum of 9 points at 0.5 degree lat/lon
    spacing surrounding a given radar. If a point is not provided it means
    that the iso0 for that point is at or below the ground level. Out of these
    points a single reference iso-0 is obtained according to the user defined
    iso0 statistic.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            arbitrary data type supported by prad
        iso0_statistic : str. Dataset keyword
            The statistic used to weight the iso0 points. Can be avg_by_dist,
            avg, min, max

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "H_ISO0"
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # debugging
    # start_time = time.time()

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    iso0_statistic = dscfg.get("iso0_statistic", "avg_by_dist")
    if iso0_statistic not in ("avg_by_dist", "avg", "min", "max"):
        warn(
            "iso0 statistic "
            + iso0_statistic
            + " unknown. Default avg_by_dist will be applied"
        )
        iso0_statistic = "avg_by_dist"

    fname = find_iso0_file(dscfg["timeinfo"], dscfg, ind_rad=ind_rad)
    if fname is None:
        return None, None

    iso0_data = read_iso0_mf_data(fname)
    if iso0_data is None:
        warn("iso0 data not found")
        return None, None

    np.argmin(abs(iso0_data["fcst_time"] - dscfg["timeinfo"]))
    iso0_field = iso2radar_data(
        radar, iso0_data, dscfg["timeinfo"], iso0_statistic=iso0_statistic
    )
    if iso0_field is None:
        warn("Unable to obtain iso0 fields")
        return None, None

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field("height_over_iso0", iso0_field)

    return new_dataset, ind_rad


def process_iso0_grib(procstatus, dscfg, radar_list=None):
    """
    Gets iso0 degree data in GRIB format and put it in radar coordinates.
    This function is meant to process data received from the MeteoFrance NWP
    model. It can output the height over the iso0 of each gate or the iso0
    height at each gate

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            arbitrary data type supported by pyrad
        time_interp : bool. Dataset keyword
            whether to perform an interpolation in time between consecutive
            model outputs. Default True
        voltype: str. Dataset keyword
            The type of data to output. Can be H_ISO0 or HZT. Default H_ISO0

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field H_ISO0
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # debugging
    # start_time = time.time()

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    time_interp = dscfg.get("time_interp", True)
    voltype = dscfg.get("voltype", "H_ISO0")
    field_name = get_fieldname_pyart(voltype)

    fname = find_iso0_grib_file(dscfg["timeinfo"], dscfg, ind_rad=ind_rad)
    if fname is None:
        return None, None

    iso0_data = read_iso0_grib_data(fname)
    if iso0_data is None:
        warn("iso0 data not found")
        return None, None

    np.argmin(abs(iso0_data["fcst_time"] - dscfg["timeinfo"]))
    iso0_field = grib2radar_data(
        radar,
        iso0_data,
        dscfg["timeinfo"],
        time_interp=time_interp,
        field_name=field_name,
    )
    if iso0_field is None:
        warn("Unable to obtain iso0 fields")
        return None, None

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()
    new_dataset["radar_out"].add_field(field_name, iso0_field)

    return new_dataset, ind_rad


# @profile
def process_icon_lookup_table(procstatus, dscfg, radar_list=None):
    """
    Gets icon data and put it in radar coordinates
    using look up tables computed or loaded when initializing

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            arbitrary data type supported by pyrad
        lookup_table : int. Dataset keyword
            if set a pre-computed look up table for the icon coordinates is
            loaded. Otherwise the look up table is computed taking the first
            radar object as reference
        regular_grid : int. Dataset keyword
            if set it is assume that the radar has a grid constant in time and
            therefore there is no need to interpolate the icon field in
            memory to the current radar grid
        icon_type : str. Dataset keyword
            name of the icon field to process. Default TEMP
        icon_variables : list of strings. Dataset keyword
            Py-art name of the icon fields. Default temperature
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields corresponding to icon_variables
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # debugging
    # start_time = time.time()

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    regular_grid = dscfg.get("regular_grid", 0)
    lookup_table = dscfg.get("lookup_table", 0)

    icon_type = dscfg.get("icon_type", "TEMP")
    if icon_type == "TEMP":
        field_names = ["temperature"]
        if "icon_variables" in dscfg:
            field_names = []
            for var in dscfg["icon_variables"]:
                field_names.append(get_fieldname_pyart(var))
        zmin = None
    elif icon_type == "WIND":
        field_names = ["wind_speed", "wind_direction", "vertical_wind_shear"]
        if "icon_variables" in dscfg:
            field_names = []
            for var in dscfg["icon_variables"]:
                field_names.append(get_fieldname_pyart(var))
        zmin = None
    else:
        warn("Unknown icon data type " + icon_type)
        return None, None

    fname = find_raw_icon_file(dscfg["timeinfo"], icon_type, dscfg, ind_rad=ind_rad)

    if fname is None:
        return None, None

    model = os.path.basename(fname)[14:26]
    if model not in ("icon-ch1-eps", "icon-ch2-eps"):
        warn("Unknown NWP model " + model)
        return None, None

    if dscfg["initialized"] == 0:
        if lookup_table:
            savedir = dscfg["iconpath"][ind_rad] + "rad2icon/"
            fname_ind = "rad2icon_icon_index_" + dscfg["procname"] + ".nc"
            fname_ind2 = glob.glob(savedir + fname_ind)
            if not fname_ind2:
                warn("File " + savedir + fname_ind + " not found")
                return None, None
            icon_radar = pyart.io.read_cfradial(fname_ind2[0])
        else:
            icon_coord = read_icon_coord(
                dscfg["iconpath"][ind_rad] + "rad2icon/" + model + "_MDR_3D_const.nc",
                zmin=zmin,
            )
            print("icon coordinates files read")
            icon_ind_field = icon2radar_coord(radar, icon_coord)
            icon_radar = deepcopy(radar)
            icon_radar.fields = dict()
            icon_radar.add_field("icon_index", icon_ind_field)
            print("icon index field added")

        dscfg["global_data"] = {
            "icon_fname": None,
            "icon_data": None,
            "icon_fields": None,
            "icon_radar": icon_radar,
            "time_index": None,
        }
        dscfg["initialized"] = 1

    if fname != dscfg["global_data"]["icon_fname"]:
        # debugging
        # start_time2 = time.time()
        icon_data = read_icon_data(fname, field_names=field_names, celsius=True)
        # print(" reading icon takes %s seconds " %
        #      (time.time() - start_time2))
        if icon_data is None:
            warn("icon data not found")
            return None, None

        dscfg["global_data"]["icon_data"] = icon_data
    else:
        print("raw icon data already in memory")
        icon_data = dscfg["global_data"]["icon_data"]

    dticon = num2date(icon_data["time"]["data"][:], icon_data["time"]["units"])

    time_index = np.argmin(abs(dticon - dscfg["timeinfo"]))

    if (
        fname != dscfg["global_data"]["icon_fname"]
        or time_index != dscfg["global_data"]["time_index"]
    ):
        # debugging
        # start_time3 = time.time()
        icon_fields = get_icon_fields(
            icon_data,
            dscfg["global_data"]["icon_radar"].fields["icon_index"],
            time_index=time_index,
            field_names=field_names,
        )
        if icon_fields is None:
            warn("Unable to obtain icon fields")
            return None, None
        # print(" getting icon data takes %s seconds "
        #      % (time.time() - start_time3))

        dscfg["global_data"]["time_index"] = time_index
        dscfg["global_data"]["icon_fields"] = icon_fields
    else:
        print("icon field already in memory")
        icon_fields = dscfg["global_data"]["icon_fields"]

    dscfg["global_data"]["icon_fname"] = fname

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()

    if not regular_grid:
        radar_aux = deepcopy(dscfg["global_data"]["icon_radar"])
        radar_aux.fields = dict()

    for field in icon_fields:
        for field_name in field:
            try:
                if regular_grid:
                    new_dataset["radar_out"].add_field(field_name, field[field_name])
                else:
                    # interpolate to current radar grid
                    radar_aux.add_field(field_name, field[field_name])
                    icon_field_interp = interpol_field(radar, radar_aux, field_name)
                    new_dataset["radar_out"].add_field(field_name, icon_field_interp)
            except Exception as ee:
                warn(str(ee))
                warn("Unable to add icon " + field_name + " field to radar object")
                return None, None

    return new_dataset, ind_rad


def process_hzt_lookup_table(procstatus, dscfg, radar_list=None):
    """
    Gets HZT data and put it in radar coordinates
    using look up tables computed or loaded when initializing

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        metranet_read_lib : str. Global keyword
            Type of METRANET reader library used to read the data.
            Can be 'C' or 'python'
        datatype : string. Dataset keyword
            arbitrary data type supported by pyrad
        lookup_table : int. Dataset keyword
            if set a pre-computed look up table for the icon coordinates is
            loaded. Otherwise the look up table is computed taking the first
            radar object as reference
        regular_grid : int. Dataset keyword
            if set it is assume that the radar has a grid constant in time and
            therefore there is no need to interpolate the icon field in
            memory to the current radar grid
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field HISO0
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # debugging
    # start_time = time.time()

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    regular_grid = dscfg.get("regular_grid", 0)
    lookup_table = dscfg.get("lookup_table", 0)

    fname = find_hzt_file(dscfg["timeinfo"], dscfg, ind_rad=ind_rad)

    if fname is None:
        return None, None

    if dscfg["initialized"] == 0:
        if lookup_table:
            savedir = dscfg["iconpath"][ind_rad] + "rad2icon/"
            fname_ind = "rad2icon_hzt_index_" + dscfg["procname"] + ".nc"
            fname_ind2 = glob.glob(savedir + fname_ind)
            if not fname_ind2:
                warn("File " + savedir + fname_ind + " not found")
                return None, None
            hzt_radar = pyart.io.read_cfradial(fname_ind2[0])
        else:
            hzt_data = read_hzt_data(fname, read_lib=dscfg["metranet_read_lib"])
            if hzt_data is None:
                warn("HZT data not found")
                return None, None
            hzt_coord = {"x": hzt_data["x"], "y": hzt_data["y"]}
            hzt_ind_field = hzt2radar_coord(radar, hzt_coord)
            hzt_radar = deepcopy(radar)
            hzt_radar.fields = dict()
            hzt_radar.add_field("hzt_index", hzt_ind_field)

        dscfg["global_data"] = {
            "hzt_fname": None,
            "hzt_data": None,
            "iso0_field": None,
            "hzt_radar": hzt_radar,
            "time_index": None,
        }
        dscfg["initialized"] = 1

    if fname != dscfg["global_data"]["hzt_fname"]:
        hzt_data = read_hzt_data(fname, read_lib=dscfg["metranet_read_lib"])
        if hzt_data is None:
            warn("HZT data not found")
            return None, None

        dscfg["global_data"]["hzt_data"] = hzt_data
    else:
        print("HZT data already in memory")
        hzt_data = dscfg["global_data"]["hzt_data"]

    dthzt = num2date(hzt_data["time"]["data"][:], hzt_data["time"]["units"])
    time_index = np.argmin(abs(dthzt - dscfg["timeinfo"]))

    if (
        fname != dscfg["global_data"]["hzt_fname"]
        or time_index != dscfg["global_data"]["time_index"]
    ):
        # debugging
        # start_time3 = time.time()
        iso0_field = get_iso0_field(
            hzt_data,
            dscfg["global_data"]["hzt_radar"].fields["hzt_index"],
            dscfg["global_data"]["hzt_radar"].gate_altitude["data"],
        )
        if iso0_field is None:
            warn("Unable to obtain height over iso0")
            return None, None

        dscfg["global_data"]["time_index"] = time_index
        dscfg["global_data"]["iso0_field"] = iso0_field
    else:
        print("HZT field already in memory")
        iso0_field = dscfg["global_data"]["iso0_field"]

    dscfg["global_data"]["hzt_fname"] = fname

    # prepare for exit
    new_dataset = {"radar_out": deepcopy(radar)}
    new_dataset["radar_out"].fields = dict()

    if not regular_grid:
        radar_aux = deepcopy(dscfg["global_data"]["hzt_radar"])
        radar_aux.fields = dict()

    try:
        if regular_grid:
            new_dataset["radar_out"].add_field("height_over_iso0", iso0_field)
        else:
            # interpolate to current radar grid
            radar_aux.add_field("height_over_iso0", iso0_field)
            hzt_field_interp = interpol_field(radar, radar_aux, "height_over_iso0")
            new_dataset["radar_out"].add_field("height_over_iso0", hzt_field_interp)
    except Exception as ee:
        warn(str(ee))
        warn("Unable to add height_over_iso0 " + " field to radar object")
        return None, None

    return new_dataset, ind_rad


# @profile
def process_icon_to_radar(procstatus, dscfg, radar_list=None):
    """
    Gets icon data and put it in radar coordinates using look up tables

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            arbitrary data type supported by pyrad
        icon_type : str. Dataset keyword
            name of the icon field to process. Default TEMP
        icon_variables : list of strings. Dataset keyword
            Py-art name of the icon fields. Default temperature
        icon_time_index_min, icon_time_index_max : int
            minimum and maximum indices of the icon data to retrieve. If a
            value is provided only data corresponding to the time indices
            within the interval will be used. If None all data will be used.
            Default None
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output fields corresponding to icon_variables
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break

    ind_rad = int(radarnr[5:8]) - 1

    icon_type = dscfg.get("icon_type", "TEMP")
    if icon_type == "TEMP":
        field_names = ["temperature"]
        if "icon_variables" in dscfg:
            field_names = []
            for var in dscfg["icon_variables"]:
                field_names.append(get_fieldname_pyart(var))
    elif icon_type == "WIND":
        field_names = ["wind_speed", "wind_direction", "vertical_wind_shear"]
        if "icon_variables" in dscfg:
            field_names = []
            for var in dscfg["icon_variables"]:
                field_names.append(get_fieldname_pyart(var))
    else:
        warn("Unknown icon data type " + icon_type)
        return None, None

    time_index_min = dscfg.get("icon_time_index_min", None)
    time_index_max = dscfg.get("icon_time_index_max", None)

    fname = find_raw_icon_file(dscfg["timeinfo"], icon_type, dscfg, ind_rad=ind_rad)

    if fname is None:
        return None, None

    model = os.path.basename(fname)[14:26]
    if model not in ("icon-ch1-eps", "icon-ch2-eps"):
        warn("Unknown NWP model " + model)
        return None, None

    if dscfg["initialized"] == 0:
        savedir = dscfg["iconpath"][ind_rad] + "rad2icon/"
        fname_ind = "rad2icon_icon_index_" + dscfg["procname"] + ".nc"
        fname_ind2 = glob.glob(savedir + fname_ind)
        if not fname_ind2:
            warn("File " + savedir + fname_ind + " not found")
            return None, None
        icon_radar = pyart.io.read_cfradial(fname_ind2[0])

        dscfg["global_data"] = {"icon_radar": icon_radar}
        dscfg["initialized"] = 1

    icon_data = read_icon_data(fname, field_names=field_names, celsius=True)
    if icon_data is None:
        warn("icon data not found")
        return None, None

    dticon = num2date(icon_data["time"]["data"][:], icon_data["time"]["units"])

    if time_index_min is None:
        time_index_min = 0
    if time_index_max is None:
        time_index_max = len(dticon) - 1

    if time_index_max > len(dticon) - 1:
        warn(
            "icon_time_index_max larger than number of icon forecasts"
            "icon forecasts " + len(dticon)
        )
        time_index_max = len(dticon) - 1

    icon_radars = []
    for time_index, dtc in enumerate(dticon):
        if time_index < time_index_min:
            continue
        if time_index > time_index_max:
            break

        icon_fields = get_icon_fields(
            icon_data,
            dscfg["global_data"]["icon_radar"].fields["icon_index"],
            time_index=time_index,
            field_names=field_names,
        )
        if icon_fields is None:
            warn("Unable to obtain icon fields")
            return None, None

        radar_out = deepcopy(dscfg["global_data"]["icon_radar"])
        radar_out.fields = dict()
        radar_out.time["units"] = dtc.strftime("seconds since %Y-%m-%dT%H:%M:%SZ")

        for field in icon_fields:
            for field_name in field:
                radar_out.add_field(field_name, field[field_name])

        icon_radars.append({"ind_rad": ind_rad, "radar_out": radar_out, "dticon": dtc})

    return icon_radars, ind_rad


def process_icon_coord(procstatus, dscfg, radar_list=None):
    """
    Gets the icon indices corresponding to each icon coordinates

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            arbitrary data type supported by pyrad
        iconpath : string. General keyword
            path where to store the look up table
        model : string. Dataset keyword
            The icon model to use. Can be icon-1, icon-1e, icon-2, icon-7
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field
        "icon_index"
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    if dscfg["initialized"] == 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    model = dscfg.get("model", "icon-ch1-eps")
    if model not in ("icon-ch1-eps", "icon-ch2-eps"):
        warn("Unknown NWP model " + model)
        return None, None

    icon_coord = read_icon_coord(
        dscfg["iconpath"][ind_rad] + "rad2icon/" + model + "_MDR_3D_const.nc", zmin=None
    )

    if icon_coord is None:
        return None, None

    icon_ind_field = icon2radar_coord(radar, icon_coord, slice_xy=True, slice_z=False)

    # prepare for exit
    radar_obj = deepcopy(radar)
    radar_obj.fields = dict()
    radar_obj.add_field("icon_index", icon_ind_field)

    new_dataset = {"ind_rad": ind_rad, "radar_out": radar_obj}

    dscfg["initialized"] = 1

    return new_dataset, ind_rad


def process_hzt_coord(procstatus, dscfg, radar_list=None):
    """
    Gets the HZT indices corresponding to each HZT coordinates

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        metranet_read_lib : str. Global keyword
            Type of METRANET reader library used to read the data.
            Can be 'C' or 'python'
        datatype : string. Dataset keyword
            arbitrary data type supported by pyrad
        iconpath : string. General keyword
            path where to store the look up table
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output field "icon_index"
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    if dscfg["initialized"] == 1:
        return None, None

    for datatypedescr in dscfg["datatype"]:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break

    ind_rad = int(radarnr[5:8]) - 1
    if radar_list[ind_rad] is None:
        warn("No valid radar")
        return None, None
    radar = radar_list[ind_rad]

    fname = find_hzt_file(dscfg["timeinfo"], dscfg, ind_rad=ind_rad)

    hzt_data = read_hzt_data(fname, read_lib=dscfg["metranet_read_lib"])
    if hzt_data is None:
        warn("HZT data not found")
        return None, None
    hzt_coord = {"x": hzt_data["x"], "y": hzt_data["y"]}

    hzt_ind_field = hzt2radar_coord(radar, hzt_coord)

    # prepare for exit
    radar_obj = deepcopy(radar)
    radar_obj.fields = dict()
    radar_obj.add_field("hzt_index", hzt_ind_field)

    new_dataset = {"ind_rad": ind_rad, "radar_out": radar_obj}

    dscfg["initialized"] = 1

    return new_dataset, ind_rad
