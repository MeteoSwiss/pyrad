"""
pyrad.proc.process_intercomp
============================

Functions used in the inter-comparison between radars

.. autosummary::
    :toctree: generated/

    process_colocated_gates
    process_intercomp
    process_intercomp_with_QC
    process_fields_diff
    process_intercomp_fields

"""

from copy import deepcopy
from ..util import warn
import numpy as np
import os
from netCDF4 import num2date

import pyart

from ..util.date_utils import cftodatetime
from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename
from ..io.read_data_other import read_colocated_gates, read_colocated_data
from ..io.read_data_other import read_colocated_data_with_QC

from ..util.radar_utils import get_range_bins_to_avg
from ..util.radar_utils import find_colocated_indexes_fast


def process_colocated_gates(procstatus, dscfg, radar_list=None):
    """
    Find colocated gates within two radars

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types to use to check colocated gates (one for every radar)
            Any datatype supported by pyrad and available in both radars is accepted.
            If visibility filtering is desired, the fields
            "visibility" or "visibility_polar" must be specified for both radars.
        h_tol : float. Dataset keyword
            Tolerance in altitude difference between radar gates [m].
            Default 100.
        latlon_tol : float. Dataset keyword
            Tolerance in latitude and longitude position between radar gates
            [deg]. Default 0.0005
        vol_d_tol : float. Dataset keyword
            Tolerance in pulse volume diameter [m]. Default 100.
        vismin : float. Dataset keyword
            Minimum visibility [percent]. Default None.
        hmin : float. Dataset keyword
            Minimum altitude [m MSL]. Default None.
        hmax : float. Dataset keyword
            Maximum altitude [m MSL]. Default None.
        rmin : float. Dataset keyword
            Minimum range [m]. Default None.
        rmax : float. Dataset keyword
            Maximum range [m]. Default None.
        elmin : float. Dataset keyword
            Minimum elevation angle [deg]. Default None.
        elmax : float. Dataset keyword
            Maximum elevation angle [deg]. Default None.
        azrad1min : float. Dataset keyword
            Minimum azimuth angle [deg] for radar 1. Default None.
        azrad1max : float. Dataset keyword
            Maximum azimuth angle [deg] for radar 1. Default None.
        azrad2min : float. Dataset keyword
            Minimum azimuth angle [deg] for radar 2. Default None.
        azrad2max : float. Dataset keyword
            Maximum azimuth angle [deg] for radar 2. Default None.

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the field "colocated_gates"
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # check how many radars are there
    radarnr_dict = dict()
    ind_radar_list = set()
    for datatypedescr in dscfg["datatype"]:
        radarnr = datatypedescr.split(":")[0]
        radarnr_dict.update({radarnr: []})
        ind_radar_list.add(int(radarnr[5:8]) - 1)

    ind_radar_list = list(ind_radar_list)

    if (len(radarnr_dict) != 2) or (len(radar_list) < 2):
        warn("Intercomparison requires data from two different radars")
        return None, None

    # create the list of data types for each radar
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if radarnr in radarnr_dict:
            radarnr_dict[radarnr].append(get_fieldname_pyart(datatype))

    radar1 = radar_list[ind_radar_list[0]]
    radar2 = radar_list[ind_radar_list[1]]

    if radar1 is None or radar2 is None:
        warn("Unable to inter-compare radars. Missing radar")

    if "instrument_name" in radar1.metadata:
        print("Radar 1: " + radar1.metadata["instrument_name"])
    if "instrument_name" in radar2.metadata:
        print("Radar 2: " + radar2.metadata["instrument_name"])

    coloc_gates_field = "colocated_gates"

    h_tol = dscfg.get("h_tol", 100.0)
    latlon_tol = dscfg.get("latlon_tol", 0.0005)
    vol_d_tol = dscfg.get("vol_d_tol", 100.0)
    vismin = dscfg.get("vismin", None)
    hmin = dscfg.get("hmin", None)
    hmax = dscfg.get("hmax", None)
    rmin = dscfg.get("rmin", None)
    rmax = dscfg.get("rmax", None)
    elmin = dscfg.get("elmin", None)
    elmax = dscfg.get("elmax", None)
    azrad1min = dscfg.get("azrad1min", None)
    azrad1max = dscfg.get("azrad1max", None)
    azrad2min = dscfg.get("azrad2min", None)
    azrad2max = dscfg.get("azrad2max", None)

    visib_field = None
    if "visibility" in radarnr_dict["RADAR" + "{:03d}".format(ind_radar_list[0] + 1)]:
        visib_field = "visibility"  # older IDL visibility codes
    elif (
        "visibility_polar"
        in radarnr_dict["RADAR" + "{:03d}".format(ind_radar_list[0] + 1)]
    ):
        visib_field = "visibility_polar"  # GECSX format

    if vismin is not None and visib_field is None:
        warn(
            "Unable to filter data according to visibility. "
            + "Visibility field for RADAR"
            + "{:03d}".format(ind_radar_list[0] + 1)
            + " not available"
        )

    gate_coloc_rad1_dict = pyart.util.intersection(
        radar1,
        radar2,
        h_tol=h_tol,
        latlon_tol=latlon_tol,
        vol_d_tol=vol_d_tol,
        vismin=vismin,
        hmin=hmin,
        hmax=hmax,
        rmin=rmin,
        rmax=rmax,
        elmin=elmin,
        elmax=elmax,
        azmin=azrad1min,
        azmax=azrad1max,
        visib_field=visib_field,
        intersec_field=coloc_gates_field,
    )

    visib_field = None
    if "visibility" in radarnr_dict["RADAR" + "{:03d}".format(ind_radar_list[1] + 1)]:
        visib_field = "visibility"  # older IDL visibility codes
    elif (
        "visibility_polar"
        in radarnr_dict["RADAR" + "{:03d}".format(ind_radar_list[1] + 1)]
    ):
        visib_field = "visibility_polar"  # GECSX format

    if vismin is not None and visib_field is None:
        warn(
            "Unable to filter data according to visibility. "
            + "Visibility field for RADAR"
            + "{:03d}".format(ind_radar_list[1] + 1)
            + " not available"
        )

    gate_coloc_rad2_dict = pyart.util.intersection(
        radar2,
        radar1,
        h_tol=h_tol,
        latlon_tol=latlon_tol,
        vol_d_tol=vol_d_tol,
        vismin=vismin,
        hmin=hmin,
        hmax=hmax,
        rmin=rmin,
        rmax=rmax,
        elmin=elmin,
        elmax=elmax,
        azmin=azrad2min,
        azmax=azrad2max,
        visib_field=visib_field,
        intersec_field=coloc_gates_field,
    )

    new_rad1 = deepcopy(radar1)
    new_rad1.fields = dict()
    new_rad1.add_field("colocated_gates", gate_coloc_rad1_dict)

    new_rad2 = deepcopy(radar2)
    new_rad2.fields = dict()
    new_rad2.add_field("colocated_gates", gate_coloc_rad2_dict)

    coloc_rad1_dict, new_rad1.fields["colocated_gates"] = pyart.util.colocated_gates(
        new_rad1,
        new_rad2,
        h_tol=h_tol,
        latlon_tol=latlon_tol,
        coloc_gates_field=coloc_gates_field,
    )

    coloc_rad2_dict, new_rad2.fields["colocated_gates"] = pyart.util.colocated_gates(
        new_rad2,
        new_rad1,
        h_tol=h_tol,
        latlon_tol=latlon_tol,
        coloc_gates_field=coloc_gates_field,
    )

    # prepare output
    rad1_dict = {"coloc_dict": coloc_rad1_dict, "radar_out": new_rad1}

    rad2_dict = {"coloc_dict": coloc_rad2_dict, "radar_out": new_rad2}

    new_dataset = {
        "RADAR" + "{:03d}".format(ind_radar_list[0] + 1): rad1_dict,
        "RADAR" + "{:03d}".format(ind_radar_list[1] + 1): rad2_dict,
    }

    return new_dataset, ind_radar_list


def process_intercomp(procstatus, dscfg, radar_list=None):
    """
    intercomparison between two radars at co-located gates. The variables
    compared must be of the same type.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types (one for every radar).
            Any arbitrary datatype supported by pyrad and available
            in both radar is accepted.
        rays_are_indexed : bool. Dataset keyword
            If True it is considered that the rays are indexed and that the
            data can be selected simply looking at the ray number.
            Default false
        azi_tol : float. Dataset keyword
            azimuth tolerance between the two radars. Default 0.5 deg
        ele_tol : float. Dataset keyword
            elevation tolerance between the two radars. Default 0.5 deg
        rng_tol : float. Dataset keyword
            range tolerance between the two radars. Default 50 m
        coloc_gates_file : string. Dataset keyword
            name of the csv file containing the list of colocated data,
            as generated with the WRITE_COLOCATED_GATES product.
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing a dictionary with intercomparison data and the
        key "final" which contains a boolean that is true when all volumes
        have been processed
    ind_rad : int
        radar index

    """
    # check how many radars have to be compared and which datatype to use
    ind_radar_list = []
    field_name_list = []
    datatype_list = []
    for datatypedescr in dscfg["datatype"]:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        field_name = get_fieldname_pyart(datatype)
        ind_radar_list.append(int(radarnr[5:8]) - 1)
        datatype_list.append(datatype)
        field_name_list.append(field_name)
    ind_radar_list = np.array(ind_radar_list)
    datatype_list = np.array(datatype_list)
    field_name_list = np.array(field_name_list)
    if (ind_radar_list.size != 2) or (np.unique(ind_radar_list).size < 2):
        warn("Intercomparison requires data from two different radars")
        return None, None
    if np.unique(field_name_list).size > 1:
        warn("Intercomparison must be performed on the same variable")
        return None, None

    if procstatus == 0:
        coloc_gates_file = dscfg.get("coloc_gates_file", None)
        if coloc_gates_file:
            if not os.path.exists(coloc_gates_file):
                raise FileNotFoundError(
                    f"Unable to find colocated gates file {coloc_gates_file}"
                )

            coloc_gates_file
            (
                rad1_ray_ind,
                rad1_rng_ind,
                rad1_ele,
                rad1_azi,
                rad1_rng,
                rad2_ray_ind,
                rad2_rng_ind,
                rad2_ele,
                rad2_azi,
                rad2_rng,
            ) = read_colocated_gates(coloc_gates_file)

            if rad1_ele is None:
                raise ValueError(
                    "Unable to intercompare radars. " + "Invalid colocated gates file"
                )
            dscfg["global_data"] = {
                "rad1_ray_ind": rad1_ray_ind,
                "rad1_rng_ind": rad1_rng_ind,
                "rad1_ele": rad1_ele,
                "rad1_azi": rad1_azi,
                "rad1_rng": rad1_rng,
                "rad2_ray_ind": rad2_ray_ind,
                "rad2_rng_ind": rad2_rng_ind,
                "rad2_ele": rad2_ele,
                "rad2_azi": rad2_azi,
                "rad2_rng": rad2_rng,
            }

        return None, None

    if procstatus == 1:
        radar1 = radar_list[ind_radar_list[0]]
        radar2 = radar_list[ind_radar_list[1]]

        # rays are indexed to regular grid
        rays_are_indexed = dscfg.get("rays_are_indexed", False)

        if radar1 is None or radar2 is None:
            warn("Unable to inter-compare radars. Missing radar")
            return None, None

        if (field_name not in radar1.fields) or (field_name not in radar2.fields):
            warn(
                "Unable to get values of field "
                + field_name
                + " at colocated range bins. "
                + "Field missing in one of the radars"
            )
            return None, None

        if not dscfg["initialized"]:
            dscfg["global_data"].update({"timeinfo": dscfg["timeinfo"]})
            dscfg["global_data"].update(
                {"rad1_name": dscfg["RadarName"][ind_radar_list[0]]}
            )
            dscfg["global_data"].update(
                {"rad2_name": dscfg["RadarName"][ind_radar_list[1]]}
            )
            dscfg["initialized"] = 1

        rad1_field = radar1.fields[field_name]["data"]
        rad2_field = radar2.fields[field_name]["data"]

        intercomp_dict = {
            "rad1_time": [],
            "rad1_ray_ind": [],
            "rad1_rng_ind": [],
            "rad1_ele": [],
            "rad1_azi": [],
            "rad1_rng": [],
            "rad1_val": [],
            "rad2_time": [],
            "rad2_ray_ind": [],
            "rad2_rng_ind": [],
            "rad2_ele": [],
            "rad2_azi": [],
            "rad2_rng": [],
            "rad2_val": [],
            "rad1_name": dscfg["global_data"]["rad1_name"],
            "rad2_name": dscfg["global_data"]["rad2_name"],
        }

        # determine if radar data has to be averaged
        avg_rad1, avg_rad2, avg_rad_lim = get_range_bins_to_avg(
            radar1.range["data"], radar2.range["data"]
        )
        if not rays_are_indexed:
            azi_tol = dscfg.get("azi_tol", 0.5)
            ele_tol = dscfg.get("ele_tol", 0.5)
            rng_tol = dscfg.get("rng_tol", 50.0)
            (
                rad1_ray_ind,
                rad1_rng_ind,
                rad2_ray_ind,
                rad2_rng_ind,
            ) = find_colocated_indexes_fast(
                radar1,
                radar2,
                dscfg["global_data"]["rad1_ele"],
                dscfg["global_data"]["rad1_azi"],
                dscfg["global_data"]["rad1_rng"],
                dscfg["global_data"]["rad2_ele"],
                dscfg["global_data"]["rad2_azi"],
                dscfg["global_data"]["rad2_rng"],
                ele_tol=ele_tol,
                azi_tol=azi_tol,
                rng_tol=rng_tol,
            )
        else:
            rad1_ray_ind = deepcopy(dscfg["global_data"]["rad1_ray_ind"])
            rad1_rng_ind = deepcopy(dscfg["global_data"]["rad1_rng_ind"])
            rad2_ray_ind = deepcopy(dscfg["global_data"]["rad2_ray_ind"])
            rad2_rng_ind = deepcopy(dscfg["global_data"]["rad2_rng_ind"])

        # keep only indices of valid gates
        val1_vec = rad1_field[rad1_ray_ind, rad1_rng_ind]
        val2_vec = rad2_field[rad2_ray_ind, rad2_rng_ind]

        mask_val1 = np.ma.getmaskarray(val1_vec)
        mask_val2 = np.ma.getmaskarray(val2_vec)

        isvalid = np.logical_not(np.logical_or(mask_val1, mask_val2))

        rad1_ray_ind = rad1_ray_ind[isvalid]
        rad1_rng_ind = rad1_rng_ind[isvalid]
        rad2_ray_ind = rad2_ray_ind[isvalid]
        rad2_rng_ind = rad2_rng_ind[isvalid]

        # if averaging required loop over valid gates and average
        if avg_rad1:
            ngates_valid = len(rad1_ray_ind)
            val1_vec = np.ma.masked_all(ngates_valid, dtype=float)
            is_valid_avg = np.zeros(ngates_valid, dtype=bool)
            for i in range(ngates_valid):
                if rad1_rng_ind[i] + avg_rad_lim[1] >= radar1.ngates:
                    continue
                if rad1_rng_ind[i] + avg_rad_lim[0] < 0:
                    continue
                ind_rng = list(
                    range(
                        rad1_rng_ind[i] + avg_rad_lim[0],
                        rad1_rng_ind[i] + avg_rad_lim[1] + 1,
                    )
                )

                if np.any(np.ma.getmaskarray(rad1_field[rad1_ray_ind[i], ind_rng])):
                    continue

                val1_vec[i] = np.ma.asarray(
                    np.ma.mean(rad1_field[rad1_ray_ind[i], ind_rng])
                )

                is_valid_avg[i] = True

            rad1_ray_ind = rad1_ray_ind[is_valid_avg]
            rad1_rng_ind = rad1_rng_ind[is_valid_avg]
            rad2_ray_ind = rad2_ray_ind[is_valid_avg]
            rad2_rng_ind = rad2_rng_ind[is_valid_avg]

            val1_vec = val1_vec[is_valid_avg]
            val2_vec = rad2_field[rad2_ray_ind, rad2_rng_ind]

        elif avg_rad2:
            ngates_valid = len(rad2_ray_ind)
            val2_vec = np.ma.masked_all(ngates_valid, dtype=float)
            is_valid_avg = np.zeros(ngates_valid, dtype=bool)
            for i in range(ngates_valid):
                if rad2_rng_ind[i] + avg_rad_lim[1] >= radar2.ngates:
                    continue
                if rad2_rng_ind[i] + avg_rad_lim[0] < 0:
                    continue
                ind_rng = list(
                    range(
                        rad2_rng_ind[i] + avg_rad_lim[0],
                        rad2_rng_ind[i] + avg_rad_lim[1] + 1,
                    )
                )

                if np.any(np.ma.getmaskarray(rad2_field[rad2_ray_ind[i], ind_rng])):
                    continue

                val2_vec[i] = np.ma.asarray(
                    np.ma.mean(rad2_field[rad2_ray_ind[i], ind_rng])
                )

                is_valid_avg[i] = True

            rad1_ray_ind = rad1_ray_ind[is_valid_avg]
            rad1_rng_ind = rad1_rng_ind[is_valid_avg]
            rad2_ray_ind = rad2_ray_ind[is_valid_avg]
            rad2_rng_ind = rad2_rng_ind[is_valid_avg]

            val2_vec = val2_vec[is_valid_avg]
            val1_vec = rad1_field[rad1_ray_ind, rad1_rng_ind]
        else:
            val1_vec = val1_vec[isvalid]
            val2_vec = val2_vec[isvalid]

        intercomp_dict["rad1_time"] = cftodatetime(
            num2date(
                radar1.time["data"][rad1_ray_ind],
                radar1.time["units"],
                radar1.time["calendar"],
            )
        )
        intercomp_dict["rad1_ray_ind"] = rad1_ray_ind
        intercomp_dict["rad1_rng_ind"] = rad1_rng_ind
        intercomp_dict["rad1_ele"] = radar1.elevation["data"][rad1_ray_ind]
        intercomp_dict["rad1_azi"] = radar1.azimuth["data"][rad1_ray_ind]
        intercomp_dict["rad1_rng"] = radar1.range["data"][rad1_rng_ind]
        intercomp_dict["rad1_val"] = val1_vec

        intercomp_dict["rad2_time"] = cftodatetime(
            num2date(
                radar2.time["data"][rad2_ray_ind],
                radar2.time["units"],
                radar2.time["calendar"],
            )
        )
        intercomp_dict["rad2_ray_ind"] = rad2_ray_ind
        intercomp_dict["rad2_rng_ind"] = rad2_rng_ind
        intercomp_dict["rad2_ele"] = radar2.elevation["data"][rad2_ray_ind]
        intercomp_dict["rad2_azi"] = radar2.azimuth["data"][rad2_ray_ind]
        intercomp_dict["rad2_rng"] = radar2.range["data"][rad2_rng_ind]
        intercomp_dict["rad2_val"] = val2_vec

        new_dataset = {
            "intercomp_dict": intercomp_dict,
            "timeinfo": dscfg["global_data"]["timeinfo"],
            "final": False,
        }
        return new_dataset, None

    if procstatus == 2:
        tseries_prod = [
            prod
            for prod in dscfg["products"]
            if "WRITE_INTERCOMP" in dscfg["products"][prod]["type"]
        ][0]

        savedir = get_save_dir(
            dscfg["basepath"],
            dscfg["procname"],
            dscfg["dsname"],
            tseries_prod,
            timeinfo=dscfg["global_data"]["timeinfo"],
            create_dir=False,
        )

        fname = make_filename(
            "colocated_data",
            dscfg["type"],
            datatype,
            ["csv"],
            timeinfo=dscfg["global_data"]["timeinfo"],
            timeformat="%Y%m%d",
        )

        fname = savedir + fname[0]
        coloc_data = read_colocated_data(fname)
        intercomp_dict = {
            "rad1_name": dscfg["global_data"]["rad1_name"],
            "rad1_time": coloc_data[0],
            "rad1_ray_ind": coloc_data[1],
            "rad1_rng_ind": coloc_data[2],
            "rad1_ele": coloc_data[3],
            "rad1_azi": coloc_data[4],
            "rad1_rng": coloc_data[5],
            "rad1_val": coloc_data[6],
            "rad2_name": dscfg["global_data"]["rad2_name"],
            "rad2_time": coloc_data[7],
            "rad2_ray_ind": coloc_data[8],
            "rad2_rng_ind": coloc_data[9],
            "rad2_ele": coloc_data[10],
            "rad2_azi": coloc_data[11],
            "rad2_rng": coloc_data[12],
            "rad2_val": coloc_data[13],
        }

        new_dataset = {
            "intercomp_dict": intercomp_dict,
            "timeinfo": dscfg["global_data"]["timeinfo"],
            "final": True,
        }

        return new_dataset, None


def process_intercomp_with_QC(procstatus, dscfg, radar_list=None):
    """
    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types, must contain
            - one "main" datatype to be intercompared (ANY pyrad-supported datatype),
            - "PhiDP" or "PhiDPc"
            - "time_avg_flag"
            for the two radars
        colocgatespath : string.
            base path to the file containing the coordinates of the co-located
            gates
        coloc_data_dir : string. Dataset keyword
            name of the directory containing the csv file with colocated data
        coloc_radars_name : string. Dataset keyword
            string identifying the radar names
        rays_are_indexed : bool. Dataset keyword
            If True it is considered that the rays are indexed and that the
            data can be selected simply looking at the ray number.
            Default false
        azi_tol : float. Dataset keyword
            azimuth tolerance between the two radars. Default 0.5 deg
        ele_tol : float. Dataset keyword
            elevation tolerance between the two radars. Default 0.5 deg
        rng_tol : float. Dataset keyword
            range tolerance between the two radars. Default 50 m
        min_nb_samples_rad1: int. Dataset keyword
            minimum number of samples required for radar 1 to be included in the
            intercomparison. Default 1.
        min_nb_samples_rad2: int. Dataset keyword
            minimum number of samples required for radar 2 to be included in the
            intercomparison. Default 1.
        clt_max : float. Dataset keyword
            maximum fraction  of samples that can be clutter contaminated.
            Default 1. i.e. all
        phi_excess_max : float. Dataset keyword
            maximum fraction of samples that can have excess instantaneous
            PhiDP. Excess phidp is computed in the FLAG_TIME_AVG dataset, it corresponds to
            an exceedance of a threshold in PhiDP at the native time resolution of the data
            (no averaging). The excess phidp value corresponds to the parameter phidpmax used in the
            FLAG_TIME_AVG dataset used as input. Default 1. i.e. all
        non_rain_max : float. Dataset keyword
            maximum fraction of samples that can be no rain. Default 1. i.e. all
        phi_avg_max : float. Dataset keyword
            maximum average (hourly) PhiDP in degrees allowed. This conditions differs from
            phi_excess_max, because it concerns only average values and is defined in degrees
            and not in percentage of samples. Default is infinite, i.e. any

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing a dictionary with intercomparison data and the
        key "final" which contains a boolean that is true when all volumes
        have been processed
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        savedir = dscfg["colocgatespath"] + dscfg["coloc_radars_name"] + "/"

        prdtype = dscfg.get("prdtype", "info")

        fname = make_filename(
            prdtype,
            "COLOCATED_GATES",
            dscfg["coloc_radars_name"],
            ["csv"],
            timeinfo=None,
        )[0]

        (
            rad1_ray_ind,
            rad1_rng_ind,
            rad1_ele,
            rad1_azi,
            rad1_rng,
            rad2_ray_ind,
            rad2_rng_ind,
            rad2_ele,
            rad2_azi,
            rad2_rng,
        ) = read_colocated_gates(savedir + fname)

        if rad1_ele is None:
            raise ValueError(
                "Unable to intercompare radars. Missing colocated gates file"
            )

        dscfg["global_data"] = {
            "rad1_ray_ind": rad1_ray_ind,
            "rad1_rng_ind": rad1_rng_ind,
            "rad1_ele": rad1_ele,
            "rad1_azi": rad1_azi,
            "rad1_rng": rad1_rng,
            "rad2_ray_ind": rad2_ray_ind,
            "rad2_rng_ind": rad2_rng_ind,
            "rad2_ele": rad2_ele,
            "rad2_azi": rad2_azi,
            "rad2_rng": rad2_rng,
        }
        return None, None

    if procstatus == 1:
        # -----------------------------
        # Determine the two radars involved
        # -----------------------------
        ind_radar_set = set()
        for datatypedescr in dscfg["datatype"]:
            radarnr = datatypedescr.split(":")[0]
            ind_radar_set.add(int(radarnr[5:8]) - 1)
        ind_radar_list = list(ind_radar_set)
        if (len(ind_radar_list) != 2) or (radar_list is None) or (len(radar_list) < 2):
            warn("Intercomparison requires data from two different radars")
            return None, None

        radarnr_list = [
            "RADAR" + "{:03d}".format(ind_radar_list[0] + 1),
            "RADAR" + "{:03d}".format(ind_radar_list[1] + 1),
        ]

        # -----------------------------
        # Parse datatypes -> field names per radar
        # -----------------------------
        # main field = any datatype that is not PhiDP/PhiDPc/time_avg_flag
        rad1_main_field = rad1_phidp_field = rad1_flag_field = None
        rad2_main_field = rad2_phidp_field = rad2_flag_field = None

        main_type = None  # for procstatus==2 filename
        for datatypedescr in dscfg["datatype"]:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            fname = get_fieldname_pyart(datatype)

            is_phidp = datatype in ("PhiDP", "PhiDPc")
            is_flag = datatype == "time_avg_flag"
            is_main = (not is_phidp) and (not is_flag)

            if is_main and main_type is None:
                main_type = datatype

            if radarnr == radarnr_list[0]:
                if is_main:
                    rad1_main_field = fname
                elif is_phidp:
                    rad1_phidp_field = fname
                elif is_flag:
                    rad1_flag_field = fname
            elif radarnr == radarnr_list[1]:
                if is_main:
                    rad2_main_field = fname
                elif is_phidp:
                    rad2_phidp_field = fname
                elif is_flag:
                    rad2_flag_field = fname

        if main_type is None:
            warn(
                "No main datatype found (datatype must include one non-PhiDP/non-time_avg_flag entry)"
            )
            return None, None

        linear_avg = False
        if "reflectivity" in get_fieldname_pyart(main_type):
            linear_avg = True

        radar1 = radar_list[ind_radar_list[0]]
        radar2 = radar_list[ind_radar_list[1]]

        dscfg.setdefault("global_data", {})
        dscfg["global_data"].update({"timeinfo": dscfg["timeinfo"]})
        dscfg["global_data"].update({"main_type": main_type})

        if radar1 is None or radar2 is None:
            warn("Unable to inter-compare radars. Missing radar")
            return None, None

        if any(
            v is None
            for v in (
                rad1_main_field,
                rad1_phidp_field,
                rad1_flag_field,
                rad2_main_field,
                rad2_phidp_field,
                rad2_flag_field,
            )
        ):
            warn(
                "Unable to compare radar time avg fields. Missing datatype definitions"
            )
            return None, None

        if (
            (rad1_main_field not in radar1.fields)
            or (rad1_phidp_field not in radar1.fields)
            or (rad1_flag_field not in radar1.fields)
            or (rad2_main_field not in radar2.fields)
            or (rad2_phidp_field not in radar2.fields)
            or (rad2_flag_field not in radar2.fields)
        ):
            warn("Unable to compare radar time avg fields. Fields missing")
            return None, None

        if not dscfg.get("initialized", False):
            dscfg["global_data"].update(
                {"rad1_name": dscfg["RadarName"][ind_radar_list[0]]}
            )
            dscfg["global_data"].update(
                {"rad2_name": dscfg["RadarName"][ind_radar_list[1]]}
            )
            dscfg["initialized"] = True

        val1 = radar1.fields[rad1_main_field]["data"]
        val2 = radar2.fields[rad2_main_field]["data"]

        phidp1 = radar1.fields[rad1_phidp_field]["data"]
        phidp2 = radar2.fields[rad2_phidp_field]["data"]

        flag1 = radar1.fields[rad1_flag_field]["data"]
        flag2 = radar2.fields[rad2_flag_field]["data"]

        intercomp_dict = {
            "rad1_time": [],
            "rad1_ray_ind": [],
            "rad1_rng_ind": [],
            "rad1_ele": [],
            "rad1_azi": [],
            "rad1_rng": [],
            "rad1_val": [],
            "rad1_PhiDPavg": [],
            "rad1_Flagavg": [],
            "rad2_time": [],
            "rad2_ray_ind": [],
            "rad2_rng_ind": [],
            "rad2_ele": [],
            "rad2_azi": [],
            "rad2_rng": [],
            "rad2_val": [],
            "rad2_PhiDPavg": [],
            "rad2_Flagavg": [],
            "rad1_name": dscfg["global_data"]["rad1_name"],
            "rad2_name": dscfg["global_data"]["rad2_name"],
        }
        # determine if radar data has to be averaged
        avg_rad1, avg_rad2, avg_rad_lim1 = get_range_bins_to_avg(
            radar1.range["data"], radar2.range["data"]
        )
        avg_rad_lim2 = avg_rad_lim1

        # rays are indexed to regular grid
        rays_are_indexed = dscfg.get("rays_are_indexed", False)
        if not rays_are_indexed:
            azi_tol = dscfg.get("azi_tol", 0.5)
            ele_tol = dscfg.get("ele_tol", 0.5)
            rng_tol = dscfg.get("rng_tol", 50.0)

            (
                rad1_ray_ind,
                rad1_rng_ind,
                rad2_ray_ind,
                rad2_rng_ind,
            ) = find_colocated_indexes_fast(
                radar1,
                radar2,
                dscfg["global_data"]["rad1_ele"],
                dscfg["global_data"]["rad1_azi"],
                dscfg["global_data"]["rad1_rng"],
                dscfg["global_data"]["rad2_ele"],
                dscfg["global_data"]["rad2_azi"],
                dscfg["global_data"]["rad2_rng"],
                ele_tol=ele_tol,
                azi_tol=azi_tol,
                rng_tol=rng_tol,
            )
        else:
            rad1_ray_ind = deepcopy(dscfg["global_data"]["rad1_ray_ind"])
            rad1_rng_ind = deepcopy(dscfg["global_data"]["rad1_rng_ind"])
            rad2_ray_ind = deepcopy(dscfg["global_data"]["rad2_ray_ind"])
            rad2_rng_ind = deepcopy(dscfg["global_data"]["rad2_rng_ind"])

        # keep only indices and data of valid gates (val + phidp must be valid on both)
        val1_vec = val1[rad1_ray_ind, rad1_rng_ind]
        phidp1_vec = phidp1[rad1_ray_ind, rad1_rng_ind]
        flag1_vec = flag1[rad1_ray_ind, rad1_rng_ind]

        val2_vec = val2[rad2_ray_ind, rad2_rng_ind]
        phidp2_vec = phidp2[rad2_ray_ind, rad2_rng_ind]
        flag2_vec = flag2[rad2_ray_ind, rad2_rng_ind]

        mask_val1 = np.ma.getmaskarray(val1_vec)
        mask_phidp1 = np.ma.getmaskarray(phidp1_vec)
        mask_val2 = np.ma.getmaskarray(val2_vec)
        mask_phidp2 = np.ma.getmaskarray(phidp2_vec)

        isvalid = ~((mask_val1 | mask_val2) | (mask_phidp1 | mask_phidp2))

        rad1_ray_ind = rad1_ray_ind[isvalid]
        rad1_rng_ind = rad1_rng_ind[isvalid]
        rad2_ray_ind = rad2_ray_ind[isvalid]
        rad2_rng_ind = rad2_rng_ind[isvalid]

        # if averaging required loop over valid gates and average
        # only if all gates valid
        if avg_rad1:
            ngates_valid = len(rad1_ray_ind)
            val1_vec = np.ma.masked_all(ngates_valid, dtype=float)
            phidp1_vec = np.ma.masked_all(ngates_valid, dtype=float)
            flag1_vec = np.ma.masked_all(ngates_valid, dtype=int)
            is_valid_avg = np.zeros(ngates_valid, dtype=bool)

            for i in range(ngates_valid):
                if rad1_rng_ind[i] + avg_rad_lim1[1] >= radar1.ngates:
                    continue
                if rad1_rng_ind[i] + avg_rad_lim1[0] < 0:
                    continue

                ind_rng = list(
                    range(
                        rad1_rng_ind[i] + avg_rad_lim1[0],
                        rad1_rng_ind[i] + avg_rad_lim1[1] + 1,
                    )
                )

                if np.any(np.ma.getmaskarray(val1[rad1_ray_ind[i], ind_rng])):
                    continue
                if np.any(np.ma.getmaskarray(phidp1[rad1_ray_ind[i], ind_rng])):
                    continue

                if linear_avg:
                    val1_vec[i] = np.ma.asarray(
                        10
                        * np.log10(
                            np.ma.mean(10 ** (0.1 * val1[rad1_ray_ind[i], ind_rng]))
                        )
                    )
                else:
                    val1_vec[i] = np.ma.asarray(
                        np.ma.mean(val1[rad1_ray_ind[i], ind_rng])
                    )

                phidp1_vec[i] = np.ma.asarray(
                    np.ma.mean(phidp1[rad1_ray_ind[i], ind_rng])
                )

                rad1_flag = flag1[rad1_ray_ind[i], ind_rng]
                rad1_nsamples = rad1_flag % int(1e3)
                rad1_excess_phi = (rad1_flag // int(1e3)) % int(1e3)
                rad1_clt = (rad1_flag // int(1e6)) % int(1e3)
                rad1_non_rain = (rad1_flag // int(1e9)) % int(1e3)

                flag1_vec[i] = (
                    int(1e9) * np.max(rad1_non_rain)
                    + int(1e6) * np.max(rad1_clt)
                    + int(1e3) * np.max(rad1_excess_phi)
                    + np.max(rad1_nsamples)
                )
                is_valid_avg[i] = True

            rad1_ray_ind = rad1_ray_ind[is_valid_avg]
            rad1_rng_ind = rad1_rng_ind[is_valid_avg]
            rad2_ray_ind = rad2_ray_ind[is_valid_avg]
            rad2_rng_ind = rad2_rng_ind[is_valid_avg]

            val1_vec = val1_vec[is_valid_avg]
            phidp1_vec = phidp1_vec[is_valid_avg]
            flag1_vec = flag1_vec[is_valid_avg]

            val2_vec = val2[rad2_ray_ind, rad2_rng_ind]
            phidp2_vec = phidp2[rad2_ray_ind, rad2_rng_ind]
            flag2_vec = flag2[rad2_ray_ind, rad2_rng_ind]

        elif avg_rad2:
            ngates_valid = len(rad2_ray_ind)
            val2_vec = np.ma.masked_all(ngates_valid, dtype=float)
            phidp2_vec = np.ma.masked_all(ngates_valid, dtype=float)
            flag2_vec = np.ma.masked_all(ngates_valid, dtype=int)
            is_valid_avg = np.zeros(ngates_valid, dtype=bool)

            for i in range(ngates_valid):
                if rad2_rng_ind[i] + avg_rad_lim2[1] >= radar2.ngates:
                    continue
                if rad2_rng_ind[i] + avg_rad_lim2[0] < 0:
                    continue

                ind_rng = list(
                    range(
                        rad2_rng_ind[i] + avg_rad_lim2[0],
                        rad2_rng_ind[i] + avg_rad_lim2[1] + 1,
                    )
                )

                if np.any(np.ma.getmaskarray(val2[rad2_ray_ind[i], ind_rng])):
                    continue
                if np.any(np.ma.getmaskarray(phidp2[rad2_ray_ind[i], ind_rng])):
                    continue

                if linear_avg:
                    val2_vec[i] = np.ma.asarray(
                        10
                        * np.log10(
                            np.ma.mean(10 ** (0.1 * val2[rad2_ray_ind[i], ind_rng]))
                        )
                    )
                else:
                    val2_vec[i] = np.ma.asarray(
                        np.ma.mean(val2[rad2_ray_ind[i], ind_rng])
                    )

                phidp2_vec[i] = np.ma.asarray(
                    np.ma.mean(phidp2[rad2_ray_ind[i], ind_rng])
                )

                rad2_flag = flag2[rad2_ray_ind[i], ind_rng]
                rad2_nsamples = rad2_flag % int(1e3)
                rad2_excess_phi = (rad2_flag // int(1e3)) % int(1e3)
                rad2_clt = (rad2_flag // int(1e6)) % int(1e3)
                rad2_non_rain = (rad2_flag // int(1e9)) % int(1e3)

                flag2_vec[i] = (
                    int(1e9) * np.max(rad2_non_rain)
                    + int(1e6) * np.max(rad2_clt)
                    + int(1e3) * np.max(rad2_excess_phi)
                    + np.max(rad2_nsamples)
                )
                is_valid_avg[i] = True

            rad1_ray_ind = rad1_ray_ind[is_valid_avg]
            rad1_rng_ind = rad1_rng_ind[is_valid_avg]
            rad2_ray_ind = rad2_ray_ind[is_valid_avg]
            rad2_rng_ind = rad2_rng_ind[is_valid_avg]

            val2_vec = val2_vec[is_valid_avg]
            phidp2_vec = phidp2_vec[is_valid_avg]
            flag2_vec = flag2_vec[is_valid_avg]

            val1_vec = val1[rad1_ray_ind, rad1_rng_ind]
            phidp1_vec = phidp1[rad1_ray_ind, rad1_rng_ind]
            flag1_vec = flag1[rad1_ray_ind, rad1_rng_ind]
        else:
            val1_vec = val1_vec[isvalid]
            phidp1_vec = phidp1_vec[isvalid]
            flag1_vec = flag1_vec[isvalid]
            val2_vec = val2_vec[isvalid]
            phidp2_vec = phidp2_vec[isvalid]
            flag2_vec = flag2_vec[isvalid]

        # time is the same for all samples in this volume -> fill with timeinfo
        intercomp_dict["rad1_time"] = cftodatetime(
            num2date(
                radar1.time["data"][rad1_ray_ind],
                radar1.time["units"],
                radar1.time["calendar"],
            )
        )

        intercomp_dict["rad1_ray_ind"] = rad1_ray_ind
        intercomp_dict["rad1_rng_ind"] = rad1_rng_ind
        intercomp_dict["rad1_ele"] = radar1.elevation["data"][rad1_ray_ind]
        intercomp_dict["rad1_azi"] = radar1.azimuth["data"][rad1_ray_ind]
        intercomp_dict["rad1_rng"] = radar1.range["data"][rad1_rng_ind]
        intercomp_dict["rad1_val"] = val1_vec
        intercomp_dict["rad1_PhiDPavg"] = phidp1_vec
        intercomp_dict["rad1_Flagavg"] = flag1_vec

        intercomp_dict["rad2_time"] = cftodatetime(
            num2date(
                radar2.time["data"][rad2_ray_ind],
                radar2.time["units"],
                radar2.time["calendar"],
            )
        )
        intercomp_dict["rad2_ray_ind"] = rad2_ray_ind
        intercomp_dict["rad2_rng_ind"] = rad2_rng_ind
        intercomp_dict["rad2_ele"] = radar2.elevation["data"][rad2_ray_ind]
        intercomp_dict["rad2_azi"] = radar2.azimuth["data"][rad2_ray_ind]
        intercomp_dict["rad2_rng"] = radar2.range["data"][rad2_rng_ind]
        intercomp_dict["rad2_val"] = val2_vec
        intercomp_dict["rad2_PhiDPavg"] = phidp2_vec
        intercomp_dict["rad2_Flagavg"] = flag2_vec

        new_dataset = {
            "intercomp_dict": intercomp_dict,
            "timeinfo": dscfg["global_data"]["timeinfo"],
            "final": False,
        }
        return new_dataset, None

    if procstatus == 2:
        # Main datatype for filename (stored during procstatus==1; fallback to parsing)
        main_type = dscfg.get("global_data", {}).get("main_type", None)
        if main_type is None:
            for datatypedescr in dscfg["datatype"]:
                _, _, datatype, _, _ = get_datatype_fields(datatypedescr)
                if datatype not in ("PhiDP", "PhiDPc", "time_avg_flag"):
                    main_type = datatype
                    break

        if main_type is None:
            warn("No main datatype found")
            return None, None

        tseries_prod = [
            prod
            for prod in dscfg["products"]
            if "WRITE_INTERCOMP" in dscfg["products"][prod]["type"]
        ][0]

        savedir = get_save_dir(
            dscfg["basepath"],
            dscfg["procname"],
            dscfg["dsname"],
            tseries_prod,
            timeinfo=dscfg["global_data"]["timeinfo"],
            create_dir=False,
        )

        fname = make_filename(
            "colocated_data",
            dscfg["type"],
            main_type,
            ["csv"],
            timeinfo=dscfg["global_data"]["timeinfo"],
            timeformat="%Y%m%d",
        )[0]

        fname = savedir + fname

        (
            rad1_time,
            rad1_ray_ind,
            rad1_rng_ind,
            rad1_ele,
            rad1_azi,
            rad1_rng,
            rad1_val,
            rad1_phi,
            rad1_flag,
            rad2_time,
            rad2_ray_ind,
            rad2_rng_ind,
            rad2_ele,
            rad2_azi,
            rad2_rng,
            rad2_val,
            rad2_phi,
            rad2_flag,
        ) = read_colocated_data_with_QC(fname)

        if rad1_time is None:
            return None, None

        # Decode flags
        rad1_nsamples = rad1_flag % int(1e3)
        rad1_excess_phi = (rad1_flag // int(1e3)) % int(1e3)
        rad1_clt = (rad1_flag // int(1e6)) % int(1e3)
        rad1_non_rain = (rad1_flag // int(1e9)) % int(1e3)

        rad2_nsamples = rad2_flag % int(1e3)
        rad2_excess_phi = (rad2_flag // int(1e3)) % int(1e3)
        rad2_clt = (rad2_flag // int(1e6)) % int(1e3)
        rad2_non_rain = (rad2_flag // int(1e9)) % int(1e3)

        clt_max = dscfg.get("clt_max", 100)
        phi_excess_max = dscfg.get("phi_excess_max", 100)
        non_rain_max = dscfg.get("non_rain_max", 100)
        phi_avg_max = dscfg.get("phi_avg_max", np.inf)
        min_nb_samples_rad1 = dscfg.get("min_nb_samples_rad1", 1)
        min_nb_samples_rad2 = dscfg.get("min_nb_samples_rad2", 1)

        # Avoid division by zero
        rad1_ns = np.maximum(rad1_nsamples.astype(float), 1.0)
        rad2_ns = np.maximum(rad2_nsamples.astype(float), 1.0)

        # Filter out invalid data (percentages are in 0..100)
        ind_val = np.where(
            np.logical_and.reduce(
                (
                    (rad1_clt / rad1_ns) <= clt_max,
                    (rad2_clt / rad2_ns) <= clt_max,
                    (rad1_excess_phi / rad1_ns) <= phi_excess_max,
                    (rad2_excess_phi / rad2_ns) <= phi_excess_max,
                    (rad1_non_rain / rad1_ns) <= non_rain_max,
                    (rad2_non_rain / rad2_ns) <= non_rain_max,
                    (rad1_phi / rad1_ns) <= phi_avg_max,
                    (rad2_phi / rad2_ns) <= phi_avg_max,
                    (rad1_nsamples >= min_nb_samples_rad1),
                    (rad2_nsamples >= min_nb_samples_rad2),
                )
            )
        )[0]
        intercomp_dict = {
            "rad1_name": dscfg["global_data"]["rad1_name"],
            "rad1_time": rad1_time[ind_val],
            "rad1_ray_ind": rad1_ray_ind[ind_val],
            "rad1_rng_ind": rad1_rng_ind[ind_val],
            "rad1_ele": rad1_ele[ind_val],
            "rad1_azi": rad1_azi[ind_val],
            "rad1_rng": rad1_rng[ind_val],
            "rad1_val": rad1_val[ind_val],
            "rad2_name": dscfg["global_data"]["rad2_name"],
            "rad2_time": rad2_time[ind_val],
            "rad2_ray_ind": rad2_ray_ind[ind_val],
            "rad2_rng_ind": rad2_rng_ind[ind_val],
            "rad2_ele": rad2_ele[ind_val],
            "rad2_azi": rad2_azi[ind_val],
            "rad2_rng": rad2_rng[ind_val],
            "rad2_val": rad2_val[ind_val],
        }

        new_dataset = {
            "intercomp_dict": intercomp_dict,
            "timeinfo": dscfg["global_data"]["timeinfo"],
            "final": True,
        }
        return new_dataset, None

    return None, None


def process_fields_diff(procstatus, dscfg, radar_list=None):
    """
    Extracts the difference between two radars fields at a given time, assuming that
    their gates are all colocated. It does not do any time post-processing and
    only extracts the data at volume time.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types for each radar,
            Any datatype supported by pyrad is supported
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing a radar object containing the field differences
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # check how many radars are there
    radarnr_dict = dict()
    ind_radar_list = set()
    for datatypedescr in dscfg["datatype"]:
        radarnr = datatypedescr.split(":")[0]
        radarnr_dict.update({radarnr: []})
        ind_radar_list.add(int(radarnr[5:8]) - 1)

    ind_radar_list = list(ind_radar_list)

    if (len(radarnr_dict) != 2) or (len(radar_list) < 2):
        warn("Intercomparison requires data from two different radars")
        return None, None

    # create the list of data types for each radar
    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    field_name_1 = get_fieldname_pyart(datatype)

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][1])
    field_name_2 = get_fieldname_pyart(datatype)

    radar1 = radar_list[ind_radar_list[0]]
    radar2 = radar_list[ind_radar_list[1]]

    if radar1 is None or radar2 is None:
        warn("Unable to inter-compare radars. Missing radar")
        return None, None

    if (field_name_1 not in radar1.fields) or (field_name_2 not in radar2.fields):
        warn(
            "Unable to compare fields "
            + field_name_1
            + "and "
            + field_name_2
            + ". Field missing in one of the radars"
        )
        return None, None

    field_diff = pyart.config.get_metadata("fields_difference")
    field_diff["data"] = (
        radar1.fields[field_name_1]["data"] - radar2.fields[field_name_2]["data"]
    )
    field_diff["long_name"] = field_name_1 + " - " + field_name_2

    rad_diff = deepcopy(radar1)
    rad_diff.fields = dict()
    rad_diff.add_field("fields_difference", field_diff)

    new_dataset = {"radar_out": rad_diff}

    return new_dataset, None


def process_intercomp_fields(procstatus, dscfg, radar_list=None):
    """
    Intercomparison between two radars fields at a given time.
    This functions extracts the data from both radars, assuming that
    their gates are all colocated. It does not do any time post-processing and
    only extracts the data at volume time.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types for each radar,
            Any datatype supported by pyrad and available in both radars is supported
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing a dictionary with intercomparison data
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # check how many radars are there
    radarnr_dict = dict()
    ind_radar_list = set()
    for datatypedescr in dscfg["datatype"]:
        radarnr = datatypedescr.split(":")[0]
        radarnr_dict.update({radarnr: []})
        ind_radar_list.add(int(radarnr[5:8]) - 1)

    ind_radar_list = list(ind_radar_list)

    if (len(radarnr_dict) != 2) or (len(radar_list) < 2):
        warn("Intercomparison requires data from two different radars")
        return None, None

    # create the list of data types for each radar
    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][0])
    field_name_1 = get_fieldname_pyart(datatype)

    radarnr, _, datatype, _, _ = get_datatype_fields(dscfg["datatype"][1])
    field_name_2 = get_fieldname_pyart(datatype)

    radar1 = radar_list[ind_radar_list[0]]
    radar2 = radar_list[ind_radar_list[1]]

    if radar1 is None or radar2 is None:
        warn("Unable to inter-compare radars. Missing radar")
        return None, None

    if (field_name_1 not in radar1.fields) or (field_name_2 not in radar2.fields):
        warn(
            "Unable to compare fields "
            + field_name_1
            + " and "
            + field_name_2
            + ". Field missing in one of the radars"
        )
        return None, None

    data1 = deepcopy(radar1.fields[field_name_1]["data"])
    data2 = deepcopy(radar2.fields[field_name_2]["data"])
    mask1 = np.ma.getmaskarray(data1)
    mask2 = np.ma.getmaskarray(data2)

    data1[mask2] = np.ma.masked
    data2[mask1] = np.ma.masked

    intercomp_dict = {
        "rad1_name": dscfg["RadarName"][ind_radar_list[0]],
        "rad1_val": data1.compressed(),
        "rad2_name": dscfg["RadarName"][ind_radar_list[1]],
        "rad2_val": data2.compressed(),
    }

    # Add times
    intercomp_dict["rad1_time"] = cftodatetime(
        num2date(
            radar1.time["data"],
            radar1.time["units"],
            radar1.time["calendar"],
        )
    )
    intercomp_dict["rad2_time"] = cftodatetime(
        num2date(
            radar1.time["data"],
            radar1.time["units"],
            radar1.time["calendar"],
        )
    )

    new_dataset = {
        "intercomp_dict": intercomp_dict,
        "timeinfo": dscfg["timeinfo"],
        "final": False,
    }

    return new_dataset, None
