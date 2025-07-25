"""
==================================
Input and output (:mod:`pyrad.io`)
==================================

.. currentmodule:: pyrad.io

Functions to read and write data and configuration files.

Reading configuration files
===========================

.. autosummary::
    :toctree: generated/

    read_config

Reading radar data
==================

.. autosummary::
    :toctree: generated/

    get_data

Reading icon data
==================

.. autosummary::
    :toctree: generated/

    icon2radar_data
    icon2radar_coord
    hzt2radar_data
    hzt2radar_coord
    get_icon_fields
    get_iso0_field
    read_icon_data
    read_icon_coord
    read_hzt_data
    read_iso0_mf_data
    read_iso0_grib_data
    iso2radar_data
    grib2radar_data
    get_iso0_ref

Reading DEM data
==================

.. autosummary::
    :toctree: generated/

    dem2radar_data
    read_idrisi_data
    read_idrisi_metadata

Reading other data
==================

.. autosummary::
    :toctree: generated/

    read_centroids_npz
    read_centroids
    read_proc_periods
    read_last_state
    read_status
    read_rad4alp_icon
    read_rad4alp_vis
    read_excess_gates
    read_colocated_gates
    read_colocated_data
    read_timeseries
    read_ts_cum
    read_monitoring_ts
    read_intercomp_scores_ts
    get_sensor_data
    read_smn
    read_smn2
    read_coord_sensors
    read_disdro_scattering
    read_disdro_parsivel
    read_knmi
    read_sun_hits
    read_sun_hits_multiple_days
    read_sun_retrieval
    read_solar_flux
    read_selfconsistency
    read_antenna_pattern
    read_meteorage
    read_lightning
    read_lightning_traj
    read_lightning_all
    read_trt_scores
    read_trt_data
    read_trt_traj_data
    read_trt_thundertracking_traj_data
    read_trt_cell_lightning
    read_trt_info_all
    read_trt_info_all2
    read_trt_info
    read_trt_info2
    read_thundertracking_info
    read_rhi_profile
    read_vpr_theo_parameters
    read_histogram
    read_quantiles
    read_profile_ts
    read_histogram_ts
    read_quantiles_ts
    read_ml_ts
    read_windmills_data
    read_radiosounding_wyoming
    read_radiosounding_igra
    read_fzl_igra
    read_meteoswiss_xml_vad

Writing data
==================

.. autosummary::
    :toctree: generated/

    write_vol_csv
    write_vol_kml
    write_centroids
    write_proc_periods
    write_ts_lightning
    send_msg
    write_alarm_msg
    write_last_state
    write_smn
    write_trt_info
    write_trt_thundertracking_data
    write_trt_cell_data
    write_trt_cell_scores
    write_trt_cell_lightning
    write_trt_rpc
    write_vpr_theo_params
    write_vpr_info
    write_rhi_profile
    write_field_coverage
    write_cdf
    write_histogram
    write_quantiles
    write_multiple_points
    write_multiple_points_grid
    write_ts_polar_data
    write_ts_grid_data
    write_ts_cum
    write_ts_stats
    write_monitoring_ts
    write_excess_gates
    write_intercomp_scores_ts
    write_colocated_gates
    write_colocated_data
    write_colocated_data_time_avg
    write_sun_hits
    write_sun_retrieval
    write_fixed_angle


Auxiliary functions
===================

.. autosummary::
    :toctree: generated/

    get_rad4alp_prod_fname
    get_datatype_skyecho
    map_hydro
    map_Doppler
    get_save_dir
    make_filename
    generate_field_name_str
    get_fieldname_pyart
    get_fieldname_icon
    get_field_unit
    get_file_list
    get_file_list_s3
    get_scan_files_to_merge
    get_scan_files_to_merge_s3
    get_rad4alp_dir
    get_rad4alp_grid_dir
    get_trtfile_list
    get_new_rainbow_file_name
    get_datatype_fields
    get_dataset_fields
    get_datetime
    find_raw_icon_file
    find_hzt_file
    find_iso0_file
    find_iso0_grib_file
    find_date_in_file_name
    _get_datetime

Trajectory
==========

.. autosummary::
    :toctree: generated/

    Trajectory

TimeSeries
==========

.. autosummary::
    :toctree: generated/

    TimeSeries

"""

from .config import read_config  # noqa
from .default_config import DEFAULT_CONFIG  # noqa
from .read_data_radar import get_data, add_field, interpol_field  # noqa

from .read_data_icon import read_icon_data, read_icon_coord  # noqa
from .read_data_icon import icon2radar_data, icon2radar_coord  # noqa
from .read_data_icon import get_icon_fields  # noqa

from .read_data_dem import read_idrisi_data, read_idrisi_metadata  # noqa
from .read_data_dem import dem2radar_data  # noqa

from .read_data_hzt import read_hzt_data, hzt2radar_data, hzt2radar_coord  # noqa
from .read_data_hzt import get_iso0_field  # noqa

from .read_data_iso0_mf import read_iso0_mf_data, iso2radar_data, get_iso0_ref  # noqa
from .read_data_iso0_mf import read_iso0_grib_data, grib2radar_data  # noqa

from .read_data_other import read_status, read_rad4alp_icon, read_rad4alp_vis  # noqa
from .read_data_other import read_timeseries, read_monitoring_ts, read_ts_cum  # noqa
from .read_data_other import read_intercomp_scores_ts, read_quantiles  # noqa
from .read_data_other import read_selfconsistency, read_colocated_gates  # noqa
from .read_data_other import read_colocated_data, read_antenna_pattern  # noqa
from .read_data_other import read_colocated_data_time_avg  # noqa
from .read_data_other import read_last_state, read_rhi_profile, read_centroids  # noqa
from .read_data_other import read_excess_gates, read_histogram  # noqa
from .read_data_other import read_profile_ts, read_histogram_ts  # noqa
from .read_data_other import read_quantiles_ts, read_ml_ts, read_proc_periods  # noqa
from .read_data_other import read_centroids_npz, read_mf_vis  # noqa
from .read_data_other import read_vpr_theo_parameters  # noqa
from .read_data_other import read_mch_xml_vad  # noqa

from .read_data_sensor import read_lightning, read_lightning_traj  # noqa
from .read_data_sensor import get_sensor_data, read_smn, read_smn2  # noqa
from .read_data_sensor import read_disdro_scattering, read_trt_data  # noqa
from .read_data_sensor import read_trt_traj_data, read_lightning_all  # noqa
from .read_data_sensor import read_trt_scores, read_trt_cell_lightning  # noqa
from .read_data_sensor import read_meteorage, read_trt_info_all, read_trt_info  # noqa
from .read_data_sensor import read_thundertracking_info, read_windmills_data  # noqa
from .read_data_sensor import read_trt_info2, read_trt_info_all2  # noqa
from .read_data_sensor import read_trt_thundertracking_traj_data  # noqa
from .read_data_sensor import read_coord_sensors, read_disdro_parsivel  # noqa
from .read_data_sensor import read_knmi  # noqa
from .read_data_sensor import (  # noqa
    read_radiosounding_wyoming,
    read_radiosounding_igra,
)
from .read_data_sensor import read_fzl_igra  # noqa

from .read_data_sun import read_sun_hits_multiple_days, read_sun_hits  # noqa
from .read_data_sun import read_sun_retrieval, read_solar_flux  # noqa

from .write_data import write_smn, write_ts_polar_data, write_ts_cum  # noqa
from .write_data import write_monitoring_ts, write_intercomp_scores_ts  # noqa
from .write_data import write_sun_hits, write_sun_retrieval  # noqa
from .write_data import write_colocated_gates, write_colocated_data  # noqa
from .write_data import write_colocated_data_time_avg, write_cdf  # noqa
from .write_data import write_rhi_profile, write_field_coverage  # noqa
from .write_data import write_last_state, write_alarm_msg, send_msg  # noqa
from .write_data import write_excess_gates, write_trt_cell_data  # noqa
from .write_data import write_histogram, write_quantiles, write_ts_lightning  # noqa
from .write_data import write_trt_cell_scores, write_trt_cell_lightning  # noqa
from .write_data import write_trt_info, write_fixed_angle, write_proc_periods  # noqa
from .write_data import write_trt_thundertracking_data, write_ts_grid_data  # noqa
from .write_data import write_trt_rpc, write_ts_stats, write_centroids  # noqa
from .write_data import write_multiple_points, write_vpr_theo_params  # noqa
from .write_data import write_multiple_points_grid, write_vpr_info  # noqa
from .write_data import write_vol_kml, write_vol_csv  # noqa

from .io_aux import get_save_dir, make_filename, get_new_rainbow_file_name  # noqa
from .io_aux import get_datetime, get_dataset_fields, map_hydro, map_Doppler  # noqa
from .io_aux import get_file_list, get_trtfile_list, get_datatype_fields  # noqa
from .io_aux import get_fieldname_pyart, get_field_unit, get_fieldname_icon  # noqa
from .io_aux import generate_field_name_str, find_raw_icon_file  # noqa
from .io_aux import find_hzt_file, _get_datetime, get_rad4alp_prod_fname  # noqa
from .io_aux import get_rad4alp_dir, get_rad4alp_grid_dir, find_iso0_file  # noqa
from .io_aux import find_iso0_grib_file, find_date_in_file_name  # noqa
from .io_aux import get_datatype_skyecho, get_file_list_s3  # noqa
from .io_aux import get_scan_files_to_merge, get_scan_files_to_merge_s3  # noqa

from .trajectory import Trajectory  # noqa

from .timeseries import TimeSeries  # noqa

__all__ = [s for s in dir() if not s.startswith("_")]
