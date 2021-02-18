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

Reading cosmo data
==================

.. autosummary::
    :toctree: generated/

    cosmo2radar_data
    cosmo2radar_coord
    hzt2radar_data
    hzt2radar_coord
    get_cosmo_fields
    get_iso0_field
    read_cosmo_data
    read_cosmo_coord
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
    dem2radar_coord
    read_idrisi_data
    read_idrisi_metadata

Reading other data
==================

.. autosummary::
    :toctree: generated/

    read_proc_periods
    read_last_state
    read_status
    read_rad4alp_cosmo
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
    read_disdro_scattering
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
    read_histogram
    read_quantiles
    read_profile_ts
    read_histogram_ts
    read_quantiles_ts
    read_ml_ts
    read_windmills_data

Writing data
==================

.. autosummary::
    :toctree: generated/

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
    write_rhi_profile
    write_field_coverage
    write_cdf
    write_histogram
    write_quantiles
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
    map_hydro
    map_Doppler
    get_save_dir
    make_filename
    generate_field_name_str
    get_fieldname_pyart
    get_fieldname_cosmo
    get_field_unit
    get_file_list
    get_rad4alp_dir
    get_rad4alp_grid_dir
    get_trtfile_list
    get_new_rainbow_file_name
    get_datatype_fields
    get_dataset_fields
    get_datetime
    find_raw_cosmo_file
    find_hzt_file
    find_iso0_file
    find_iso0_grib_file
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

from .config import read_config

from .read_data_radar import get_data, add_field, interpol_field

from .read_data_cosmo import read_cosmo_data, read_cosmo_coord
from .read_data_cosmo import cosmo2radar_data, cosmo2radar_coord
from .read_data_cosmo import get_cosmo_fields

from .read_data_dem import read_idrisi_data, read_idrisi_metadata

from .read_data_hzt import read_hzt_data, hzt2radar_data, hzt2radar_coord
from .read_data_hzt import get_iso0_field

from .read_data_iso0_mf import read_iso0_mf_data, iso2radar_data, get_iso0_ref
from .read_data_iso0_mf import read_iso0_grib_data, grib2radar_data

from .read_data_other import read_status, read_rad4alp_cosmo, read_rad4alp_vis
from .read_data_other import read_timeseries, read_monitoring_ts, read_ts_cum
from .read_data_other import read_intercomp_scores_ts, read_quantiles
from .read_data_other import read_selfconsistency, read_colocated_gates
from .read_data_other import read_colocated_data, read_antenna_pattern
from .read_data_other import read_last_state, read_rhi_profile
from .read_data_other import read_excess_gates, read_histogram
from .read_data_other import read_profile_ts, read_histogram_ts
from .read_data_other import read_quantiles_ts, read_ml_ts, read_proc_periods

from .read_data_sensor import read_lightning, read_lightning_traj
from .read_data_sensor import get_sensor_data, read_smn, read_smn2
from .read_data_sensor import read_disdro_scattering, read_trt_data
from .read_data_sensor import read_trt_traj_data, read_lightning_all
from .read_data_sensor import read_trt_scores, read_trt_cell_lightning
from .read_data_sensor import read_meteorage, read_trt_info_all, read_trt_info
from .read_data_sensor import read_thundertracking_info, read_windmills_data
from .read_data_sensor import read_trt_info2, read_trt_info_all2
from .read_data_sensor import read_trt_thundertracking_traj_data

from .read_data_sun import read_sun_hits_multiple_days, read_sun_hits
from .read_data_sun import read_sun_retrieval, read_solar_flux

from .write_data import write_smn, write_ts_polar_data, write_ts_cum
from .write_data import write_monitoring_ts, write_intercomp_scores_ts
from .write_data import write_sun_hits, write_sun_retrieval
from .write_data import write_colocated_gates, write_colocated_data
from .write_data import write_colocated_data_time_avg, write_cdf
from .write_data import write_rhi_profile, write_field_coverage
from .write_data import write_last_state, write_alarm_msg, send_msg
from .write_data import write_excess_gates, write_trt_cell_data
from .write_data import write_histogram, write_quantiles, write_ts_lightning
from .write_data import write_trt_cell_scores, write_trt_cell_lightning
from .write_data import write_trt_info, write_fixed_angle, write_proc_periods
from .write_data import write_trt_thundertracking_data, write_ts_grid_data
from .write_data import write_trt_rpc, write_ts_stats

from .io_aux import get_save_dir, make_filename, get_new_rainbow_file_name
from .io_aux import get_datetime, get_dataset_fields, map_hydro, map_Doppler
from .io_aux import get_file_list, get_trtfile_list, get_datatype_fields
from .io_aux import get_fieldname_pyart, get_field_unit, get_fieldname_cosmo
from .io_aux import generate_field_name_str, find_raw_cosmo_file
from .io_aux import find_hzt_file, _get_datetime, get_rad4alp_prod_fname
from .io_aux import get_rad4alp_dir, get_rad4alp_grid_dir, find_iso0_file
from .io_aux import find_iso0_grib_file

from .trajectory import Trajectory

from .timeseries import TimeSeries

__all__ = [s for s in dir() if not s.startswith('_')]
