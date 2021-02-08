"""
pyrad.prod.process_vol_products
===============================

Functions for obtaining Pyrad products from a radar volume dataset

.. autosummary::
    :toctree: generated/

    generate_vol_products

"""

from copy import deepcopy
from warnings import warn

import numpy as np
import pyart
from netCDF4 import num2date

from ..io.io_aux import get_save_dir, make_filename, get_fieldname_pyart
from ..io.io_aux import generate_field_name_str

from ..io.write_data import write_cdf, write_rhi_profile, write_field_coverage
from ..io.write_data import write_last_state, write_histogram, write_quantiles
from ..io.write_data import write_fixed_angle, write_monitoring_ts
from ..io.write_data import write_alarm_msg, write_timeseries_point, send_msg

from ..io.read_data_other import read_monitoring_ts

from ..graph.plots_vol import plot_ppi, plot_ppi_map, plot_rhi, plot_cappi
from ..graph.plots_vol import plot_bscope, plot_rhi_profile, plot_along_coord
from ..graph.plots_vol import plot_field_coverage, plot_time_range
from ..graph.plots_vol import plot_rhi_contour, plot_ppi_contour
from ..graph.plots_vol import plot_fixed_rng, plot_fixed_rng_span
from ..graph.plots_vol import plot_roi_contour, plot_ray
from ..graph.plots import plot_quantiles, plot_histogram
from ..graph.plots import plot_selfconsistency_instrument
from ..graph.plots import plot_selfconsistency_instrument2
from ..graph.plots_aux import get_colobar_label, get_field_name
from ..graph.plots_timeseries import plot_monitoring_ts

from ..util.radar_utils import get_ROI, compute_profile_stats
from ..util.radar_utils import compute_histogram, compute_quantiles
from ..util.radar_utils import get_data_along_rng, get_data_along_azi
from ..util.radar_utils import get_data_along_ele
from ..util.stat_utils import quantiles_weighted


def generate_vol_products(dataset, prdcfg):
    """
    Generates radar volume products. Accepted product types:
        'CDF': plots and writes the cumulative density function of data
            User defined parameters:
                quantiles: list of floats
                    The quantiles to compute in percent. Default None
                sector: dict
                    dictionary defining the sector where to compute the CDF.
                    Default is None and the CDF is computed over all the data
                    May contain:
                        rmin, rmax: float
                            min and max range [m]
                        azmin, azmax: float
                            min and max azimuth angle [deg]
                        elmin, elmax: float
                            min and max elevation angle [deg]
                        hmin, hmax: float
                            min and max altitude [m MSL]
                vismin: float
                    The minimum visibility to use the data. Default None
                absolute: Bool
                    If true the absolute values of the data will be used.
                    Default False
                use_nans: Bool
                    If true NaN values will be used. Default False
                nan_value: Bool
                    The value by which the NaNs are substituted if NaN values
                    are to be used in the computation
                filterclt: Bool
                    If True the gates containing clutter are filtered
                filterprec: list of ints
                    The hydrometeor types that are filtered from the analysis.
                    Default empty list.
        'BSCOPE_IMAGE': Creates a B-scope image (azimuth, range)
            User defined parameters:
                anglenr : int
                    The elevation angle number to use
                ray_dim : str
                    the ray dimension. Can be 'ang' or 'time'. Default 'ang'
                xaxis_rng : bool
                    if True the range will be in the x-axis. Otherwise it will
                    be in the y-axis. Default True
                vmin, vmax: float or None
                    The minimum and maximum values of the color scale. If None
                    the scale is going to be set according to the Py-ART
                    config file
        'CAPPI_IMAGE': Creates a CAPPI image
            User defined parameters:
                altitude: flt
                    CAPPI altitude [m MSL]
                wfunc: str
                    The function used to produce the CAPPI as defined in
                    pyart.map.grid_from_radars. Default 'NEAREST'
                cappi_res: float
                    The CAPPI resolution [m]. Default 500.
        'FIELD_COVERAGE': Gets the field coverage over a certain sector
            User defined parameters:
                threshold: float or None
                    Minimum value to consider the data valid. Default None
                nvalid_min: float
                    Minimum number of valid gates in the ray to consider it
                    valid. Default 5
                ele_res, azi_res: float
                    Elevation and azimuth resolution of the sectors [deg].
                    Default 1. and 2.
                ele_min, ele_max: float
                    Min and max elevation angle defining the sector [deg].
                    Default 0. and 30.
                ele_step: float
                    Elevation step [deg]. Default 5.
                ele_sect_start, ele_sect_stop: float or None
                    start and stop angles of the sector coverage. Default None
                quantiles: list of floats
                    The quantiles to compute in the sector. Default 10. to 90.
                    by steps of 10.
                AngTol: float
                    The tolerance in elevation angle when putting the data in
                    a fixed grid
        'FIXED_RNG_IMAGE': Plots a fixed range image
            User defined parameters:
                AngTol : float
                    The tolerance between the nominal angles and the actual
                    radar angles. Default 1.
                ele_res, azi_res: float or None
                    The resolution of the fixed grid [deg]. If None it will be
                    obtained from the separation between angles
                vmin, vmax : float or None
                    Min and Max values of the color scale. If None the values
                    are taken from the Py-ART config file
        'FIXED_RNG_SPAN_IMAGE': Plots a user-defined statistic over a fixed
            range image
            User defined parameters:
                AngTol : float
                    The tolerance between the nominal angles and the actual
                    radar angles. Default 1.
                ele_res, azi_res: float or None
                    The resolution of the fixed grid [deg]. If None it will be
                    obtained from the separation between angles
                stat : str
                    The statistic to compute. Can be 'min', 'max', 'mean',
                    'mode'. Default 'max'
        'HISTOGRAM': Computes a histogram of the radar volum data
            User defined parameters:
                step: float or None
                    the data quantization step. If none it will be obtained
                    from the Py-ART configuration file
                write_data: Bool
                    If true the histogram data is written in a csv file
        'PLOT_ALONG_COORD': Plots the radar volume data along a particular
            coordinate
            User defined parameters:
                colors: list of str or None
                    The colors of each ploted line
                mode: str
                    Ploting mode. Can be 'ALONG_RNG', 'ALONG_AZI' or
                    'ALONG_ELE'
                value_start, value_stop: float
                    The starting and ending points of the data to plot.
                    According to the mode it may refer to the range, azimuth
                    or elevation. If not specified the minimum and maximum
                    possible values are used
                fix_elevations, fix_azimuths, fix_ranges: list of floats
                    The elevations, azimuths or ranges to plot for each mode.
                    'ALONG_RNG' would use fix_elevations and fix_azimuths
                    'ALONG_AZI' fix_ranges and fix_elevations
                    'ALONG_ELE' fix_ranges and fix_azimuths
                AngTol: float
                    The tolerance to match the radar angle to the fixed angles
                    Default 1.
                RngTol: float
                    The tolerance to match the radar range to the fixed ranges
                    Default 50.
        'PLOT_TXH': Plots the transmitted signal power (H) for a standard sunscan.
        'PPI_CONTOUR': Plots a PPI countour plot
            User defined parameters:
                contour_values: list of floats or None
                    The list of contour values to plot. If None the contour
                    values are going to be obtained from the Py-ART config
                    file either with the dictionary key 'contour_values' or
                    from the minimum and maximum values of the field with an
                    assumed division of 10 levels.
                anglenr: float
                    The elevation angle number
        'PPI_CONTOUR_OVERPLOT': Plots a PPI of a field with another field
            overplotted as a contour plot.
            User defined parameters:
                contour_values: list of floats or None
                    The list of contour values to plot. If None the contour
                    values are going to be obtained from the Py-ART config
                    file either with the dictionary key 'contour_values' or
                    from the minimum and maximum values of the field with an
                    assumed division of 10 levels.
                anglenr: float
                    The elevation angle number
        'PPI_IMAGE': Plots a PPI image. It can also plot the histogram and the
            quantiles of the data in the PPI.
            User defined parameters:
                anglenr: float
                    The elevation angle number
                plot_type: str
                    The type of plot to perform. Can be 'PPI', 'QUANTILES' or
                    'HISTOGRAM'
                step: float or None
                    If the plot type is 'HISTOGRAM', the width of the
                    histogram bin. If None it will be obtained from the Py-ART
                    config file
                quantiles: list of float or None
                    If the plot type is 'QUANTILES', the list of quantiles to
                    compute. If None a default list of quantiles will be
                    computed
                vmin, vmax: float or None
                    The minimum and maximum values of the color scale. If None
                    the scale is going to be set according to the Py-ART
                    config file
        'PPI_MAP': Plots a PPI image over a map. The map resolution and the
            type of maps used are defined in the variables 'mapres' and 'maps'
            in 'ppiMapImageConfig' in the loc config file.
            User defined parameters:
                anglenr: float
                    The elevation angle number
        'PPIMAP_ROI_OVERPLOT': Over plots a polygon delimiting a region of
            interest on a PPI map. The map resolution and the
            type of maps used are defined in the variables 'mapres' and 'maps'
            in 'ppiMapImageConfig' in the loc config file.
            User defined parameters:
                anglenr: float
                    The elevation angle number
        'PROFILE_STATS': Computes and plots a vertical profile statistics.
            The statistics are saved in a csv file
            User defined parameters:
                heightResolution: float
                    The height resolution of the profile [m]. Default 100.
                heightMin, heightMax: float or None
                    The minimum and maximum altitude of the profile [m MSL].
                    If None the values will be obtained from the minimum and
                    maximum gate altitude.
                quantity: str
                    The type of statistics to plot. Can be 'quantiles',
                    'mode', 'reqgression_mean' or 'mean'.
                quantiles: list of floats
                    If quantity type is 'quantiles' the list of quantiles to
                    compute. Default 25., 50., 75.
                nvalid_min: int
                    The minimum number of valid points to consider the
                    statistic valid. Default 4
                make_linear: Bool
                    If true the data is converted from log to linear before
                    computing the stats
                include_nans: Bool
                    If true NaN values are included in the statistics
                fixed_span: Bool
                    If true the profile plot has a fix X-axis
                vmin, vmax: float or None
                    If fixed_span is set, the minimum and maximum values of
                    the X-axis. If None, they are obtained from the Py-ART
                    config file
        'PSEUDOPPI_CONTOUR': Plots a pseudo-PPI countour plot
            User defined parameters:
                contour_values: list of floats or None
                    The list of contour values to plot. If None the contour
                    values are going to be obtained from the Py-ART config
                    file either with the dictionary key 'contour_values' or
                    from the minimum and maximum values of the field with an
                    assumed division of 10 levels.
                angle: float
                    The elevation angle at which compute the PPI
                EleTol: float
                    The tolerance between the actual radar elevation angle and
                    the nominal pseudo-PPI elevation angle.
        'PSEUDOPPI_CONTOUR_OVERPLOT': Plots a pseudo-PPI of a field with
            another field over-plotted as a contour plot
            User defined parameters:
                contour_values: list of floats or None
                    The list of contour values to plot. If None the contour
                    values are going to be obtained from the Py-ART config
                    file either with the dictionary key 'contour_values' or
                    from the minimum and maximum values of the field with an
                    assumed division of 10 levels.
                angle: float
                    The elevation angle at which compute the PPI
                EleTol: float
                    The tolerance between the actual radar elevation angle and
                    the nominal pseudo-PPI elevation angle.
        'PSEUDOPPI_IMAGE': Plots a pseudo-PPI image. It can also plot the
            histogram and the quantiles of the data in the pseudo-PPI.
            User defined parameters:
                angle: float
                    The elevation angle of the pseudo-PPI
                EleTol: float
                    The tolerance between the actual radar elevation angle and
                    the nominal pseudo-PPI elevation angle.
                plot_type: str
                    The type of plot to perform. Can be 'PPI', 'QUANTILES' or
                    'HISTOGRAM'
                step: float or None
                    If the plot type is 'HISTOGRAM', the width of the
                    histogram bin. If None it will be obtained from the Py-ART
                    config file
                quantiles: list of float or None
                    If the plot type is 'QUANTILES', the list of quantiles to
                    compute. If None a default list of quantiles will be
                    computed
                vmin, vmax : float or None
                    Min and Max values of the color scale. If None the values
                    are taken from the Py-ART config file
        'PSEUDOPPI_MAP': Plots a pseudo-PPI image over a map. The map
            resolution and the type of maps used are defined in the variables
            'mapres' and 'maps' in 'ppiMapImageConfig' in the loc config file.
            User defined parameters:
                angle: float
                    The elevation angle of the pseudo-PPI
                EleTol: float
                    The tolerance between the actual radar elevation angle and
                    the nominal pseudo-PPI elevation angle.
        'PSEUDORHI_CONTOUR': Plots a pseudo-RHI countour plot
            User defined parameters:
                contour_values: list of floats or None
                    The list of contour values to plot. If None the contour
                    values are going to be obtained from the Py-ART config
                    file either with the dictionary key 'contour_values' or
                    from the minimum and maximum values of the field with an
                    assumed division of 10 levels.
                angle: float
                    The azimuth angle at which to compute the RPI
                AziTol: float
                    The tolerance between the actual radar azimuth angle and
                    the nominal pseudo-RHI azimuth angle.
        'PSEUDORHI_CONTOUR_OVERPLOT': Plots a pseudo-RHI of a field with
            another field over-plotted as a contour plot
            User defined parameters:
                contour_values: list of floats or None
                    The list of contour values to plot. If None the contour
                    values are going to be obtained from the Py-ART config
                    file either with the dictionary key 'contour_values' or
                    from the minimum and maximum values of the field with an
                    assumed division of 10 levels.
                angle: float
                    The azimuth angle at which to compute the RPI
                AziTol: float
                    The tolerance between the actual radar azimuth angle and
                    the nominal pseudo-RHI azimuth angle.
        'PSEUDORHI_IMAGE': Plots a pseudo-RHI image. It can also plot the
            histogram and the quantiles of the data in the pseudo-RHI.
            User defined parameters:
                angle: float
                    The azimuth angle at which to compute the RPI
                AziTol: float
                    The tolerance between the actual radar azimuth angle and
                    the nominal pseudo-RHI azimuth angle.
                plot_type: str
                    The type of plot to perform. Can be 'RHI', 'QUANTILES' or
                    'HISTOGRAM'
                step: float or None
                    If the plot type is 'HISTOGRAM', the width of the
                    histogram bin. If None it will be obtained from the Py-ART
                    config file
                quantiles: list of float or None
                    If the plot type is 'QUANTILES', the list of quantiles to
                    compute. If None a default list of quantiles will be
                    computed
                vmin, vmax : float or None
                    Min and Max values of the color scale. If None the values
                    are taken from the Py-ART config file
        'QUANTILES': Plots and writes the quantiles of a radar volume
            User defined parameters:
                quantiles: list of floats or None
                    the list of quantiles to compute. If None a default list
                    of quantiles will be computed.
                write_data: Bool
                    If True the computed data will be also written in a csv
                    file
                fixed_span: Bool
                    If true the quantile plot has a fix Y-axis
                vmin, vmax: float or None
                    If fixed_span is set, the minimum and maximum values of
                    the Y-axis. If None, they are obtained from the Py-ART
                    config file
        'RHI_CONTOUR': Plots an RHI countour plot
            User defined parameters:
                contour_values: list of floats or None
                    The list of contour values to plot. If None the contour
                    values are going to be obtained from the Py-ART config
                    file either with the dictionary key 'contour_values' or
                    from the minimum and maximum values of the field with an
                    assumed division of 10 levels.
                anglenr: int
                    The azimuth angle number
        'RHI_CONTOUR_OVERPLOT': Plots an RHI of a field with another field
            over-plotted as a contour plot
            User defined parameters:
                contour_values: list of floats or None
                    The list of contour values to plot. If None the contour
                    values are going to be obtained from the Py-ART config
                    file either with the dictionary key 'contour_values' or
                    from the minimum and maximum values of the field with an
                    assumed division of 10 levels.
                anglenr: int
                    The azimuth angle number
        'RHI_IMAGE': Plots an RHI image. It can also plot the
            histogram and the quantiles of the data in the RHI.
            User defined parameters:
                anglenr: int
                    The azimuth angle number
                plot_type: str
                    The type of plot to perform. Can be 'RHI', 'QUANTILES' or
                    'HISTOGRAM'
                step: float or None
                    If the plot type is 'HISTOGRAM', the width of the
                    histogram bin. If None it will be obtained from the Py-ART
                    config file
                quantiles: list of float or None
                    If the plot type is 'QUANTILES', the list of quantiles to
                    compute. If None a default list of quantiles will be
                    computed
                vmin, vmax: float or None
                    The minimum and maximum values of the color scale. If None
                    the scale is going to be set according to the Py-ART
                    config file
        'RHI_PROFILE': Computes and plots a vertical profile statistics out of
            an RHI.
            The statistics are saved in a csv file
            User defined parameters:
                rangeStart, rangeStop: float
                    The range start and stop of the data to extract from the
                    RHI to compute the statistics [m]. Default 0., 25000.
                heightResolution: float
                    The height resolution of the profile [m]. Default 100.
                heightMin, heightMax: float or None
                    The minimum and maximum altitude of the profile [m MSL].
                    If None the values will be obtained from the minimum and
                    maximum gate altitude.
                quantity: str
                    The type of statistics to plot. Can be 'quantiles',
                    'mode', 'reqgression_mean' or 'mean'.
                quantiles: list of floats
                    If quantity type is 'quantiles' the list of quantiles to
                    compute. Default 25., 50., 75.
                nvalid_min: int
                    The minimum number of valid points to consider the
                    statistic valid. Default 4
                make_linear: Bool
                    If true the data is converted from log to linear before
                    computing the stats
                include_nans: Bool
                    If true NaN values are included in the statistics
                fixed_span: Bool
                    If true the profile plot has a fix X-axis
                vmin, vmax: float or None
                    If fixed_span is set, the minimum and maximum values of
                    the X-axis. If None, they are obtained from the Py-ART
                    config file
        'SAVEALL': Saves radar volume data including all or a list of user-
            defined fields in a C/F radial or ODIM file
            User defined parameters:
                file_type: str
                    The type of file used to save the data. Can be 'nc' or
                    'h5'. Default 'nc'
                datatypes: list of str or None
                    The list of data types to save. If it is None, all fields
                    in the radar object will be saved
                physical: Bool
                    If True the data will be saved in physical units (floats).
                    Otherwise it will be quantized and saved as binary
                compression: str
                    For ODIM file formats, the type of compression. Can be any
                    of the allowed compression types for hdf5 files. Default
                    gzip
                compression_opts: any
                    The compression options allowed by the hdf5. Depends on
                    the type of compression. Default 6 (The gzip compression
                    level).
        'SAVEALL_VOL' : Same as before but can be used in a mixed GRID/VOL
            dataset, as there is no ambiguity with SAVEALL for VOL datasets
        'SAVESTATE': Saves the last processed data in a file. Used for real-
            time data processing
        'SAVEVOL': Saves one field of a radar volume data in a C/F radial or
            ODIM file
            User defined parameters:
                file_type: str
                    The type of file used to save the data. Can be 'nc' or
                    'h5'. Default 'nc'
                physical: Bool
                    If True the data will be saved in physical units (floats).
                    Otherwise it will be quantized and saved as binary
                compression: str
                    For ODIM file formats, the type of compression. Can be any
                    of the allowed compression types for hdf5 files. Default
                    gzip
                compression_opts: any
                    The compression options allowed by the hdf5. Depends on
                    the type of compression. Default 6 (The gzip compression
                    level).
        'SAVEVOL_VOL' : Same as before but can be used in a mixed GRID/VOL
            dataset, as there is no ambiguity with SAVEVOL for GRID datasets
        'SAVE_FIXED_ANGLE': Saves the position of the first fix angle in a
            csv file
        'SELFCONSISTENCY': Plots a ZDR versus KDP/ZH histogram of data.
            User defined parameters:
                retrieve_relation : bool
                    If True plots also the retrieved relationship. Default
                    True
                plot_theoretical : bool
                    If True plots also the theoretical relationship. Default
                    True
                normalize : bool
                    If True the occurrence density of ZK/KDP for each ZDR bin
                    is going to be represented. Otherwise it will show the
                    number of gates at each bin. Default True
        'SELFCONSISTENCY2': Plots a ZH measured versus ZH inferred from a self-consistency
            relation histogram of data.
            User defined parameters:
                normalize : bool
                    If True the occurrence density of ZK/KDP for each ZDR bin
                    is going to be represented. Otherwise it will show the
                    number of gates at each bin. Default True
        'TIME_RANGE': Plots a time-range/azimuth/elevation plot
            User defined parameters:
                anglenr: float
                    The number of the fixed angle to plot
                vmin, vmax: float or None
                    The minimum and maximum values of the color scale. If None
                    the scale is going to be set according to the Py-ART
                    config file
        'VOL_TS': Writes and plots a value corresponding to a time series.
            Meant primarily for writing and plotting the results of the
            SELFCONSISTENCY2 algorithm
            User defined parameters:
                ref_value: float
                    The reference value. Default 0
                sort_by_date: Bool
                    If true when reading the csv file containing the
                    statistics the data is sorted by date. Default False
                rewrite: Bool
                    If true the csv file containing the statistics is
                    rewritten
                add_data_in_fname: Bool
                    If true and the data used is cumulative the year is
                    written in the csv file name and the plot file name
                npoints_min: int
                    Minimum number of points to use the data point in the
                    plotting and to send an alarm. Default 0
                vmin, vmax: float or None
                    Limits of the Y-axis (data value). If None the limits
                    are obtained from the Py-ART config file
                alarm: Bool
                    If true an alarm is sent
                tol_abs: float
                    Margin of tolerance from the reference value. If the
                    current value is above this margin an alarm is sent. If
                    the margin is not specified it is not possible to send any
                    alarm
                tol_trend: float
                    Margin of tolerance from the reference value. If the
                    trend of the last X events is above this margin an alarm
                    is sent. If the margin is not specified it is not possible
                    to send any alarm
                nevents_min: int
                    Minimum number of events with sufficient points to send an
                    alarm related to the trend. If not specified it is not
                    possible to send any alarm
                sender: str
                    The mail of the alarm sender. If not specified it is not
                    possible to send any alarm
                receiver_list: list of str
                    The list of emails of the people that will receive the
                    alarm.. If not specified it is not possible to send any
                    alarm
        'WIND_PROFILE': Plots vertical profile of wind data (U, V, W
            components and wind velocity and direction) out of a radar
            volume containing the retrieved U,V and W components of the wind,
            the standard deviation of the retrieval and the velocity
            difference between the estimated radial velocity (assuming the
            wind to be uniform) and the actual measured radial velocity.
            User defined parameters:
                heightResolution: float
                    The height resolution of the profile [m]. Default 100.
                heightMin, heightMax: float or None
                    The minimum and maximum altitude of the profile [m MSL].
                    If None the values will be obtained from the minimum and
                    maximum gate altitude.
                min_ele: float
                    The minimum elevation to be used in the computation of the
                    vertical velocities. Default 5.
                max_ele: float
                    The maximum elevation to be used in the computation of the
                    horizontal velocities. Default 85.
                fixed_span: Bool
                    If true the profile plot has a fix X-axis
                vmin, vmax: float or None
                    If fixed_span is set, the minimum and maximum values of
                    the X-axis. If None, they are obtained from the span of
                    the U component defined in the Py-ART config file

    Parameters
    ----------
    dataset : dict
        dictionary with key radar_out containing a radar object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    The list of created fields or None

    """

    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    prdsavedir = prdcfg['prdname']
    if 'prdsavedir' in prdcfg:
        prdsavedir = prdcfg['prdsavedir']

    if prdcfg['type'] == 'PPI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        el_vec = np.sort(dataset['radar_out'].fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(dataset['radar_out'].fixed_angle['data'] == el)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'ppi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        step = prdcfg.get('step', None)
        quantiles = prdcfg.get('quantiles', None)
        plot_type = prdcfg.get('plot_type', 'PPI')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        fname_list = plot_ppi(
            dataset['radar_out'], field_name, ind_el, prdcfg, fname_list,
            plot_type=plot_type, vmin=vmin, vmax=vmax, step=step,
            quantiles=quantiles)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PSEUDOPPI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        plot_type = prdcfg.get('plot_type', 'PPI')
        step = prdcfg.get('step', None)
        quantiles = prdcfg.get('quantiles', None)
        ele_tol = prdcfg.get('EleTol', 1.)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        try:
            xsect = pyart.util.cross_section_rhi(
                dataset['radar_out'], [prdcfg['angle']], el_tol=ele_tol)

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdsavedir, timeinfo=prdcfg['timeinfo'])

            fname_list = make_filename(
                'ppi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['imgformat'],
                prdcfginfo='el'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

            for i, fname in enumerate(fname_list):
                fname_list[i] = savedir+fname

            plot_ppi(
                xsect, field_name, 0, prdcfg, fname_list, plot_type=plot_type,
                vmin=vmin, vmax=vmax, step=step, quantiles=quantiles)

            print('----- save to '+' '.join(fname_list))

            return fname_list
        except EnvironmentError:
            warn(
                'No data found at elevation ' + str(prdcfg['angle']) +
                '. Skipping product ' + prdcfg['type'])

            return None

    if prdcfg['type'] == 'PPI_MAP':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        el_vec = np.sort(dataset['radar_out'].fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(
            dataset['radar_out'].fixed_angle['data'] == el)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'ppi_map', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_ppi_map(
            dataset['radar_out'], field_name, ind_el, prdcfg, fname_list)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PPI_CONTOUR_OVERPLOT':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_name = get_fieldname_pyart(prdcfg['contourtype'])
        if contour_name not in dataset['radar_out'].fields:
            warn(
                'Contour type ' + contour_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_values = prdcfg.get('contour_values', None)

        el_vec = np.sort(dataset['radar_out'].fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(
            dataset['radar_out'].fixed_angle['data'] == el)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'ppi', prdcfg['dstype'],
            prdcfg['voltype']+'-'+prdcfg['contourtype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        titl = (
            pyart.graph.common.generate_title(
                dataset['radar_out'], field_name, ind_el) +
            ' - ' +
            pyart.graph.common.generate_field_name(
                dataset['radar_out'], contour_name))

        fig, ax = plot_ppi(
            dataset['radar_out'], field_name, ind_el, prdcfg, fname_list,
            titl=titl, save_fig=False)

        fname_list = plot_ppi_contour(
            dataset['radar_out'], contour_name, ind_el, prdcfg, fname_list,
            contour_values=contour_values, ax=ax, fig=fig, save_fig=True)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PPIMAP_ROI_OVERPLOT':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        el_vec = np.sort(dataset['radar_out'].fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(
            dataset['radar_out'].fixed_angle['data'] == el)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'ppi_map', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        fig, ax, _ = plot_ppi_map(
            dataset['radar_out'], field_name, ind_el, prdcfg, fname_list,
            save_fig=False)

        fname_list = plot_roi_contour(
            dataset['roi_dict'], prdcfg, fname_list, ax=ax, fig=fig,
            save_fig=True)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PSEUDOPPI_CONTOUR_OVERPLOT':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_name = get_fieldname_pyart(prdcfg['contourtype'])
        if contour_name not in dataset['radar_out'].fields:
            warn(
                'Contour type ' + contour_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_values = prdcfg.get('contour_values', None)

        try:
            xsect = pyart.util.cross_section_rhi(
                dataset['radar_out'], [prdcfg['angle']],
                el_tol=prdcfg['EleTol'])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdsavedir, timeinfo=prdcfg['timeinfo'])

            fname_list = make_filename(
                'ppi', prdcfg['dstype'],
                prdcfg['voltype']+'-'+prdcfg['contourtype'],
                prdcfg['imgformat'],
                prdcfginfo='el'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

            for i, fname in enumerate(fname_list):
                fname_list[i] = savedir+fname

            titl = (
                pyart.graph.common.generate_title(xsect, field_name, 0) +
                ' - ' +
                pyart.graph.common.generate_field_name(xsect, contour_name))

            fig, ax = plot_ppi(xsect, field_name, 0, prdcfg, fname_list,
                               titl=titl, save_fig=False)

            fname_list = plot_ppi_contour(
                xsect, contour_name, 0, prdcfg, fname_list,
                contour_values=contour_values, ax=ax, fig=fig, save_fig=True)

            print('----- save to '+' '.join(fname_list))

            return fname_list
        except EnvironmentError:
            warn(
                'No data found at elevation ' + str(prdcfg['angle']) +
                '. Skipping product ' + prdcfg['type'])

            return None

    if prdcfg['type'] == 'PPI_CONTOUR':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_values = prdcfg.get('contour_values', None)

        el_vec = np.sort(dataset['radar_out'].fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(dataset['radar_out'].fixed_angle['data'] == el)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'ppi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        fname_list = plot_ppi_contour(
            dataset['radar_out'], field_name, ind_el, prdcfg, fname_list,
            contour_values=contour_values)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PSEUDOPPI_CONTOUR':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_values = prdcfg.get('contour_values', None)

        try:
            xsect = pyart.util.cross_section_rhi(
                dataset['radar_out'], [prdcfg['angle']],
                el_tol=prdcfg['EleTol'])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdsavedir, timeinfo=prdcfg['timeinfo'])

            fname_list = make_filename(
                'ppi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['imgformat'],
                prdcfginfo='el'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

            for i, fname in enumerate(fname_list):
                fname_list[i] = savedir+fname

            fname_list = plot_ppi_contour(
                xsect, field_name, 0, prdcfg, fname_list,
                contour_values=contour_values)

            print('----- save to '+' '.join(fname_list))

            return fname_list
        except EnvironmentError:
            warn(
                'No data found at elevation ' + str(prdcfg['angle']) +
                '. Skipping product ' + prdcfg['type'])

            return None

    if prdcfg['type'] == 'RHI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        az_vec = np.sort(dataset['radar_out'].fixed_angle['data'])
        az = az_vec[prdcfg['anglenr']]
        ind_az = np.where(dataset['radar_out'].fixed_angle['data'] == az)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'rhi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='az'+'{:.1f}'.format(az),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        step = prdcfg.get('step', None)
        quantiles = prdcfg.get('quantiles', None)
        plot_type = prdcfg.get('plot_type', 'RHI')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        plot_rhi(
            dataset['radar_out'], field_name, ind_az, prdcfg, fname_list,
            plot_type=plot_type, step=step, vmin=vmin, vmax=vmax,
            quantiles=quantiles)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PSEUDORHI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        plot_type = prdcfg.get('plot_type', 'RHI')
        step = prdcfg.get('step', None)
        quantiles = prdcfg.get('quantiles', None)
        azi_tol = prdcfg.get('AziTol', 1.)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        try:
            xsect = pyart.util.cross_section_ppi(
                dataset['radar_out'], [prdcfg['angle']], az_tol=azi_tol)

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdsavedir, timeinfo=prdcfg['timeinfo'])

            fname_list = make_filename(
                'rhi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['imgformat'],
                prdcfginfo='az'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

            for i, fname in enumerate(fname_list):
                fname_list[i] = savedir+fname

            plot_rhi(
                xsect, field_name, 0, prdcfg, fname_list, vmin=vmin,
                vmax=vmax, plot_type=plot_type, step=step,
                quantiles=quantiles)

            print('----- save to '+' '.join(fname_list))

            return fname_list
        except EnvironmentError:
            warn(
                ' No data found at azimuth ' +
                str(prdcfg['angle'])+'. Skipping product ' +
                prdcfg['type'])
            return None

    if prdcfg['type'] == 'RHI_CONTOUR_OVERPLOT':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_name = get_fieldname_pyart(prdcfg['contourtype'])
        if contour_name not in dataset['radar_out'].fields:
            warn(
                'Contour type ' + contour_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_values = prdcfg.get('contour_values', None)

        az_vec = np.sort(dataset['radar_out'].fixed_angle['data'])
        az = az_vec[prdcfg['anglenr']]
        ind_az = np.where(dataset['radar_out'].fixed_angle['data'] == az)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'rhi', prdcfg['dstype'],
            prdcfg['voltype']+'-'+prdcfg['contourtype'],
            prdcfg['imgformat'], prdcfginfo='az'+'{:.1f}'.format(az),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        titl = (
            pyart.graph.common.generate_title(dataset['radar_out'],
                                              field_name, ind_az)+' - ' +
            pyart.graph.common.generate_field_name(dataset['radar_out'],
                                                   contour_name))

        fig, ax = plot_rhi(
            dataset['radar_out'], field_name, ind_az, prdcfg, fname_list,
            titl=titl, save_fig=False)

        fname_list = plot_rhi_contour(
            dataset['radar_out'], contour_name, ind_az, prdcfg, fname_list,
            contour_values=contour_values, ax=ax, fig=fig, save_fig=True)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PSEUDORHI_CONTOUR_OVERPLOT':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_name = get_fieldname_pyart(prdcfg['contourtype'])
        if contour_name not in dataset['radar_out'].fields:
            warn(
                'Contour type ' + contour_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_values = prdcfg.get('contour_values', None)

        try:
            xsect = pyart.util.cross_section_ppi(
                dataset['radar_out'], [prdcfg['angle']],
                az_tol=prdcfg['AziTol'])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdsavedir, timeinfo=prdcfg['timeinfo'])

            fname_list = make_filename(
                'rhi', prdcfg['dstype'],
                prdcfg['voltype']+'-'+prdcfg['contourtype'],
                prdcfg['imgformat'],
                prdcfginfo='az'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

            for i, fname in enumerate(fname_list):
                fname_list[i] = savedir+fname

            titl = (
                pyart.graph.common.generate_title(xsect, field_name, 0) +
                ' - ' +
                pyart.graph.common.generate_field_name(xsect, contour_name))

            fig, ax = plot_rhi(xsect, field_name, 0, prdcfg, fname_list,
                               titl=titl, save_fig=False)

            fname_list = plot_rhi_contour(
                xsect, contour_name, 0, prdcfg, fname_list,
                contour_values=contour_values, ax=ax, fig=fig, save_fig=True)

            print('----- save to '+' '.join(fname_list))

            return fname_list
        except EnvironmentError:
            warn(
                ' No data found at azimuth ' +
                str(prdcfg['angle'])+'. Skipping product ' +
                prdcfg['type'])
            return None

    if prdcfg['type'] == 'RHI_CONTOUR':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_values = prdcfg.get('contour_values', None)

        az_vec = np.sort(dataset['radar_out'].fixed_angle['data'])
        az = az_vec[prdcfg['anglenr']]
        ind_az = np.where(
            dataset['radar_out'].fixed_angle['data'] == az)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'rhi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='az'+'{:.1f}'.format(az),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        fname_list = plot_rhi_contour(
            dataset['radar_out'], field_name, ind_az, prdcfg, fname_list,
            contour_values=contour_values)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PSEUDORHI_CONTOUR':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        contour_values = prdcfg.get('contour_values', None)

        try:
            xsect = pyart.util.cross_section_ppi(
                dataset['radar_out'], [prdcfg['angle']],
                az_tol=prdcfg['AziTol'])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdsavedir, timeinfo=prdcfg['timeinfo'])

            fname_list = make_filename(
                'rhi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['imgformat'],
                prdcfginfo='az'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

            for i, fname in enumerate(fname_list):
                fname_list[i] = savedir+fname

            fname_list = plot_rhi_contour(
                xsect, field_name, 0, prdcfg, fname_list,
                contour_values=contour_values)

            print('----- save to '+' '.join(fname_list))

            return fname_list
        except EnvironmentError:
            warn(
                ' No data found at azimuth ' +
                str(prdcfg['angle'])+'. Skipping product ' +
                prdcfg['type'])
            return None

    if prdcfg['type'] == 'RHI_PROFILE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined parameters
        rangeStart = prdcfg.get('rangeStart', 0.)
        rangeStop = prdcfg.get('rangeStop', 25000.)
        heightResolution = prdcfg.get('heightResolution', 500.)
        hmin_user = prdcfg.get('heightMin', None)
        hmax_user = prdcfg.get('heightMax', None)
        quantity = prdcfg.get('quantity', 'quantiles')
        quantiles = prdcfg.get('quantiles', np.array([25., 50., 75.]))
        nvalid_min = prdcfg.get('nvalid_min', 4)
        make_linear = prdcfg.get('make_linear', 0)
        include_nans = prdcfg.get('include_nans', 0)

        fixed_span = prdcfg.get('fixed_span', 1)
        vmin = None
        vmax = None
        if fixed_span:
            vmin, vmax = pyart.config.get_field_limits(field_name)
            if 'vmin' in prdcfg:
                vmin = prdcfg['vmin']
            if 'vmax' in prdcfg:
                vmax = prdcfg['vmax']

        # create new radar object with only data for the given rhi and range
        az_vec = np.sort(dataset['radar_out'].fixed_angle['data'])
        az = az_vec[prdcfg['anglenr']]
        ind_az = np.where(dataset['radar_out'].fixed_angle['data'] == az)[0][0]

        new_dataset = dataset['radar_out'].extract_sweeps([ind_az])
        field = new_dataset.fields[field_name]
        rng_mask = np.logical_and(new_dataset.range['data'] >= rangeStart,
                                  new_dataset.range['data'] <= rangeStop)
        field['data'] = field['data'][:, rng_mask]
        new_dataset.range['data'] = new_dataset.range['data'][rng_mask]
        new_dataset.ngates = len(new_dataset.range['data'])
        new_dataset.init_gate_x_y_z()
        new_dataset.init_gate_longitude_latitude()
        new_dataset.init_gate_altitude()

        new_dataset.fields = dict()
        new_dataset.add_field(field_name, field)

        # compute quantities
        if hmin_user is None:
            minheight = (
                round(np.min(dataset['radar_out'].gate_altitude['data']) /
                      heightResolution)*heightResolution-heightResolution)
        else:
            minheight = hmin_user
        if hmax_user is None:
            maxheight = (
                round(np.max(dataset['radar_out'].gate_altitude['data']) /
                      heightResolution)*heightResolution+heightResolution)
        else:
            maxheight = hmax_user
        nlevels = int((maxheight-minheight)/heightResolution)

        h_vec = (
            minheight+np.arange(nlevels)*heightResolution+heightResolution/2.)
        vals, val_valid = compute_profile_stats(
            field['data'], new_dataset.gate_altitude['data'], h_vec,
            heightResolution, quantity=quantity, quantiles=quantiles/100.,
            nvalid_min=nvalid_min, make_linear=make_linear,
            include_nans=include_nans)

        # plot data
        if quantity == 'mean':
            data = [vals[:, 0], vals[:, 1], vals[:, 2]]
            labels = ['Mean', 'Min', 'Max']
            colors = ['b', 'k', 'k']
            linestyles = ['-', '--', '--']
        elif quantity == 'mode':
            data = [vals[:, 0], vals[:, 2], vals[:, 4]]
            labels = ['Mode', '2nd most common', '3rd most common']
            colors = ['b', 'k', 'r']
            linestyles = ['-', '--', '--']
        else:
            data = [vals[:, 1], vals[:, 0], vals[:, 2]]
            labels = [
                str(quantiles[1])+'-percentile',
                str(quantiles[0])+'-percentile',
                str(quantiles[2])+'-percentile']
            colors = ['b', 'k', 'k']
            linestyles = ['-', '--', '--']

        labelx = get_colobar_label(
            dataset['radar_out'].fields[field_name], field_name)
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(dataset['radar_out'].fields[field_name],
                           field_name))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'az'+'{:.1f}'.format(az)+'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'rhi_profile', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles)

        print('----- save to '+' '.join(fname_list))

        fname = make_filename(
            'rhi_profile', prdcfg['dstype'], prdcfg['voltype'],
            ['csv'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        if quantity == 'mode':
            data.append(vals[:, 1])
            labels.append('% points mode')
            data.append(vals[:, 3])
            labels.append('% points 2nd most common')
            data.append(vals[:, 5])
            labels.append('% points 3rd most common')

        sector = {
            'rmin': rangeStart,
            'rmax': rangeStop,
            'az': az
        }
        write_rhi_profile(
            h_vec, data, val_valid, labels, fname, datatype=labelx,
            timeinfo=prdcfg['timeinfo'], sector=sector)

        print('----- save to '+fname)

        # TODO: add Cartesian interpolation option

        return fname

    if prdcfg['type'] == 'PROFILE_STATS':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # mask unclassified data
        field = deepcopy(dataset['radar_out'].fields[field_name]['data'])
        if prdcfg['voltype'] in ('hydro', 'hydroc'):
            field = np.ma.masked_equal(field, 1)

        # user defined parameters
        heightResolution = prdcfg.get('heightResolution', 100.)
        hmin_user = prdcfg.get('heightMin', None)
        hmax_user = prdcfg.get('heightMax', None)
        quantity = prdcfg.get('quantity', 'quantiles')
        quantiles = prdcfg.get('quantiles', np.array([25., 50., 75.]))
        nvalid_min = prdcfg.get('nvalid_min', 4)
        make_linear = prdcfg.get('make_linear', 0)
        include_nans = prdcfg.get('include_nans', 0)

        fixed_span = prdcfg.get('fixed_span', 1)
        vmin = None
        vmax = None
        if fixed_span:
            vmin, vmax = pyart.config.get_field_limits(field_name)
            if 'vmin' in prdcfg:
                vmin = prdcfg['vmin']
            if 'vmax' in prdcfg:
                vmax = prdcfg['vmax']

        # compute quantities
        if hmin_user is None:
            minheight = (round(
                np.min(dataset['radar_out'].gate_altitude['data']) /
                heightResolution)*heightResolution-heightResolution)
        else:
            minheight = hmin_user
        if hmax_user is None:
            maxheight = (round(
                np.max(dataset['radar_out'].gate_altitude['data']) /
                heightResolution)*heightResolution+heightResolution)
        else:
            maxheight = hmax_user
        nlevels = int((maxheight-minheight)/heightResolution)

        h_vec = (
            minheight+np.arange(nlevels)*heightResolution+heightResolution/2.)
        vals, val_valid = compute_profile_stats(
            field, dataset['radar_out'].gate_altitude['data'], h_vec,
            heightResolution, quantity=quantity, quantiles=quantiles/100.,
            nvalid_min=nvalid_min, make_linear=make_linear,
            include_nans=include_nans)

        # plot data
        if quantity == 'mean':
            data = [vals[:, 0], vals[:, 1], vals[:, 2]]
            labels = ['Mean', 'Min', 'Max']
            colors = ['b', 'k', 'k']
            linestyles = ['-', '--', '--']
        elif quantity == 'mode':
            data = [vals[:, 0], vals[:, 2], vals[:, 4]]
            labels = ['Mode', '2nd most common', '3rd most common']
            colors = ['b', 'k', 'r']
            linestyles = ['-', '--', '--']
        else:
            data = [vals[:, 1], vals[:, 0], vals[:, 2]]
            labels = [
                str(quantiles[1])+'-percentile',
                str(quantiles[0])+'-percentile',
                str(quantiles[2])+'-percentile']
            colors = ['b', 'k', 'k']
            linestyles = ['-', '--', '--']

        labelx = get_colobar_label(
            dataset['radar_out'].fields[field_name], field_name)
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(dataset['radar_out'].fields[field_name],
                           field_name))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'rhi_profile', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=vmin, vmax=vmax,
            hmin=minheight, hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        fname = make_filename(
            'rhi_profile', prdcfg['dstype'], prdcfg['voltype'],
            ['csv'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        if quantity == 'mode':
            data.append(vals[:, 1])
            labels.append('% points mode')
            data.append(vals[:, 3])
            labels.append('% points 2nd most common')
            data.append(vals[:, 5])
            labels.append('% points 3rd most common')

        write_rhi_profile(
            h_vec, data, val_valid, labels, fname, datatype=labelx,
            timeinfo=prdcfg['timeinfo'])

        print('----- save to '+fname)

        # TODO: add Cartesian interpolation option

        return fname

    if prdcfg['type'] == 'WIND_PROFILE':
        # user defined parameters
        heightResolution = prdcfg.get('heightResolution', 100.)
        hmin_user = prdcfg.get('heightMin', None)
        hmax_user = prdcfg.get('heightMax', None)
        min_ele = prdcfg.get('min_ele', 5.)
        max_ele = prdcfg.get('max_ele', 85.)

        fixed_span = prdcfg.get('fixed_span', 1)
        vmin = None
        vmax = None
        if fixed_span:
            vmin, vmax = pyart.config.get_field_limits(
                'eastward_wind_component')
            if 'vmin' in prdcfg:
                vmin = prdcfg['vmin']
            if 'vmax' in prdcfg:
                vmax = prdcfg['vmax']

        u_vel = deepcopy(
            dataset['radar_out'].fields['eastward_wind_component']['data'])
        v_vel = deepcopy(
            dataset['radar_out'].fields['northward_wind_component']['data'])
        w_vel = deepcopy(
            dataset['radar_out'].fields['vertical_wind_component']['data'])
        std_vel = deepcopy(
            dataset['radar_out'].fields['retrieved_velocity_std']['data'])
        diff_vel = deepcopy(
            dataset['radar_out'].fields['velocity_difference']['data'])

        # remove azimuth information
        u_vel_aux = np.ma.masked_all(
            (dataset['radar_out'].nsweeps, dataset['radar_out'].ngates),
            dtype=float)
        v_vel_aux = np.ma.masked_all(
            (dataset['radar_out'].nsweeps, dataset['radar_out'].ngates),
            dtype=float)
        w_vel_aux = np.ma.masked_all(
            (dataset['radar_out'].nsweeps, dataset['radar_out'].ngates),
            dtype=float)
        std_vel_aux = np.ma.masked_all(
            (dataset['radar_out'].nsweeps, dataset['radar_out'].ngates),
            dtype=float)
        ngates_aux = np.zeros(
            (dataset['radar_out'].nsweeps, dataset['radar_out'].ngates),
            dtype=int)
        gate_altitude_aux = np.empty(
            (dataset['radar_out'].nsweeps, dataset['radar_out'].ngates),
            dtype=float)
        for ind_sweep in range(dataset['radar_out'].nsweeps):
            ind_start = (
                dataset['radar_out'].sweep_start_ray_index['data'][ind_sweep])
            ind_end = (
                dataset['radar_out'].sweep_end_ray_index['data'][ind_sweep])
            for ind_rng in range(dataset['radar_out'].ngates):
                u_vel_aux[ind_sweep, ind_rng] = u_vel[ind_start, ind_rng]
                v_vel_aux[ind_sweep, ind_rng] = v_vel[ind_start, ind_rng]
                w_vel_aux[ind_sweep, ind_rng] = w_vel[ind_start, ind_rng]
                std_vel_aux[ind_sweep, ind_rng] = std_vel[ind_start, ind_rng]
                gate_altitude_aux[ind_sweep, ind_rng] = (
                    dataset['radar_out'].gate_altitude['data'][
                        ind_start, ind_rng])
                ngates_aux[ind_sweep, ind_rng] = (
                    diff_vel[ind_start:ind_end, ind_rng].compressed().size)

        # exclude low elevations in the computation of vertical velocities
        std_w_vel_aux = deepcopy(std_vel_aux)
        ngates_w_aux = deepcopy(ngates_aux)
        ind = np.where(dataset['radar_out'].fixed_angle['data'] < min_ele)[0]
        if ind.size > 0:
            w_vel_aux[ind, :] = np.ma.masked
            std_w_vel_aux[ind, :] = np.ma.masked
            ngates_w_aux[ind, :] = 0

        # exclude hig elevations in the computation of horizontal velocities
        ind = np.where(dataset['radar_out'].fixed_angle['data'] > max_ele)[0]
        if ind.size > 0:
            u_vel_aux[ind, :] = np.ma.masked
            v_vel_aux[ind, :] = np.ma.masked
            std_vel_aux[ind, :] = np.ma.masked
            ngates_aux[ind, :] = 0

        # compute quantities
        if hmin_user is None:
            minheight = (round(
                np.min(dataset['radar_out'].gate_altitude['data']) /
                heightResolution)*heightResolution-heightResolution)
        else:
            minheight = hmin_user
        if hmax_user is None:
            maxheight = (round(
                np.max(dataset['radar_out'].gate_altitude['data']) /
                heightResolution)*heightResolution+heightResolution)
        else:
            maxheight = hmax_user
        nlevels = int((maxheight-minheight)/heightResolution)

        h_vec = (
            minheight+np.arange(nlevels)*heightResolution+heightResolution/2.)

        u_vals, val_valid = compute_profile_stats(
            u_vel_aux, gate_altitude_aux, h_vec, heightResolution,
            quantity='regression_mean', std_field=std_vel_aux,
            np_field=ngates_aux)
        v_vals, val_valid = compute_profile_stats(
            v_vel_aux, gate_altitude_aux, h_vec, heightResolution,
            quantity='regression_mean', std_field=std_vel_aux,
            np_field=ngates_aux)
        w_vals, w_val_valid = compute_profile_stats(
            w_vel_aux, gate_altitude_aux, h_vec, heightResolution,
            quantity='regression_mean', std_field=std_w_vel_aux,
            np_field=ngates_w_aux)

        # plot u wind data
        u_data = [u_vals[:, 0], u_vals[:, 0]+u_vals[:, 1],
                  u_vals[:, 0]-u_vals[:, 1]]
        labels = ['Regression mean', '+std', '-std']
        colors = ['b', 'k', 'k']
        linestyles = ['-', '--', '--']

        labelx = get_colobar_label(
            dataset['radar_out'].fields['eastward_wind_component'],
            'eastward_wind_component')
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(
                dataset['radar_out'].fields['eastward_wind_component'],
                'eastward_wind_component'))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'wind_profile', prdcfg['dstype'], 'wind_vel_h_u',
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            u_data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=vmin, vmax=vmax,
            hmin=minheight, hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        # plot v wind data
        v_data = [v_vals[:, 0], v_vals[:, 0]+v_vals[:, 1],
                  v_vals[:, 0]-v_vals[:, 1]]
        labels = ['Regression mean', '+std', '-std']
        colors = ['b', 'k', 'k']
        linestyles = ['-', '--', '--']

        labelx = get_colobar_label(
            dataset['radar_out'].fields['northward_wind_component'],
            'northward_wind_component')
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(
                dataset['radar_out'].fields['northward_wind_component'],
                'northward_wind_component'))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'wind_profile', prdcfg['dstype'], 'wind_vel_h_v',
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            v_data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=vmin, vmax=vmax,
            hmin=minheight, hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        # plot vertical wind data
        w_data = [w_vals[:, 0], w_vals[:, 0]+w_vals[:, 1],
                  w_vals[:, 0]-w_vals[:, 1]]
        labels = ['Regression mean', '+std', '-std']
        colors = ['b', 'k', 'k']
        linestyles = ['-', '--', '--']

        labelx = get_colobar_label(
            dataset['radar_out'].fields['vertical_wind_component'],
            'vertical_wind_component')
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(
                dataset['radar_out'].fields['vertical_wind_component'],
                'vertical_wind_component'))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'wind_profile', prdcfg['dstype'], 'wind_vel_v',
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            w_data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=vmin, vmax=vmax,
            hmin=minheight, hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        # plot horizontal wind magnitude
        mag = np.ma.sqrt(
            np.ma.power(u_vals[:, 0], 2.)+np.ma.power(v_vals[:, 0], 2.))
        mag_data = [mag]
        labels = ['Regression mean']
        colors = ['b']
        linestyles = ['-']

        field_dict = pyart.config.get_metadata('wind_speed')
        labelx = get_colobar_label(field_dict, 'wind_speed')
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(field_dict, 'wind_speed'))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'wind_profile', prdcfg['dstype'], 'WIND_SPEED',
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            mag_data, h_vec, fname_list, labelx=labelx,
            labely='Height (m MSL)', labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=vmin, vmax=vmax, hmin=minheight,
            hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        # plot horizontal wind direction
        wind_dir = (
            90.-np.ma.arctan2(u_vals[:, 0]/mag, v_vals[:, 0]/mag)*180. /
            np.pi+180.)
        wind_dir[wind_dir >= 360.] = wind_dir[wind_dir >= 360.]-360.
        dir_data = [wind_dir]
        labels = ['Regression mean']
        colors = ['b']
        linestyles = ['-']

        field_dict = pyart.config.get_metadata('wind_direction')
        labelx = get_colobar_label(field_dict, 'wind_direction')
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(field_dict, 'wind_direction'))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'wind_profile', prdcfg['dstype'], 'WIND_DIRECTION',
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            dir_data, h_vec, fname_list, labelx=labelx,
            labely='Height (m MSL)', labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=0., vmax=360.,
            hmin=minheight, hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        fname = make_filename(
            'wind_profile', prdcfg['dstype'], 'WIND',
            ['csv'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        data = [
            u_vals[:, 0], u_vals[:, 1], np.ma.asarray(val_valid),
            v_vals[:, 0], v_vals[:, 1], np.ma.asarray(val_valid),
            w_vals[:, 0], w_vals[:, 1], np.ma.asarray(w_val_valid),
            mag, wind_dir]
        labels = [
            'u_wind', 'std_u_wind', 'np_u_wind',
            'v_wind', 'std_v_wind', 'np_v_wind',
            'w_wind', 'std_w_wind', 'np_w_wind',
            'mag_h_wind', 'dir_h_wind']

        write_rhi_profile(
            h_vec, data, val_valid, labels, fname, datatype=None,
            timeinfo=prdcfg['timeinfo'])

        print('----- save to '+fname)

        return fname

    if prdcfg['type'] == 'PSEUDOPPI_MAP':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        try:
            xsect = pyart.util.cross_section_rhi(
                dataset['radar_out'], [prdcfg['angle']],
                el_tol=prdcfg['EleTol'])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdsavedir, timeinfo=prdcfg['timeinfo'])

            fname_list = make_filename(
                'ppi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['imgformat'],
                prdcfginfo='el'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

            for i, fname in enumerate(fname_list):
                fname_list[i] = savedir+fname

            plot_ppi_map(xsect, field_name, 0, prdcfg, fname_list)

            print('----- save to '+' '.join(fname_list))

            return fname_list
        except EnvironmentError:
            warn(
                'No data found at elevation ' + str(prdcfg['angle']) +
                '. Skipping product ' + prdcfg['type'])

            return None

    if prdcfg['type'] == 'CAPPI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'cappi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='alt'+'{:.1f}'.format(prdcfg['altitude']),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_cappi(
            dataset['radar_out'], field_name, prdcfg['altitude'], prdcfg,
            fname_list)
        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'FIXED_RNG_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined parameters
        ang_tol = prdcfg.get('AngTol', 1.)
        azi_res = prdcfg.get('azi_res', None)
        ele_res = prdcfg.get('ele_res', None)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'fixed_rng', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='rng'+'{:.1f}'.format(
                dataset['radar_out'].range['data'][0]),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_fixed_rng(
            dataset['radar_out'], field_name, prdcfg, fname_list,
            azi_res=azi_res, ele_res=ele_res, ang_tol=ang_tol, vmin=vmin,
            vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'FIXED_RNG_SPAN_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined parameters
        ang_tol = prdcfg.get('AngTol', 1.)
        azi_res = prdcfg.get('azi_res', None)
        ele_res = prdcfg.get('ele_res', None)
        stat = prdcfg.get('stat', 'max')

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            stat, prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='rng' +
            '{:.1f}'.format(dataset['radar_out'].range['data'][0])+'-' +
            '{:.1f}'.format(dataset['radar_out'].range['data'][-1]),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_fixed_rng_span(
            dataset['radar_out'], field_name, prdcfg, fname_list,
            azi_res=azi_res, ele_res=ele_res, ang_tol=ang_tol, stat=stat)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PLOT_ALONG_COORD':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        if ((dataset['radar_out'].scan_type != 'ppi') and
                (dataset['radar_out'].scan_type != 'rhi')):
            warn('This product is only available for PPI or RHI volumes')
            return None

        colors = prdcfg.get('colors', None)
        if prdcfg['mode'] == 'ALONG_RNG':
            value_start = prdcfg.get('value_start', 0.)
            value_stop = prdcfg.get(
                'value_stop', np.max(dataset['radar_out'].range['data']))
            ang_tol = prdcfg.get('AngTol', 1.)

            xvals, yvals, valid_azi, valid_ele = get_data_along_rng(
                dataset['radar_out'], field_name, prdcfg['fix_elevations'],
                prdcfg['fix_azimuths'], ang_tol=ang_tol, rmin=value_start,
                rmax=value_stop)

            if not yvals:
                warn('No data found')
                return None

            labelx = 'Range (m)'

            labels = list()
            for ele, azi in zip(valid_ele, valid_azi):
                labels.append(
                    'azi '+'{:.1f}'.format(azi)+' ele '+'{:.1f}'.format(ele))

        elif prdcfg['mode'] == 'ALONG_AZI':
            value_start = prdcfg.get(
                'value_start', np.min(dataset['radar_out'].azimuth['data']))
            value_stop = prdcfg.get(
                'value_stop', np.max(dataset['radar_out'].azimuth['data']))
            ang_tol = prdcfg.get('AngTol', 1.)
            rng_tol = prdcfg.get('RngTol', 50.)

            xvals, yvals, valid_rng, valid_ele = get_data_along_azi(
                dataset['radar_out'], field_name, prdcfg['fix_ranges'],
                prdcfg['fix_elevations'], rng_tol=rng_tol, ang_tol=ang_tol,
                azi_start=value_start, azi_stop=value_stop)

            if not yvals:
                warn('No data found')
                return None

            labelx = 'Azimuth Angle (deg)'

            labels = list()
            for ele, rng in zip(valid_ele, valid_rng):
                labels.append(
                    'rng '+'{:.1f}'.format(rng)+' ele '+'{:.1f}'.format(ele))

        elif prdcfg['mode'] == 'ALONG_ELE':
            value_start = prdcfg.get(
                'value_start', np.min(dataset['radar_out'].elevation['data']))
            value_stop = prdcfg.get(
                'value_stop', np.max(dataset['radar_out'].elevation['data']))
            ang_tol = prdcfg.get('AngTol', 1.)
            rng_tol = prdcfg.get('RngTol', 50.)

            xvals, yvals, valid_rng, valid_azi = get_data_along_ele(
                dataset['radar_out'], field_name, prdcfg['fix_ranges'],
                prdcfg['fix_azimuths'], rng_tol=rng_tol, ang_tol=ang_tol,
                ele_min=value_start, ele_max=value_stop)

            if not yvals:
                warn('No data found')
                return None

            labelx = 'Elevation Angle (deg)'

            labels = list()
            for azi, rng in zip(valid_azi, valid_rng):
                labels.append(
                    'rng '+'{:.1f}'.format(rng)+' azi '+'{:.1f}'.format(azi))
        else:
            warn('Unknown plotting mode '+prdcfg['mode'])
            return None

        labely = get_colobar_label(
            dataset['radar_out'].fields[field_name], field_name)
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(
                dataset['radar_out'].fields[field_name], field_name))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            prdcfg['mode'], prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=prdcfg['timeinfo'],
            runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_along_coord(
            xvals, yvals, fname_list, labelx=labelx, labely=labely,
            labels=labels, title=titl, colors=colors)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PLOT_TXH':
        radar = dataset['radar_out']

        field_name = 'transmitted_signal_power_h'
        if field_name not in radar.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined parameters
        azi_res = prdcfg.get('azi_res', None)
        ele_res = prdcfg.get('ele_res', None)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)
        angtol = prdcfg.get('ang_tol', 0.5)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], prdcfg['timeinfo'])

        fname_list = make_filename(
            'constr', prdcfg['dstype'], prdcfg['dsname'],
            prdcfg['imgformat'],
            prdcfginfo='rng'+'{:.1f}'.format(
                dataset['radar_out'].range['data'][0]),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_fixed_rng(radar, field_name, prdcfg, fname_list, azi_res=None,
                       ele_res=None, ang_tol=angtol, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'BSCOPE_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        ray_dim = prdcfg.get('ray_dim', 'ang')
        xaxis_rng = prdcfg.get('xaxis_rng', True)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        ang_vec = np.sort(dataset['radar_out'].fixed_angle['data'])
        ang = ang_vec[prdcfg['anglenr']]
        ind_ang = np.where(
            dataset['radar_out'].fixed_angle['data'] == ang)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'b-scope', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='ang'+'{:.1f}'.format(ang),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if dataset['radar_out'].rays_per_sweep['data'][ind_ang] > 1:
            plot_bscope(
                dataset['radar_out'], field_name, ind_ang, prdcfg, fname_list,
                vmin=vmin, vmax=vmax, ray_dim=ray_dim, xaxis_rng=xaxis_rng)
        else:
            plot_ray(
                dataset['radar_out'].extract_sweeps([ind_ang]), field_name, 0,
                prdcfg, fname_list)
        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'TIME_RANGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        ang_vec = np.sort(dataset['radar_out'].fixed_angle['data'])
        ang = ang_vec[prdcfg['anglenr']]
        ind_ang = np.where(
            dataset['radar_out'].fixed_angle['data'] == ang)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'time-range', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='ang'+'{:.1f}'.format(ang),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if dataset['radar_out'].scan_type == 'ALONG_AZI':
            ylabel = 'azimuth (deg)'
        elif dataset['radar_out'].scan_type == 'ALONG_RNG':
            ylabel = 'range (km)'
        elif dataset['radar_out'].scan_type == 'ALONG_ELE':
            ylabel = 'elevation (deg)'
        else:
            ylabel = 'range (km)'

        plot_time_range(
            dataset['radar_out'], field_name, ind_ang, prdcfg, fname_list,
            vmin=vmin, vmax=vmax, ylabel=ylabel)
        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'HISTOGRAM':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        step = prdcfg.get('step', None)
        write_data = prdcfg.get('write_data', 0)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'histogram', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        bin_edges, values = compute_histogram(
            dataset['radar_out'].fields[field_name]['data'], field_name,
            step=step)

        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(
                dataset['radar_out'].fields[field_name], field_name))

        labelx = get_colobar_label(
            dataset['radar_out'].fields[field_name], field_name)

        plot_histogram(bin_edges, values, fname_list, labelx=labelx,
                       labely='Number of Samples', titl=titl)

        print('----- save to '+' '.join(fname_list))

        if write_data:
            fname = savedir+make_filename(
                'histogram', prdcfg['dstype'], prdcfg['voltype'],
                ['csv'], timeinfo=prdcfg['timeinfo'],
                runinfo=prdcfg['runinfo'])[0]

            hist, _ = np.histogram(values, bins=bin_edges)
            write_histogram(
                bin_edges, hist, fname, datatype=prdcfg['voltype'], step=step)
            print('----- save to '+fname)

            return fname

        return fname_list

    if prdcfg['type'] == 'SAMPLES_HISTOGRAM':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        if 'samples' not in dataset[field_name]:
            warn('Unable to compute histogram. No samples')
            return None

        step = prdcfg.get('step', None)
        write_data = prdcfg.get('write_data', 0)

        timeformat = '%Y%m%d'
        if dataset[field_name]['bias_type'] == 'instant':
            timeformat = '%Y%m%d%H%M%S'
        field_metadata = pyart.config.get_metadata(field_name)
        titl = (
            dataset[field_name]['timeinfo'].strftime(timeformat) + '\n' +
            get_field_name(field_metadata, field_name))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=dataset[field_name]['timeinfo'])

        fname_list = make_filename(
            'histogram', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=dataset[field_name]['timeinfo'],
            runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        bin_edges, values = compute_histogram(
            dataset[field_name]['samples'], field_name, step=step)

        labelx = get_colobar_label(field_metadata, field_name)

        plot_histogram(bin_edges, values, fname_list, labelx=labelx,
                       labely='Number of Samples', titl=titl)

        print('----- save to '+' '.join(fname_list))

        if write_data:
            fname = savedir+make_filename(
                'histogram', prdcfg['dstype'], prdcfg['voltype'],
                ['csv'], timeinfo=prdcfg['timeinfo'],
                runinfo=prdcfg['runinfo'])[0]

            hist, _ = np.histogram(values, bins=bin_edges)
            write_histogram(
                bin_edges, hist, fname, datatype=prdcfg['voltype'], step=step)
            print('----- save to '+fname)

            return fname

        return fname_list

    if prdcfg['type'] == 'QUANTILES':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # mask unclassified data
        field = deepcopy(dataset['radar_out'].fields[field_name]['data'])
        if prdcfg['voltype'] in ('hydro', 'hydroc'):
            field = np.ma.masked_equal(field, 1)

        # user defined variables
        quantiles = prdcfg.get('quantiles', None)
        write_data = prdcfg.get('write_data', 0)

        fixed_span = prdcfg.get('fixed_span', 1)
        vmin = None
        vmax = None
        if fixed_span:
            vmin, vmax = pyart.config.get_field_limits(field_name)
            if 'vmin' in prdcfg:
                vmin = prdcfg['vmin']
            if 'vmax' in prdcfg:
                vmax = prdcfg['vmax']

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'quantiles', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        quantiles, values = compute_quantiles(field, quantiles=quantiles)

        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(
                dataset['radar_out'].fields[field_name], field_name))

        labely = get_colobar_label(
            dataset['radar_out'].fields[field_name], field_name)

        plot_quantiles(quantiles, values, fname_list, labelx='quantile',
                       labely=labely, titl=titl, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        if write_data:
            fname = savedir+make_filename(
                'quantiles', prdcfg['dstype'], prdcfg['voltype'],
                ['csv'], timeinfo=prdcfg['timeinfo'],
                runinfo=prdcfg['runinfo'])[0]

            write_quantiles(
                quantiles, values, fname, datatype=prdcfg['voltype'])
            print('----- save to '+fname)

            return fname

        return fname_list

    if prdcfg['type'] == 'WRITE_MEAN':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        ele_min = prdcfg.get('ele_min', 0.)
        ele_max = prdcfg.get('ele_max', 90.)
        azi_min = prdcfg.get('azi_min', 0.)
        azi_max = prdcfg.get('azi_max', 360.)
        rmin = prdcfg.get('rmin', 0.)
        rmax = prdcfg.get('rmax', 50.)

        field_data = dataset['radar_out'].fields[field_name]['data']
        el_data = dataset['radar_out'].elevation['data']
        az_data = dataset['radar_out'].azimuth['data']
        rng_data = dataset['radar_out'].range['data']
        units = dataset['radar_out'].fields[field_name]['units']

        ind_ele = np.where((el_data >= ele_min) & (el_data <= ele_max))
        ind_azi = np.where((az_data >= azi_min) & (az_data <= azi_max))
        ind_rng = np.where((rng_data >= rmin) & (rng_data <= rmax))

        azmin = dataset['radar_out'].get_azimuth(0)[0]
        azmax = dataset['radar_out'].get_azimuth(0)[-1]
        elmin = dataset['radar_out'].fixed_angle['data'][0]
        elmax = dataset['radar_out'].fixed_angle['data'][-1]

        field_data_process = field_data[np.intersect1d(ind_ele, ind_azi), ind_rng]
        countzero = np.size(np.where(field_data_process == 0.))

        meanval = np.nanmean(field_data_process)
        stdval = np.nanstd(field_data_process)
        medianval = np.nanmedian(field_data_process)
        ntotal = np.size(field_data_process)
        nvalid = ntotal - countzero

        text = ["Statistics of a region in a volume.",
                "  Name      : " + field_name,
                "  Azimuth   : " + str(azmin) +" - "+str(azmax)+" [deg]",
                "  Elevation : " + str(elmin) +" - "+str(elmax)+" [deg]",
                "  Range     : " + str(rmin) +" - "+str(rmax)+" [m]"
               ]

        data = {'dstype': prdcfg['dstype'],
                'unit': units,
                'time':   prdcfg['timeinfo'],
                'label': ["Mean", "Median", "Stddev", "Nsamples", "Nvalid"],
                'value':  [meanval, medianval, stdval, ntotal, nvalid]
                }

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'],
            ['csv'], timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d', runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        write_timeseries_point(fname, data, field_name, text,
                               timeformat=None, timeinfo=prdcfg['timeinfo'])

        print('----- save to '+fname)

        return fname

    if prdcfg['type'] == 'FIELD_COVERAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        threshold = prdcfg.get('threshold', None)
        nvalid_min = prdcfg.get('nvalid_min', 5.)
        ele_res = prdcfg.get('ele_res', 1.)
        azi_res = prdcfg.get('azi_res', 2.)
        ele_min = prdcfg.get('ele_min', 0.)
        ele_max = prdcfg.get('ele_max', 30.)
        ele_step = prdcfg.get('ele_step', 5.)
        ele_sect_start = prdcfg.get('ele_sect_start', None)
        ele_sect_stop = prdcfg.get('ele_sect_stop', None)
        quantiles = prdcfg.get(
            'quantiles',
            np.array([10., 20., 30., 40., 50., 60., 70., 80., 90.]))

        # get coverage per ray
        field_coverage = np.ma.masked_all(dataset['radar_out'].nrays)

        for i in range(dataset['radar_out'].nrays):
            mask = np.ma.getmaskarray(
                dataset['radar_out'].fields[field_name]['data'][i, :])
            if threshold is not None:
                ind = np.where(np.logical_and(
                    ~mask,
                    dataset['radar_out'].fields[field_name]['data'][i, :] >=
                    threshold))[0]
            else:
                ind = np.where(~mask)[0]
            if len(ind) > nvalid_min:
                field_coverage[i] = (
                    dataset['radar_out'].range['data'][ind[-1]] -
                    dataset['radar_out'].range['data'][ind[0]])

        # group coverage per elevation sectors
        nsteps = int((ele_max-ele_min)/ele_step)  # number of steps
        nele = int(ele_step/ele_res)  # number of elev per step
        ele_steps_vec = np.arange(nsteps)*ele_step+ele_min

        yval = []
        xval = []
        labels = []
        for i in range(nsteps-1):
            yval_aux = np.ma.array([])
            xval_aux = np.array([])
            for j in range(nele):
                ele_target = ele_steps_vec[i]+j*ele_res
                d_ele = np.abs(
                    dataset['radar_out'].elevation['data']-ele_target)
                ind_ele = np.where(d_ele < prdcfg['AngTol'])[0]
                if ind_ele.size == 0:
                    continue
                yval_aux = np.ma.concatenate(
                    [yval_aux, field_coverage[ind_ele]])
                xval_aux = np.concatenate(
                    [xval_aux, dataset['radar_out'].azimuth['data'][ind_ele]])
            yval.append(yval_aux)
            xval.append(xval_aux)
            labels.append('ele '+'{:.1f}'.format(ele_steps_vec[i])+'-' +
                          '{:.1f}'.format(ele_steps_vec[i+1])+' deg')

        # get mean value per azimuth for a specified elevation sector
        xmeanval = None
        ymeanval = None
        quantval = None
        labelmeanval = None
        if ele_sect_start is not None and ele_sect_stop is not None:
            ind_ele = np.where(np.logical_and(
                dataset['radar_out'].elevation['data'] >= ele_sect_start,
                dataset['radar_out'].elevation['data'] <= ele_sect_stop))
            field_coverage_sector = field_coverage[ind_ele]
            azi_sector = dataset['radar_out'].azimuth['data'][ind_ele]
            nazi = int((np.max(dataset['radar_out'].azimuth['data']) -
                        np.min(dataset['radar_out'].azimuth['data'])) /
                       azi_res+1)

            xmeanval = np.arange(nazi)*azi_res+np.min(
                dataset['radar_out'].azimuth['data'])
            ymeanval = np.ma.masked_all(nazi)
            for i in range(nazi):
                d_azi = np.abs(azi_sector-xmeanval[i])
                ind_azi = np.where(d_azi < prdcfg['AngTol'])[0]
                if ind_azi.size == 0:
                    continue
                ymeanval[i] = np.ma.mean(field_coverage_sector[ind_azi])
            labelmeanval = ('ele '+'{:.1f}'.format(ele_sect_start)+'-' +
                            '{:.1f}'.format(ele_sect_stop)+' deg mean val')

            _, quantval, _ = quantiles_weighted(
                field_coverage_sector, quantiles=quantiles/100.)

        # plot field coverage
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'coverage', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(
                dataset['radar_out'].fields[field_name], field_name))

        plot_field_coverage(
            xval, yval, fname_list, labels=labels, title=titl, ymin=0.,
            ymax=np.max(dataset['radar_out'].range['data'])+60000.,
            xmeanval=xmeanval, ymeanval=ymeanval, labelmeanval=labelmeanval)

        print('----- save to '+' '.join(fname_list))

        fname = make_filename(
            'coverage', prdcfg['dstype'], prdcfg['voltype'],
            ['csv'], timeinfo=prdcfg['timeinfo'],
            runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        if quantval is not None:
            data_type = get_colobar_label(
                dataset['radar_out'].fields[field_name], field_name)
            write_field_coverage(
                quantiles, quantval, ele_sect_start, ele_sect_stop,
                np.min(xmeanval), np.max(xmeanval), threshold, nvalid_min,
                data_type, prdcfg['timeinfo'], fname)

            print('----- save to '+fname)

        return fname

    if prdcfg['type'] == 'CDF':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        quantiles = prdcfg.get('quantiles', None)
        sector = {
            'rmin': None,
            'rmax': None,
            'azmin': None,
            'azmax': None,
            'elmin': None,
            'elmax': None,
            'hmin': None,
            'hmax': None}

        if 'sector' in prdcfg:
            if 'rmin' in prdcfg['sector']:
                sector['rmin'] = prdcfg['sector']['rmin']
            if 'rmax' in prdcfg['sector']:
                sector['rmax'] = prdcfg['sector']['rmax']
            if 'azmin' in prdcfg['sector']:
                sector['azmin'] = prdcfg['sector']['azmin']
            if 'azmax' in prdcfg['sector']:
                sector['azmax'] = prdcfg['sector']['azmax']
            if 'elmin' in prdcfg['sector']:
                sector['elmin'] = prdcfg['sector']['elmin']
            if 'elmax' in prdcfg['sector']:
                sector['elmax'] = prdcfg['sector']['elmax']
            if 'hmin' in prdcfg['sector']:
                sector['hmin'] = prdcfg['sector']['hmin']
            if 'hmax' in prdcfg['sector']:
                sector['hmax'] = prdcfg['sector']['hmax']

        vismin = prdcfg.get('vismin', None)
        absolute = prdcfg.get('absolute', False)
        use_nans = prdcfg.get('use_nans', False)
        nan_value = prdcfg.get('nan_value', 0.)
        filterclt = prdcfg.get('filterclt', False)
        filterprec = prdcfg.get('filterprec', np.array([], dtype=int))

        data = deepcopy(dataset['radar_out'].fields[field_name]['data'])

        # define region of interest
        roi_flag = get_ROI(dataset['radar_out'], field_name, sector)
        data = data[roi_flag == 1]

        ntot = np.size(roi_flag[roi_flag == 1])

        if ntot == 0:
            warn('No radar gates found in sector')
            return None

        # get number of gates with clutter and mask them
        nclut = -1
        if filterclt:
            echoID_field = get_fieldname_pyart('echoID')
            if echoID_field in dataset['radar_out'].fields:
                echoID_ROI = (
                    dataset['radar_out'].fields[echoID_field]['data'][
                        roi_flag == 1])
                nclut = len(echoID_ROI[echoID_ROI == 2])
                data[echoID_ROI == 2] = np.ma.masked

        # get number of blocked gates and filter according to visibility
        nblocked = -1
        if vismin is not None:
            vis_field = get_fieldname_pyart('VIS')
            if vis_field in dataset['radar_out'].fields:
                vis_ROI = (
                    dataset['radar_out'].fields[vis_field]['data'][
                        roi_flag == 1])
                nblocked = len(vis_ROI[vis_ROI < vismin])
                data[vis_ROI < vismin] = np.ma.masked

        # filter according to precip type
        nprec_filter = -1
        if filterprec.size > 0:
            hydro_field = get_fieldname_pyart('hydro')
            if hydro_field in dataset['radar_out'].fields:
                hydro_ROI = (
                    dataset['radar_out'].fields[hydro_field]['data'][
                        roi_flag == 1])
                nprec_filter = 0
                for ind_hydro in filterprec:
                    nprec_filter += len(hydro_ROI[hydro_ROI == ind_hydro])
                    data[hydro_ROI == ind_hydro] = np.ma.masked

        if absolute:
            data = np.ma.abs(data)

        mask = np.ma.getmaskarray(data)
        nnan = np.count_nonzero(mask)

        if nnan == ntot:
            warn('No valid radar gates found in sector')
            return None

        if use_nans:
            data[mask] = nan_value

        # count and filter outliers
        _, values_lim = compute_quantiles(data, quantiles=[0.2, 99.8])
        if values_lim.mask[0] or values_lim.mask[1]:
            warn('No valid radar gates found in sector')
            return None

        nsmall = np.count_nonzero(data.compressed() < values_lim[0])
        nlarge = np.count_nonzero(data.compressed() > values_lim[1])
        noutliers = nlarge+nsmall
        data = data[np.logical_and(
            data >= values_lim[0], data <= values_lim[1])]

        # number of values used for cdf computation
        ncdf = np.size(data.compressed())

        quantiles, values = compute_quantiles(data, quantiles=quantiles)

        # plot CDF
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'cdf', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset['radar_out']).isoformat() + 'Z' + '\n' +
            get_field_name(
                dataset['radar_out'].fields[field_name], field_name))

        labelx = get_colobar_label(
            dataset['radar_out'].fields[field_name], field_name)

        plot_quantiles(values, quantiles/100., fname_list, labelx=labelx,
                       labely='Cumulative probability', titl=titl)

        print('----- save to '+' '.join(fname_list))

        # store cdf values
        fname = make_filename(
            'cdf', prdcfg['dstype'], prdcfg['voltype'],
            ['txt'], timeinfo=prdcfg['timeinfo'],
            runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        write_cdf(
            quantiles, values, ntot, nnan, nclut, nblocked, nprec_filter,
            noutliers, ncdf, fname, use_nans=use_nans, nan_value=nan_value,
            filterprec=filterprec, vismin=vismin, sector=sector,
            datatype=labelx, timeinfo=prdcfg['timeinfo'])

        print('----- save to '+fname)

        return fname

    if prdcfg['type'] == 'SELFCONSISTENCY':
        if 'selfconsistency_points' not in dataset:
            return None

        retrieve_relation = prdcfg.get('retrieve_relation', True)
        plot_theoretical = prdcfg.get('plot_theoretical', True)
        normalize = prdcfg.get('normalize', True)

        timeinfo = dataset['selfconsistency_points']['timeinfo']
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=timeinfo)

        fname_list = make_filename(
            'selfconsistency', prdcfg['dstype'], 'selfconsistency',
            prdcfg['imgformat'], timeinfo=timeinfo, runinfo=prdcfg['runinfo'],
            timeformat='%Y%m%d')

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        fname_list = plot_selfconsistency_instrument(
            np.array(dataset['selfconsistency_points']['zdr']),
            np.array(dataset['selfconsistency_points']['kdp']),
            np.array(dataset['selfconsistency_points']['zh']), fname_list,
            parametrization=dataset['selfconsistency_points']['parametrization'],
            zdr_kdpzh_dict=dataset['selfconsistency_points']['zdr_kdpzh_dict'],
            normalize=normalize, retrieve_relation=retrieve_relation,
            plot_theoretical=plot_theoretical)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'SELFCONSISTENCY2':
        if 'selfconsistency_points' not in dataset:
            return None

        normalize = prdcfg.get('normalize', True)

        timeinfo = dataset['selfconsistency_points']['timeinfo']
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=timeinfo)

        fname_list = make_filename(
            'selfconsistency2', prdcfg['dstype'], 'selfconsistency2',
            prdcfg['imgformat'], timeinfo=timeinfo, runinfo=prdcfg['runinfo'],
            timeformat='%Y%m%d')

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        fname_list = plot_selfconsistency_instrument2(
            np.array(dataset['selfconsistency_points']['zdr']),
            np.array(dataset['selfconsistency_points']['kdp']),
            10.*np.ma.log10(np.array(dataset['selfconsistency_points']['zh'])),
            fname_list,
            parametrization=dataset['selfconsistency_points']['parametrization'],
            zdr_kdpzh_dict=dataset['selfconsistency_points']['zdr_kdpzh_dict'],
            normalize=normalize)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'VOL_TS':
        # Time series of a value
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # check the type of dataset required
        bias_type = prdcfg.get('bias_type', 'cumulative')
        if dataset[field_name]['bias_type'] != bias_type:
            return None

        # put time info in file path and name
        csvtimeinfo_path = None
        csvtimeinfo_file = None
        timeformat = None
        if bias_type == 'instant':
            csvtimeinfo_path = dataset[field_name]['timeinfo']
            csvtimeinfo_file = dataset[field_name]['timeinfo']
            timeformat = '%Y%m%d'
        if prdcfg.get('add_date_in_fname', False):
            csvtimeinfo_file = dataset[field_name]['timeinfo']
            timeformat = '%Y'

        quantiles = prdcfg.get('quantiles', np.array([25., 75.]))
        ref_value = prdcfg.get('ref_value', 0.)
        sort_by_date = prdcfg.get('sort_by_date', False)
        rewrite = prdcfg.get('rewrite', False)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=csvtimeinfo_path)

        csvfname = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'], ['csv'],
            timeinfo=csvtimeinfo_file, timeformat=timeformat,
            runinfo=prdcfg['runinfo'])[0]

        csvfname = savedir+csvfname

        start_time = dataset[field_name]['timeinfo']

        values = np.ma.masked_all(3)
        values[1] = dataset[field_name]['value']
        if 'samples' in dataset[field_name]:
            values[0] = np.percentile(
                dataset[field_name]['samples'].compressed(), quantiles[0])
            values[2] = np.percentile(
                dataset[field_name]['samples'].compressed(), quantiles[1])

        write_monitoring_ts(
            start_time, dataset[field_name]['npoints'], values,
            [quantiles[0], 'value', quantiles[1]],
            prdcfg['voltype'], csvfname)
        print('saved CSV file: '+csvfname)

        date, np_t_vec, cquant_vec, lquant_vec, hquant_vec = (
            read_monitoring_ts(csvfname, sort_by_date=sort_by_date))

        if date is None:
            warn(
                'Unable to plot time series. No valid data')
            return None

        if rewrite:
            val_vec = np.ma.asarray(
                [lquant_vec, cquant_vec, hquant_vec]).T
            write_monitoring_ts(
                date, np_t_vec, val_vec, quantiles, prdcfg['voltype'],
                csvfname, rewrite=True)

        figtimeinfo = None
        titldate = ''
        if bias_type == 'instant':
            figtimeinfo = date[0]
            titldate = date[0].strftime('%Y-%m-%d')
        else:
            titldate = (date[0].strftime('%Y%m%d')+'-' +
                        date[-1].strftime('%Y%m%d'))
            if prdcfg.get('add_date_in_fname', False):
                figtimeinfo = date[0]
                timeformat = '%Y'

        figfname_list = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=figtimeinfo, timeformat=timeformat,
            runinfo=prdcfg['runinfo'])

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        titl = (prdcfg['runinfo']+' Monitoring '+titldate)

        labely = generate_field_name_str(prdcfg['voltype'])

        np_min = prdcfg.get('npoints_min', 0)
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        plot_monitoring_ts(
            date, np_t_vec, cquant_vec, lquant_vec, hquant_vec, field_name,
            figfname_list, ref_value=ref_value, vmin=vmin, vmax=vmax,
            np_min=np_min, labelx='Time UTC', labely=labely, titl=titl)
        print('----- save to '+' '.join(figfname_list))

        # generate alarms if needed
        alarm = prdcfg.get('alarm', False)
        if not alarm:
            return figfname_list

        if 'tol_abs' not in prdcfg:
            warn('unable to send alarm. Missing tolerance on target')
            return None

        if 'tol_trend' not in prdcfg:
            warn('unable to send alarm. Missing tolerance in trend')
            return None

        if 'nevents_min' not in prdcfg:
            warn('unable to send alarm. ' +
                 'Missing minimum number of events to compute trend')
            return None

        if 'sender' not in prdcfg:
            warn('unable to send alarm. Missing email sender')
            return None
        if 'receiver_list' not in prdcfg:
            warn('unable to send alarm. Missing email receivers')
            return None

        tol_abs = prdcfg['tol_abs']
        tol_trend = prdcfg['tol_trend']
        nevents_min = prdcfg['nevents_min']
        sender = prdcfg['sender']
        receiver_list = prdcfg['receiver_list']

        np_last = np_t_vec[-1]
        value_last = cquant_vec[-1]

        if np_last < np_min:
            warn('No valid data on day '+date[-1].strftime('%d-%m-%Y'))
            return None

        # check if absolute value exceeded
        abs_exceeded = False
        if ((value_last > ref_value+tol_abs) or
                (value_last < ref_value-tol_abs)):
            warn('Value '+str(value_last)+' exceeds target '+str(ref_value) +
                 ' +/- '+str(tol_abs))
            abs_exceeded = True

        # compute trend and check if last value exceeds it
        mask = np.ma.getmaskarray(cquant_vec)
        ind = np.where(np.logical_and(
            np.logical_not(mask), np_t_vec >= np_min))[0]
        nvalid = len(ind)
        if nvalid <= nevents_min:
            warn('Not enough points to compute reliable trend')
            np_trend = 0
            value_trend = np.ma.masked
        else:
            np_trend_vec = np_t_vec[ind][-(nevents_min+1):-1]
            data_trend_vec = cquant_vec[ind][-(nevents_min+1):-1]

            np_trend = np.sum(np_trend_vec)
            value_trend = np.sum(data_trend_vec*np_trend_vec)/np_trend

        trend_exceeded = False
        if np_trend > 0:
            if ((value_last > value_trend+tol_trend) or
                    (value_last < value_trend-tol_trend)):
                warn('Value '+str(value_last)+'exceeds trend ' +
                     str(value_trend)+' +/- '+str(tol_trend))
                trend_exceeded = True

        if abs_exceeded is False and trend_exceeded is False:
            return None

        alarm_dir = savedir+'/alarms/'
        if not os.path.isdir(alarm_dir):
            os.makedirs(alarm_dir)
        alarm_fname = make_filename(
            'alarm', prdcfg['dstype'], prdcfg['voltype'], ['txt'],
            timeinfo=start_time, timeformat='%Y%m%d')[0]
        alarm_fname = alarm_dir+alarm_fname

        field_dict = pyart.config.get_metadata(field_name)
        param_name = get_field_name(field_dict, field_name)
        param_name_unit = param_name+' ['+field_dict['units']+']'

        write_alarm_msg(
            prdcfg['RadarName'][0], param_name_unit, start_time, ref_value,
            tol_abs, np_trend, value_trend, tol_trend, nevents_min, np_last,
            value_last, alarm_fname)

        print('----- saved monitoring alarm to '+alarm_fname)

        subject = ('NO REPLY: '+param_name+' monitoring alarm for radar ' +
                   prdcfg['RadarName'][0]+' on day ' +
                   start_time.strftime('%d-%m-%Y'))
        send_msg(sender, receiver_list, subject, alarm_fname)

        return alarm_fname


    if prdcfg['type'] == 'SAVEVOL' or prdcfg['type'] == 'SAVEVOL_VOL':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        file_type = prdcfg.get('file_type', 'nc')
        physical = prdcfg.get('physical', True)
        compression = prdcfg.get('compression', 'gzip')
        compression_opts = prdcfg.get('compression_opts', 6)

        new_dataset = deepcopy(dataset['radar_out'])
        new_dataset.fields = dict()
        new_dataset.add_field(
            field_name, dataset['radar_out'].fields[field_name])

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'savevol', prdcfg['dstype'], prdcfg['voltype'], [file_type],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        if file_type == 'nc':
            pyart.io.write_cfradial(fname, new_dataset, physical=physical)
        elif file_type == 'h5':
            pyart.aux_io.write_odim_h5(
                fname, new_dataset, physical=physical,
                compression=compression, compression_opts=compression_opts)
        else:
            warn('Data could not be saved. ' +
                 'Unknown saving file type '+file_type)
            return None

        print('saved file: '+fname)

        return fname


    if prdcfg['type'] == 'SAVEALL' or prdcfg['type'] == 'SAVEALL_VOL':
        file_type = prdcfg.get('file_type', 'nc')
        datatypes = prdcfg.get('datatypes', None)
        physical = prdcfg.get('physical', True)
        compression = prdcfg.get('compression', 'gzip')
        compression_opts = prdcfg.get('compression_opts', 6)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'savevol', prdcfg['dstype'], 'all_fields', [file_type],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        field_names = None
        if datatypes is not None:
            field_names = []
            for datatype in datatypes:
                field_names.append(get_fieldname_pyart(datatype))

        if file_type == 'nc':
            if field_names is not None:
                radar_aux = deepcopy(dataset['radar_out'])
                radar_aux.fields = dict()
                for field_name in field_names:
                    if field_name not in dataset['radar_out'].fields:
                        warn(field_name+' not in radar object')
                    else:
                        radar_aux.add_field(
                            field_name,
                            dataset['radar_out'].fields[field_name])
            else:
                radar_aux = dataset['radar_out']
            pyart.io.write_cfradial(fname, radar_aux, physical=physical)
        elif file_type == 'h5':
            pyart.aux_io.write_odim_h5(
                fname, dataset['radar_out'], field_names=field_names,
                physical=physical, compression=compression,
                compression_opts=compression_opts)
        else:
            warn('Data could not be saved. ' +
                 'Unknown saving file type '+file_type)

        print('saved file: '+fname)

        return fname

    if prdcfg['type'] == 'SAVESTATE':
        if prdcfg['lastStateFile'] is None:
            warn('Unable to save last state file. File name not specified')
            return None

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        max_time = np.max(dataset['radar_out'].time['data'])
        units = dataset['radar_out'].time['units']
        calendar = dataset['radar_out'].time['calendar']
        last_date = num2date(max_time, units, calendar)

        write_last_state(last_date, prdcfg['lastStateFile'])
        print('saved file: '+prdcfg['lastStateFile'])

        return prdcfg['lastStateFile']

    if prdcfg['type'] == 'SAVE_FIXED_ANGLE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdsavedir, timeinfo=None)

        fname = make_filename(
            'ts', prdcfg['dstype'], 'fixed_angle', ['csv'],
            timeinfo=None, runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        write_fixed_angle(
            prdcfg['timeinfo'], dataset['radar_out'].fixed_angle['data'][0],
            dataset['radar_out'].latitude['data'][0],
            dataset['radar_out'].longitude['data'][0],
            dataset['radar_out'].altitude['data'][0],
            fname)
        print('saved file: '+fname)

        return fname

    warn(' Unsupported product type: ' + prdcfg['type'])
    return None
