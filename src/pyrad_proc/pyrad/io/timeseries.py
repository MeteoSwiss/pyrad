"""
pyrad.io.timeseries
===================

TimeSeries class implementation for holding timeseries data.

.. autosummary::
    :toctree: generated/

    TimeSeries

"""

from datetime import datetime

import numpy as np

from ..graph.plots import plot_histogram
from ..graph.plots_timeseries import plot_timeseries
from ..util.radar_utils import compute_histogram
from ..io.io_aux import get_fieldname_pyart


class TimeSeries(object):
    """
    Holding timeseries data and metadata.

    Attributes
    ----------
    description : array of str
        Description of the data of the time series.
    time_vector : array of datetime objects
    timeformat : how to print the time (default:
                 'Date, UTC [seconds since midnight]'
    dataseries : List of _dataSeries object holding the
                 data

    Methods:
    --------
    add_dataseries : Add a data series to the object
    write : Write time series to a file
    plot : Plot a figure of a time series

    """

    def __init__(
        self, desc, timevec=None, timeformat=None, maxlength=None, datatype=""
    ):
        """
        Initalize the object.

        Parameters
        ----------
        desc : array of str
        timevec : array of datetime
        timeformat : specifies time format
        maxlength : Maximal length of the time series
        num_el : Number of values in the time series
        """

        self.description = desc
        if timevec is None:
            if maxlength is None:
                raise Exception(
                    "ERROR: Either 'timevec' or 'maxlength'" " must be defined"
                )
            self.maxlength = maxlength
            self.time_vector = np.empty(maxlength, dtype=datetime)
            self.num_el = 0
        else:
            self.time_vector = timevec
            self.maxlength = timevec.size
            self.num_el = timevec.size
        self.timeformat = timeformat
        self.dataseries = []
        self.datatype = datatype

    def add_dataseries(
        self,
        label,
        unit_name,
        unit,
        dataseries=None,
        plot=True,
        color=None,
        linestyle=None,
    ):
        """
        Add a new data series to the timeseries object.
        The length of the data vector must be the same as the
        length of the time vector.
        """

        if dataseries is not None:
            if len(dataseries) != self.num_el:
                raise Exception(
                    "ERROR: Number of data series sample do "
                    "not correspond to time vector ('%s')" % label
                )
        else:
            dataseries = np.ma.empty(self.maxlength)

        ds = _DataSeries(
            label,
            unit_name,
            unit,
            dataseries,
            plot=plot,
            color=color,
            linestyle=linestyle,
        )
        self.dataseries.append(ds)

    def add_timesample(self, dt, values):
        """
        Add a new sample to the time series.
        """

        if self.num_el + 1 > self.maxlength:  # jgr changed from >=
            raise Exception(
                "ERROR: Cannot add time series sample. Max" " length reached."
            )

        self.time_vector[self.num_el] = dt
        for val, ds in zip(values, self.dataseries):
            ds.set_value(self.num_el, val)
        self.num_el += 1

    def write(self, fname):
        """
        Write time series output
        """

        print("----- write to '%s'" % fname)

        try:
            tsfile = open(fname, "w")
        except Exception:
            raise Exception("ERROR: Could not create file '%s'" % fname)

        print("# Weather radar timeseries data file", file=tsfile)
        print("# Project: MALSplus", file=tsfile)
        print(
            "# Start : %s UTC" % self.time_vector[0].strftime("%Y-%m-%d %H:%M:%S"),
            file=tsfile,
        )
        print(
            "# End   : %s UTC"
            % self.time_vector[self.num_el - 1].strftime("%Y-%m-%d %H:%M:%S"),
            file=tsfile,
        )
        print("# Header lines with comments are preceded by '#'", file=tsfile)
        for line in self.description:
            print("# %s" % line, file=tsfile)
        print("#", file=tsfile)

        # Make raw header
        if self.timeformat is None:
            print("# Date, UTC [seconds since midnight]", end="", file=tsfile)
        else:
            print("# Date [%s]" % self.timeformat, end="", file=tsfile)

        for ds in self.dataseries:
            print(", %s [%s]" % (ds.label, ds.unit), end="", file=tsfile)
        print("", file=tsfile)

        # Store the data
        for i in range(self.num_el):
            if self.timeformat is None:
                dt = self.time_vector[i]
                daystr = dt.strftime("%d-%b-%Y")
                secs = (
                    dt.hour * 3600.0
                    + dt.minute * 60.0
                    + dt.second
                    + dt.microsecond / 1000000.0
                )
                print("%s, %14.4f" % (daystr, secs), end="", file=tsfile)
            else:
                print(
                    self.time_vector[i].strftime(self.timeformat), end="", file=tsfile
                )

            for ds in self.dataseries:
                print(
                    ", %14.4f" % (np.ma.filled(ds.data, fill_value=np.nan)[i]),
                    end="",
                    file=tsfile,
                )

            print("", file=tsfile)

        tsfile.close()

    def plot(self, fname, ymin=None, ymax=None):
        """
        Make a figure of a time series
        """

        found = False
        labely = None
        ds_list = []
        color_list = []
        lstyle_list = []
        for ds in self.dataseries:
            if ds.plot:
                found = True
                ds_list.append(ds.data[: self.num_el])
                color_list.append(ds.color)
                lstyle_list.append(ds.linestyle)
                if labely is None:
                    labely = "%s [%s]" % (ds.unit_name, ds.unit)

        if not found:
            raise Exception("ERROR: Undefined time series '%s'" % ds.label)

        print("----- plot to '%s'" % fname)

        title = "Trajectory Time Series %s" % self.time_vector[0].strftime("%Y-%m-%d")
        plot_timeseries(
            self.time_vector[: self.num_el],
            ds_list,
            [fname],
            title=title,
            labels=None,
            labely=labely,
            timeformat="%H:%M",
            colors=color_list,
            linestyles=lstyle_list,
            ymin=ymin,
            ymax=ymax,
        )

    def plot_hist(self, fname, step=None):
        """
        Make histograms of time series
        """
        for ds in self.dataseries:
            if ds.plot:
                bins, values = compute_histogram(
                    ds.data[: self.num_el],
                    get_fieldname_pyart(self.datatype),
                    step=step,
                )
                fname2 = fname.replace(".", "_" + ds.label + ".")
                plot_histogram(
                    bins,
                    values,
                    [fname2],
                    labelx="%s [%s]" % (ds.unit_name, ds.unit),
                    titl=(
                        "Trajectory Histogram %s"
                        % self.time_vector[0].strftime("%Y-%m-%d")
                    ),
                )
                print("----- plot to '%s'" % fname2)


class _DataSeries(object):
    """
    Hold a data vector and some meta information.
    """

    def __init__(
        self, label, unit_name, unit, data, plot=True, color=None, linestyle=None
    ):
        """
        Initalize the object.
        """

        self.label = label
        self.unit_name = unit_name
        self.unit = unit
        self.data = data
        self.plot = plot
        self.color = color
        self.linestyle = linestyle

    def set_value(self, i, val):
        """
        Append value to array
        """
        self.data[i] = val
