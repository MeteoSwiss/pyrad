#!/home/daniel/anaconda3/bin/python
# -*- coding: utf-8 -*-

"""
================================================
main_retrieve_dwh
================================================

This program retrieves parameters from selected SwissMetNet stations, stores
them in a file and plots them and computes the average and standard deviation
over a period and stores it and plots it.

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
from subprocess import run

from pyrad.io import read_smn2, write_smn
from pyrad.graph import plot_timeseries
from pyrad.util import time_series_statistics, join_time_series

print(__doc__)


def main():
    """
    """

    file_path = '/data/FLORAKO/fm_transmitters_measurements/10min/'
    img_ext = 'png'
    avg_time = 3600
    base_time = 0
    
    # tstart+t_init and tend+t_final make the start and stop times for the 
    # retrieval of the data
    # t_init = '0101000000'
    # t_final = '0101000000'
    # 
    # tstart_vec = ['2007', '2008', '2009', '2010', '2011', '2012', '2013',
    #               '2014', '2015', '2016']
    # tend_vec = ['2008', '2009', '2010', '2011', '2012', '2013', '2014',
    #             '2015', '2016', '2017']

    t_init = '160000'
    t_final = '100000'
    
    tstart_vec = ['20140421']
    tend_vec = ['20160425']

    # SMN stations to retrieve
    smn_station_vec = ['BER']
    #    smn_station_vec = ['PIL', 'WFJ', 'MTR', 'GUE', 'DIA', 'VSTSN']
    
    # meteorological parameters
    # tre200s0: air temperature  deg C at 2 m [10 min resolution]
    # tde200s0: dew point temperature deg C at 2 m [10 min resolution]
    # ure200s0: relative humidity % at 2 m [10 min resolution]
    # gor000z0: global radiation W/m2[10 min average]
    # fkl010z0: mean wind speed in 10 min m/s
    # fkl010z1: max wind speed in 10 min m/s
    # prestas0: air pressure at station height [hPa]
    # rre150z0: precipitation 10 min accumulation
    param_vec = ['tre200s0', 'gor000z0', 'ure200s0', 'prestas0', 'fkl010z0', 'fkl010z1', 'rre150z0'] 
    ymin_vec = [-40., -40., 0., 600., 0., 0., 0., 0.]
    ymax_vec = [30., 1500., 110., 1300., 20., 50., 30.]

    
 
#    param_vec = ['tre200s0', 'tde200s0', 'ure200s0', 'gor000z0', 'fkl010z0', 'prestas0']
#    ymin_vec = [-40., -40., 0., -40., 0., 600.]
#    ymax_vec = [30., 30., 110., 1500., 20., 900.]

#    param_vec = ['fkl010z1']
#    ymin_vec = [0.]
#    ymax_vec = [50.]

    print("====== retrieval from DWH started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== retrieval from DWH finished: ")

    for station in smn_station_vec:
        i = 0
        for param in param_vec:
            ymin = ymin_vec[i]
            ymax = ymax_vec[i]
            i += 1
            for time in range(len(tstart_vec)):
                tstart = tstart_vec[time]+t_init
                tend = tend_vec[time]+t_final

                print('\n--- Processing station '+station+' parameter ' +
                      param+' between '+tstart+' and '+tend)

                base_name = station+'_'+param+'_'+tstart+'-'+tend
                file_name = file_path+base_name+'.csv'
                with open(file_name, 'wb') as out_file:
                    result = run(
                        ['jretrievedwh.sh',  '-d', 'COMMA', '-s', 'surface',
                         '-i', 'nat_abbr,'+station, '-t', tstart+','+tend,
                         '-n', param], stdout=out_file)
                id, date, value = read_smn2(file_name)

                if date is None:
                    continue

                plot_timeseries(
                    date, [value], [file_path+base_name+'.'+img_ext],
                    labelx='Time [UTC]', labely=param, labels=[station],
                    title=station+' '+tstart+'-'+tend, period=0,
                    timeformat=None, colors=None, linestyles=None, ymin=ymin,
                    ymax=ymax)

                date_avg, value_avg = time_series_statistics(
                    date, value, avg_time=avg_time, base_time=base_time,
                    method='mean', dropnan=True)
                date_std, value_std = time_series_statistics(
                    date, value, avg_time=avg_time, base_time=base_time,
                    method='std', dropnan=True)

                date_series, value_avg, value_std = join_time_series(
                    date_avg, value_avg, date_std, value_std, dropnan=True)

                plot_timeseries(
                    date_series,
                    [value_avg, value_avg+value_std, value_avg-value_std],
                    [file_path+base_name+'_avg'+str(avg_time)+'s.'+img_ext],
                    labelx='Time [UTC]', labely=param, labels=None,
                    title=station+' '+tstart+'-'+tend+' avg '+str(
                        avg_time)+' s', period=0, timeformat=None,
                    colors=['b', 'r', 'r'], linestyles=['-', '--', '--'],
                    ymin=ymin, ymax=ymax)

                write_smn(
                    date_series, value_avg, value_std,
                    file_path+base_name+'_avg'+str(avg_time)+'s.csv')


def _print_end_msg(text):
    """
    prints end message

    Parameters
    ----------
    text : str
        the text to be printed

    Returns
    -------
    Nothing

    """
    print(text + datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))


# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
if __name__ == "__main__":
    main()
