#!/home/daniel/anaconda3/bin/python
# -*- coding: utf-8 -*-

"""
================================================
rewrite_monitoring
================================================

This program rewrites a intercomp time series files into the correct
time order

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import os

from pyrad.io import read_intercomp_scores_ts, write_intercomp_scores_ts
from pyrad.graph import plot_intercomp_scores_ts
from pyrad.io import get_fieldname_pyart


print(__doc__)


def main():
    """
    """

    input_base = (
        '/store/msrad/radar/pyrad_products/rad4alp_intercomp/')
    output_base = (
        '/store/msrad/radar/pyrad_products/rad4alp_intercomp/')
    rad1_vec = [
        'A', 'A', 'A', 'A', 'D', 'D', 'D', 'D', 'L', 'L', 'L', 'L', 'P', 'P',
        'P', 'P', 'W', 'W', 'W', 'W']
    rad2_vec = [
        'D', 'L', 'P', 'W', 'A', 'L', 'P', 'W', 'A', 'D', 'P', 'W', 'A', 'D',
        'L', 'W', 'A', 'D', 'L', 'P']
    var_vec = ['dBZ', 'dBZv']
    year_vec = [datetime.datetime(2018, 1, 1)]

    plot_data = True
    np_min = 1000
    corr_min = 0.7

    print("====== Monitoring rewriting started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Monitoring rewriting finished: ")

    for i, rad1 in enumerate(rad1_vec):
        rad2 = rad2_vec[i]

        print('Processing Radars '+rad1+'-'+rad2)
        for j, var in enumerate(var_vec):

            input_path = (
                input_base+'PL'+rad1+'_PL'+rad2+'_'+var+'_avg_intercomp/' +
                'PL'+rad1+'_PL'+rad2+'_INTERCOMP_TS/')
            output_path = (
                output_base+'PL'+rad1+'_PL'+rad2+'_'+var+'_avg_intercomp/' +
                'PL'+rad1+'_PL'+rad2+'_INTERCOMP_TS/')
            if not os.path.isdir(output_path):
                os.makedirs(output_path)

            print('- Processing Variable '+var)
            for k, year in enumerate(year_vec):
                print('-- Processing Year '+year.strftime('%Y'))
                fname_input = (
                    input_path+year.strftime('%Y')+'_ts_INTERCOMP_TIME_AVG_' +
                    var+'_'+rad1+'-'+rad2+'.csv')
                fname_output = (
                    output_path+year.strftime('%Y')+'_ts_INTERCOMP_TIME_AVG_' +
                    var+'_'+rad1+'-'+rad2+'.csv')
                figfname = [
                    output_path+year.strftime('%Y')+'_ts_INTERCOMP_TIME_AVG_' +
                    var+'_'+rad1+'-'+rad2+'.png']

                (date_vec, np_vec, meanbias_vec, medianbias_vec,
                 quant25bias_vec, quant75bias_vec, modebias_vec, corr_vec,
                 slope_vec, intercep_vec, intercep_slope1_vec) = (
                    read_intercomp_scores_ts(fname_input, sort_by_date=True))

                if date_vec is None:
                    continue

                stats = {
                    'npoints': np_vec,
                    'meanbias': meanbias_vec,
                    'medianbias': medianbias_vec,
                    'quant25bias': quant25bias_vec,
                    'quant75bias': quant75bias_vec,
                    'modebias': modebias_vec,
                    'corr': corr_vec,
                    'slope': slope_vec,
                    'intercep': intercep_vec,
                    'intercep_slope_1': intercep_slope1_vec
                }
                field_name = get_fieldname_pyart(var)
                fname = write_intercomp_scores_ts(
                    date_vec, stats, field_name, fname_output,
                    rad1_name=rad1, rad2_name=rad2, rewrite=True)

                print('written file '+fname)

                if not plot_data:
                    continue

                titldate = (date_vec[0].strftime('%Y%m%d')+'-' +
                            date_vec[-1].strftime('%Y%m%d'))
                titl = (
                    rad1+'-'+rad2+' '+field_name+' intercomparison '+titldate)

                fname = plot_intercomp_scores_ts(
                    date_vec, np_vec, meanbias_vec, medianbias_vec,
                    quant25bias_vec, quant75bias_vec, modebias_vec, corr_vec,
                    slope_vec, intercep_vec, intercep_slope1_vec, figfname,
                    ref_value=0., np_min=np_min, corr_min=corr_min,
                    labelx='Time UTC', titl=titl)
                print('plotted file '+' '.join(fname))


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
