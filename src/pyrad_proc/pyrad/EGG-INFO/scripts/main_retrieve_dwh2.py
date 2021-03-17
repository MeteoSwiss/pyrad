#!/home/daniel/anaconda3/bin/python
# -*- coding: utf-8 -*-

"""
================================================
main_retrieve_dwh
================================================

This program puts together stored SwissMetNet data in a single file

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
from warnings import warn

import pandas as pd

print(__doc__)


def main():
    """
    """

    file_path_in = '/data/FLORAKO/v2/'
    file_path_out = '/data/FLORAKO/'    
    avg_time = 3600   
    fill_value = -99999.

    smn_station_vec = ['PIL', 'WFJ', 'MTR', 'GUE', 'DIA']
    tstart_vec = ['2007', '2008', '2009', '2010', '2011', '2012', '2013',
                  '2014', '2015', '2016']
    tend_vec = ['2008', '2009', '2010', '2011', '2012', '2013', '2014',
                '2015', '2016', '2017']
    param_vec = ['tre200s0', 'tde200s0', 'ure200s0', 'gor000z0', 'fkl010z0']

    print("====== retrieval from DWH started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== retrieval from DWH finished: ")

    for station in smn_station_vec:        
        df_out = None
        for param in param_vec:
            df_param = None
            for time in range(len(tstart_vec)):
                tstart = tstart_vec[time]+'0101000000'
                tend = tend_vec[time]+'0101000000'                        
                
                print('\n--- Processing station '+station+' parameter ' +
                      param+' between '+tstart+' and '+tend)

                base_name = station+'_'+param+'_'+tstart+'-'+tend
                file_name = file_path_in+base_name+'_avg'+str(avg_time)+'s.csv'                
                
                try:
                    df = pd.read_csv(file_name, parse_dates=['datetime'], index_col=0, names=['datetime', param+'_avg', param+'_std'], header=0)
                
                    if df_param is None:
                        df_param = df
                    else:
                        df_param = pd.concat([df_param, df])
                except OSError:
                    warn('Unable to read file '+file_name)        
        
            if df_out is None:
                df_out = df_param
            else:
                df_out = pd.concat([df_out, df_param], join='outer', axis=1)    
        
        if fill_value is not None:
            df_out = df_out.fillna(value=fill_value)
        
        file_name = file_path_out+station+'_'+tstart_vec[0]+'0101000000'+'-'+tend_vec[-1]+'0101000000'+'_avg'+str(avg_time)+'s.csv'
        df_out.to_csv(file_name)     
                
                


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
