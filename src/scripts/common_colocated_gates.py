#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
common_colocated_gates
================================================

This program reads colocated gates files from two radars
and creates a new file with the gates that are common to
both radars

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import numpy as np
import pandas as pd
from copy import deepcopy

from pyrad.io import read_colocated_gates, write_colocated_gates

print(__doc__)


def main():
    """
    """

    file_path = (
        '/srn/analysis/pyrad_products/rad4alp_intercomp/colocated_gates/')
    rad1_vec = ['A', 'A', 'A', 'A', 'D', 'D', 'D', 'L', 'L', 'P']
    rad2_vec = ['D', 'L', 'P', 'W', 'L', 'P', 'W', 'P', 'W', 'W']

    print("====== common colocated gates started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== common colocated gates finished: ")

    for i, rad1 in enumerate(rad1_vec):
        rad2 = rad2_vec[i]

        print('Radars: '+rad1+' '+rad2)
        fname1 = (file_path+'PL'+rad1+'_'+'PL'+rad2 +
                  '/info_COLOCATED_GATES_PL'+rad1+'_PL'+rad2+'.csv')
        fname2 = (file_path+'PL'+rad2+'_'+'PL'+rad1 +
                  '/info_COLOCATED_GATES_PL'+rad2+'_PL'+rad1+'.csv')

        (rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
         rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng) = (
            read_colocated_gates(fname1))

        print('Number of gates rad1-rad2 ', np.shape(rad1_ray_ind))

        (rad2_ray_ind_aux, rad2_rng_ind_aux, rad2_ele_aux, rad2_azi_aux,
         rad2_rng_aux, rad1_ray_ind_aux, rad1_rng_ind_aux, rad1_ele_aux,
         rad1_azi_aux, rad1_rng_aux) = read_colocated_gates(fname2)

        print('Number of gates rad2-rad1 ', np.shape(rad2_ray_ind_aux))

        # make a pool of data
        rad1_ray_ind = np.ma.concatenate((rad1_ray_ind, rad1_ray_ind_aux))
        rad1_rng_ind = np.ma.concatenate((rad1_rng_ind, rad1_rng_ind_aux))
        rad1_ele = np.ma.concatenate((rad1_ele, rad1_ele_aux))
        rad1_azi = np.ma.concatenate((rad1_azi, rad1_azi_aux))
        rad1_rng = np.ma.concatenate((rad1_rng, rad1_rng_aux))
        rad2_ray_ind = np.ma.concatenate((rad2_ray_ind, rad2_ray_ind_aux))
        rad2_rng_ind = np.ma.concatenate((rad2_rng_ind, rad2_rng_ind_aux))
        rad2_ele = np.ma.concatenate((rad2_ele, rad2_ele_aux))
        rad2_azi = np.ma.concatenate((rad2_azi, rad2_azi_aux))
        rad2_rng = np.ma.concatenate((rad2_rng, rad2_rng_aux))

        print('Total number of gates ', np.shape(rad1_ray_ind))

        # create dictionary and put it in pandas framework
        coloc_dict = {
            'rad1_ray_ind': rad1_ray_ind,
            'rad1_rng_ind': rad1_rng_ind,
            'rad1_ele': rad1_ele,
            'rad1_azi': rad1_azi,
            'rad1_rng': rad1_rng,
            'rad2_ray_ind': rad2_ray_ind,
            'rad2_rng_ind': rad2_rng_ind,
            'rad2_ele': rad2_ele,
            'rad2_azi': rad2_azi,
            'rad2_rng': rad2_rng}
        df = pd.DataFrame(data=coloc_dict)

        # keep only duplicated data
        df_common = df[df.duplicated(keep=False)].drop_duplicates()
        common_dict = df_common.to_dict(orient='list')

        print('Number of common gates', df_common.shape)
        print('rad1 elev min/max', np.min(common_dict['rad1_ele']),
              np.max(common_dict['rad1_ele']))
        print('rad2 elev min/max', np.min(common_dict['rad2_ele']),
              np.max(common_dict['rad2_ele']))

        # write resultant output
        fname1_out = (
            file_path+'PL'+rad1+'_'+'PL'+rad2 +
            '/info_common_COLOCATED_GATES_PL'+rad1+'_PL'+rad2+'.csv')
        fname2_out = (
            file_path+'PL'+rad2+'_'+'PL'+rad1 +
            '/info_common_COLOCATED_GATES_PL'+rad2+'_PL'+rad1+'.csv')

        rad1_dict = {
            'rad1_ray_ind': np.asarray(common_dict['rad1_ray_ind']),
            'rad1_rng_ind': np.asarray(common_dict['rad1_rng_ind']),
            'rad1_ele': np.asarray(common_dict['rad1_ele']),
            'rad1_azi': np.asarray(common_dict['rad1_azi']),
            'rad1_rng': np.asarray(common_dict['rad1_rng']),
            'rad2_ray_ind': np.asarray(common_dict['rad2_ray_ind']),
            'rad2_rng_ind': np.asarray(common_dict['rad2_rng_ind']),
            'rad2_ele': np.asarray(common_dict['rad2_ele']),
            'rad2_azi': np.asarray(common_dict['rad2_azi']),
            'rad2_rng': np.asarray(common_dict['rad2_rng'])}

        rad2_dict = {
            'rad1_ray_ind': np.asarray(common_dict['rad2_ray_ind']),
            'rad1_rng_ind': np.asarray(common_dict['rad2_rng_ind']),
            'rad1_ele': np.asarray(common_dict['rad2_ele']),
            'rad1_azi': np.asarray(common_dict['rad2_azi']),
            'rad1_rng': np.asarray(common_dict['rad2_rng']),
            'rad2_ray_ind': np.asarray(common_dict['rad1_ray_ind']),
            'rad2_rng_ind': np.asarray(common_dict['rad1_rng_ind']),
            'rad2_ele': np.asarray(common_dict['rad1_ele']),
            'rad2_azi': np.asarray(common_dict['rad1_azi']),
            'rad2_rng': np.asarray(common_dict['rad1_rng'])}

        write_colocated_gates(rad1_dict, fname1_out)
        write_colocated_gates(rad2_dict, fname2_out)


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
