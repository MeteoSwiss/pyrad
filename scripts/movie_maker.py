#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
movie_maker
================================================

This program produces a movie out of all files present in a folder

"""

# Author: fvj
# License: BSD 3 clause

from moviepy.editor import *
import glob
import os
import datetime
import atexit

print(__doc__)


def main():
    """
    """
    print("====== Movie maker started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Movie maker finished: ")

    file_type = 'png'
    movie_type = 'avi'
    codec = 'png'
    frames_per_second = 1
    movie_path = '/users/jfigui/movies/'

    file_path_list = [
        '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/2017-07-19/reflectivity_traj/ALT5000_flash_time/',
        '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/2017-07-19/reflectivity_traj/ALT5000_flash_rel_alt/',
        '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/2017-07-30/reflectivity_traj/ALT5000_flash_time/',
        '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/2017-07-30/reflectivity_traj/ALT5000_flash_rel_alt/',
        '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/2017-08-01/reflectivity_traj/ALT5000_flash_time/',
        '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/2017-08-01/reflectivity_traj/ALT5000_flash_rel_alt/']

    movie_name_list = [
        '20170719_allflash_cappi_time_TRAJ_LIGHTNING_dBZc_alt5000.0',
        '20170719_allflash_cappi_rel_alt_TRAJ_LIGHTNING_dBZc_alt5000.0',
        '20170730_allflash_cappi_time_TRAJ_LIGHTNING_dBZc_alt5000.0',
        '20170730_allflash_cappi_rel_alt_TRAJ_LIGHTNING_dBZc_alt5000.0',
        '20170801_allflash_cappi_time_TRAJ_LIGHTNING_dBZc_alt5000.0',
        '20170801_allflash_cappi_rel_alt_TRAJ_LIGHTNING_dBZc_alt5000.0']

    for i, file_path in enumerate(file_path_list):
        movie_name = movie_name_list[i]
        file_list = sorted(glob.glob(file_path+'*.'+file_type))

        print(file_list)

        # Generate clip
        clip = ImageSequenceClip(file_list, fps=frames_per_second)
        # Write out clip
        if not os.path.isdir(movie_path):
            os.makedirs(movie_path)
        clip.write_videofile(movie_path+movie_name+'.'+movie_type, codec=codec)

        print('Created movie '+movie_path+movie_name)


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
