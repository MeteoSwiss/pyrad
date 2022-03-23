#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
movie_maker
================================================

This program produces a movie and a gif out of all files present in a folder

"""

# Author: fvj
# License: BSD 3 clause

import glob
import os
import datetime
import atexit

from moviepy.editor import ImageSequenceClip

print(__doc__)


def main():
    """
    main programme
    """
    print("====== Movie maker started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Movie maker finished: ")

    file_type = 'png'
    movie_type = 'mp4'
    codec = 'mpeg4'
    frames_per_second = 1
    movie_path = '/utemp/mdso/figuerasiventuraj/movies/'
    file_path_list = [        
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/dBuZ/PPI_EL00/',
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/dBuZ/PPI_EL01/',
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/dBuZ/PPI_EL02/',
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/dBuZ/PPI_EL03/',
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/dBZ/PPI_EL00/',
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/dBZ/PPI_EL01/',
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/dBZ/PPI_EL02/',
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/dBZ/PPI_EL03/',
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/hydroMF_oper/PPI_EL00/',
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/hydroMF_oper/PPI_EL01/',
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/hydroMF_oper/PPI_EL02/',
        '/utemp/mdso/figuerasiventuraj/pyrad_products/MF_ODIM_OPOU_HAIL/2020-07-01/hydroMF_oper/PPI_EL03/']
    movie_name_list = [
        '20200701_OPOU_ppi_RAW_dBuZ_el0.6',
        '20200701_OPOU_ppi_RAW_dBuZ_el1.0',
        '20200701_OPOU_ppi_RAW_dBuZ_el1.4',
        '20200701_OPOU_ppi_RAW_dBuZ_el1.8',
        '20200701_OPOU_ppi_RAW_dBZ_el0.6',
        '20200701_OPOU_ppi_RAW_dBZ_el1.0',
        '20200701_OPOU_ppi_RAW_dBZ_el1.4',
        '20200701_OPOU_ppi_RAW_dBZ_el1.8',
        '20200701_OPOU_ppi_RAW_hydroMF_oper_el0.6',
        '20200701_OPOU_ppi_RAW_hydroMF_oper_el1.0',
        '20200701_OPOU_ppi_RAW_hydroMF_oper_el1.4',
        '20200701_OPOU_ppi_RAW_hydroMF_oper_el1.8']

    create_movie(
        file_path_list, movie_name_list, movie_path, file_type=file_type,
        fps=frames_per_second, movie_type=movie_type, codec=codec)


def create_movie(file_path_list, movie_name_list, movie_path, file_type='png',
                 fps=1, movie_type='mp4', codec='mpeg4'):
    """
    creates the movie.

    can support any type supported by ffmpeg
    some examples:
    movie type / codec
    .avi / rawvideo, png
    .mp4 / libx264, mpeg4
    avi/rawvideo supported by libreoffice
    mp4 supported by windows media player

    Parameters
    ----------
    file_path_list : list of str
        List of folders where to find the images for the movies
    movie_name_list : list of str
        List of movies to create_movie
    movie_path : str
        path where to store the movies
    file_type : str
        the individual images file type
    fps : int
        the frames per second
    movie_type : str
        the type of movie file
    codec : str
        the codec used for the movie

    Returns
    -------
    Nothing

    """
    for movie_name, file_path in zip(movie_name_list, file_path_list):
        file_list = sorted(glob.glob(file_path+'*.'+file_type))
        print(file_list)

        # Generate clip
        clip = ImageSequenceClip(file_list, fps=fps)
        # Write out clip
        if not os.path.isdir(movie_path):
            os.makedirs(movie_path)
        clip.write_videofile(movie_path+movie_name+'.'+movie_type, codec=codec)
        clip.write_gif(movie_path+movie_name+'.gif')

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
