#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 11:30:55 2021

@author: wolfensb
"""
import os
import datetime
from pathlib import Path
 
from pyrad.flow import main

def test_rad4alp_hydro():
    cfgfile = Path(os.environ['PYRAD_EXAMPLES_PATH'], 'config', 'processing',
                   'rad4alp_hydro_PLL.txt')
    os.environ['PYART_CONFIG'] = str(Path(os.environ['PYRAD_EXAMPLES_PATH'],
                                     'mch_config.py'))
    starttime = datetime.datetime(2022,6,28, 7, 20)
    endtime = datetime.datetime(2022,6,28, 7, 25)
    main(str(cfgfile), starttime=starttime, endtime=endtime)

def test_cpc():
    cfgfile = Path(os.environ['PYRAD_EXAMPLES_PATH'], 'config', 'processing',
                   'cpc_gif.txt')
    starttime = datetime.datetime(2017,8,1,17,0)
    endtime = datetime.datetime(2017,8,1,17,10)
    main(str(cfgfile), starttime=starttime, endtime=endtime)

def test_rainbow_vol():
    cfgfile = Path(os.environ['PYRAD_EXAMPLES_PATH'], 'config', 'processing',
                   'rainbow_vol.txt')
    starttime = datetime.datetime(2019,10,21, 8, 20)
    endtime = datetime.datetime(2019,10,21, 8, 25)
    main(str(cfgfile), starttime=starttime, endtime=endtime)
  
def test_rainbow_rhi():
    cfgfile = Path(os.environ['PYRAD_EXAMPLES_PATH'], 'config', 'processing',
                   'rainbow_rhi.txt')
    starttime = datetime.datetime(2018,5,2, 2, 25)
    endtime = datetime.datetime(2018,5,2, 2, 30)
    main(str(cfgfile), starttime=starttime, endtime=endtime)
    
def test_rainbow_poi():
    cfgfile = Path(os.environ['PYRAD_EXAMPLES_PATH'], 'config', 'processing',
                   'rainbow_poi.txt')
    starttime = datetime.datetime(2019,10,15, 9, 35)
    endtime = datetime.datetime(2019,10,15, 9, 40)
    main(str(cfgfile), starttime=starttime, endtime=endtime)

def test_nexrad2():
    cfgfile = Path(os.environ['PYRAD_EXAMPLES_PATH'], 'config', 'processing',
                   'nexrad2.txt')
    starttime = datetime.datetime(2013, 7, 17, 19, 50)
    endtime = datetime.datetime(2013, 7, 17, 19, 55)
    main(str(cfgfile), starttime=starttime, endtime=endtime)

def test_mf_pag_mtcy():
    cfgfile = Path(os.environ['PYRAD_EXAMPLES_PATH'], 'config', 'processing',
                   'mf_pag_mtcy.txt')
    starttime = datetime.datetime(2021, 1, 14, 23, 45)
    endtime = datetime.datetime(2021, 1, 14, 23, 50)
    main(str(cfgfile), starttime=starttime, endtime=endtime)
    
def test_mf_pag_mtcy():
    cfgfile = Path(os.environ['PYRAD_EXAMPLES_PATH'], 'config', 'processing',
                   'mf_odim_mtcy.txt')
    starttime = datetime.datetime(2021, 1, 14, 0, 0)
    endtime = datetime.datetime(2021, 1, 14, 0, 5)
    main(str(cfgfile), starttime=starttime, endtime=endtime)
    