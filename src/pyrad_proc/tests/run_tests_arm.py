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
    
