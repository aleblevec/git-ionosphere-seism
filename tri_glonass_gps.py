#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 17:13:37 2019

@author: antoineleblevec
"""

import numpy as np
import os.path
import os 
from os import scandir
import matplotlib.pyplot as plt
import csv
from os import chdir 
from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap
import fonction_1 as f1 
from scipy.signal import savgol_filter
import shutil

year = str('2016')
day = str('318')
station = 'mqzg'
name_dir = os.path.join('/Users/antoineleblevec/Desktop/tec-suite-master_2/tec', year, day, station)

rep = os.path.abspath(os.path.expanduser(name_dir))
files = os.listdir(rep)
try:
    os.mkdir('Glonass')
except OSError:
    pass
#for i in range(0,len(files)) :
#    name_data = files [i]
#    x = files[i].split("_")
#    sat = x [1]

#    if sat[0] == "R" :
#        shutil.copy("./{0}".format(name_data), './Glonass')
#        
#    if sat [0] == "G" :
#        shutil.copy("./{0}".format(name_data), './GPS')

    