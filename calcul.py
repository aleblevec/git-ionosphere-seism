#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 10:44:11 2019

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
import matplotlib.animation as animation


#Re = 6371032 
#H = 350e3
#
#elevation = (80.65114 * np.pi) / 180 
#az = (217.74742 * np.pi) / 180
#
#rlat = (-41.06256873321562 * np.pi) / 180
#rlon = (174.34658955 * np.pi) / 180
#print (rlon)
#
#x = np.arcsin((Re*np.cos(elevation))/(Re + H))
#
#ksi = np.pi/2 - elevation - x
#
#lat = np.arcsin(np.sin(rlat)*np.cos(ksi) + np.cos(rlat)*np.sin(ksi)*np.cos(az))
#lon = rlon + np.arcsin(np.sin(ksi)*np.sin(az)/np.cos(lat))
##print (lon)

for a in range(11) : 
    print (2a)