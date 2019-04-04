#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 17:20:04 2019

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
from labellines import labelLine, labelLines

# =============================================================================
# Constantes
# =============================================================================
Re = 6371032 
H = 350e3
lllat = 31.7; urlat = 43.7; lllon = 130.6; urlon = 145.6
time_inf = 11
time_sup = 11.6 
epoque_inf = 1320
epoque_sup = 1393
min_ele = 10 
window = 27
order = 4

name_dir = os.path.join('/Users/antoineleblevec/Desktop/G20')
print (name_dir)
rep = os.path.abspath(os.path.expanduser(name_dir))
files = os.listdir(rep) 
files.sort()
print (files)

lat = []
lon = []
reduce_tec = []
#lattotal = []
#lontotal = []
#reduce_tectotal = []
#latt = []
#lont = []
#reduce_tect = []

for j in range(0,len(files)) : 
    name_file = files [j] #nom du fichier lu
    data = np.loadtxt(rep + '/' + name_file) #chargement des données
    rlon, rlat = f1.lecture_lat_lon_sat(rep, name_file) #calcul longitude et latitude station lue
    epoque = data[:,0] # stockage des époques dans "epoque"
    ep  = epoque [(epoque_inf < epoque) & (epoque < epoque_sup)] #stockage des époques intéressantes dans ep
    ele = data [:,2]    
    elevation = ele [(epoque_inf < epoque) & (epoque < epoque_sup)]
    az = data [:,3]
    az = az [(epoque_inf < epoque) & (epoque < epoque_sup)]
    az = np.radians(az)
    elevation = np.radians(elevation)
    x = np.arcsin((Re*np.cos(elevation))/(Re + H))
    ksi = np.pi/2 - elevation - x
    lat = np.arcsin(np.sin(rlat)*np.cos(ksi) + np.cos(rlat)*np.sin(ksi)*np.cos(az))
    lon = rlon + np.arcsin(np.sin(ksi)*np.sin(az)/np.cos(lat))
    b = data [:,4]
    e = b[(epoque_inf < epoque) & (epoque < epoque_sup)]
    tec = e - min(e)
    vtec = tec * np.cos(x)
    vtec_polyn = savgol_filter(vtec, window, order)
    reduce_tec = vtec - vtec_polyn 
#    print (type(lat))
print (lon[3], lat[3], reduce_tec[3])

    
    
    
    
    
    
    
    
    
    
    
    
    