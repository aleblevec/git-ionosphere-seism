#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 10:44:15 2019

@author: antoineleblevec
"""

# =============================================================================
# Modules 
# =============================================================================
import numpy as np
import os.path
from os import scandir
import matplotlib.pyplot as plt
import csv
from os import chdir 
from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap
from scipy.signal import butter, lfilter

# =============================================================================
# Lecture du fichier 
# =============================================================================
datat = np.loadtxt('/Users/antoineleblevec/Desktop/2016_seism/coord_stations.txt', usecols = (1,2,3,4,5))
data = open('/Users/antoineleblevec/Desktop/2016_seism/coord_stations.txt', "r")

# =============================================================================
# Stockage des noms des stations dans : station 
# =============================================================================
contenu = data.readlines()
station = []
for j in range (36) : 
    list = contenu[j].split()
    station.append(list [0])

# =============================================================================
# Stockage des latitudes et longitudes des stations dans lon et lat
# =============================================================================
lon = []
lat = []
for i in datat : 
    lon.append(i[4])
    lat.append(i[3])

# =============================================================================
# Epicentre longitude latitude 
# =============================================================================
elon = 173.054
elat = -42.737

# =============================================================================
# Plot du fond de carte
# =============================================================================
m = Basemap(width=1500000,height=1500000,projection='lcc',
                resolution='c',lat_1=-80.,lat_2=-20,lat_0=-40,lon_0=176.)
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.bluemarble()

# =============================================================================
# Plot des stations sur le fond de carte  
# =============================================================================
for j in range (len(lon)) : 
    xpt,ypt = m(lon[j],lat[j])
    m.plot(xpt,ypt,'wo', ms = 1) 
    plt.text(xpt+10000,ypt+10000,'%s' %station[j], color = 'white', fontsize = 5)
    expt, eypt = m(elon,elat)
    m.plot(expt,eypt,'r*', ms = 3)

plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/map_stations_seism2016.png')
plt.show()


    
    


