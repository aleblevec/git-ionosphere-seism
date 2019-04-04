#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 11:44:37 2019

@author: antoineleblevec
"""

# =============================================================================
# Modules
# =============================================================================
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

# =============================================================================
# Definition des constantes
# =============================================================================
Re = 6371032 
H = 350e3
lllat = 31.7; urlat = 43.7; lllon = 130.6; urlon = 145.6
                
# =============================================================================
# Calcul des lon et lat du satellite 
# =============================================================================

year = str('2016')
day = str('318')
station = 'lkta'
type_sat = 'Glonass'
name_dir = os.path.join('/Users/antoineleblevec/Desktop/tec-suite-master_2/tec', year, day, station, type_sat)
print (name_dir)
borne_inf_time = 10.30
borne_sup_time = 12.30 
sueil_elevation = 10 
rep = os.path.abspath(os.path.expanduser(name_dir))
files = os.listdir(rep)
print(files[14])
data = np.loadtxt(name_dir + '/' + files[14])

# calcul de la longitude et lattitude de la station 
rlon, rlat = f1.lecture_lat_lon_sat_anim(name_dir,files)

a = data[:,1] 
# temps d'observation du séisme 
d = a [(int(borne_inf_time) < a) & (a < int(borne_sup_time))] 
# elevation  
ele = data [:,2]
# elevation au moment du séisme    
elevation = ele [(int(borne_inf_time) < a) & (a < int(borne_sup_time))]
az = data [:,3]
az = az [(int(borne_inf_time) < a) & (a < int(borne_sup_time))]
az = np.radians(az)
elevation  = np.radians(elevation)

# Calcul latitude et longitude du satelite au cours de son passage    
x = np.arcsin((Re*np.cos(elevation))/(Re + H))	
ksi = np.pi/2 - elevation - x
lat = np.arcsin(np.sin(rlat)*np.cos(ksi) + np.cos(rlat)*np.sin(ksi)*np.cos(az))
lon = rlon + np.arcsin(np.sin(ksi)*np.sin(az)/np.cos(lat))

# Trouver la lon et lat du sat au moment du séisme 
indice = np.where(d == 11.2)
t_seism = d[indice]
lat_sat_seism = lat[indice]
lon_sat_seism = lon[indice]

b = data [:,4]
e = b[(int(borne_inf_time) < a) & (a < int(borne_sup_time))]
tec = e - min(e)
vtec = tec * np.cos(x)

# =============================================================================
# Animation  
# =============================================================================
N = len(lat)


graph = f1.tec_anim(d, tec, vtec)

fig = plt.figure()
ax = plt.axes(xlim =( 10, 12), ylim = (0, 12))
line, = ax.plot([], [], lw = 2)


def init():
    line.set_data([],[])
    return line, 

# animation function.  This is called sequentially
def animate(i):
    line.set_data(d[i],vtec[i])
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func = init,
                               frames = N, interval = N, blit = True)

plt.show()
    

