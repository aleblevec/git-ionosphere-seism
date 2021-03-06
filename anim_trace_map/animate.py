#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 12:00:51 2019

@author: antoineleblevec
"""

# =============================================================================
# Description du code 
# =============================================================================

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
# Choix des données à observer 
# =============================================================================
#year = input("Enter date of the year chosen : ")
#day = input("Enter the julian day chosen from 0 to 365 : ")
#station = input("Enter the name of the station chosen : ")
#type_sat = input("Choisir GPS ou Glonass :")
year = str('2016')
day = str('318')
station = 'lkta'
type_sat = 'Glonass'
name_dir = os.path.join('/Users/antoineleblevec/Desktop/tec-suite-master_2/tec', year, day, station, type_sat)
print (name_dir)
#print (" \n Choisir les bornes sup et inf de temps d'observation du satellite :")
#borne_inf_time = input("borne_inf_time :")
#borne_sup_time = input("borne_sup_time :")
borne_inf_time = 10.30
borne_sup_time = 12.30 
sueil_elevation = 10 

# =============================================================================
# Lecture des fichiers 
# =============================================================================

rep = os.path.abspath(os.path.expanduser(name_dir))
files = os.listdir(rep)
for i in range(0,len(files)) :
    name_data = files [i]
    x = files[i].split("_")
    sat = x [1]
    data = np.loadtxt(rep + '/' + name_data)
    
    # calcul de la longitude et lattitude de la station 
    rlon, rlat = f1.lecture_lat_lon_sat(rep, name_data)
    
    a = data[:,1] 
    # temps d'observation du séisme 
    d = a [(int(borne_inf_time) < a) & (a < int(borne_sup_time))] 
    # elevation  
    ele = data [:,2]
    # elevation au moment du séisme    
    elevation = ele [(int(borne_inf_time) < a) & (a < int(borne_sup_time))]

    if len(d) == 0 or np.any(elevation < 10) :    
        print("%s wasn't above the seism or didn't have the requested elevation" %sat)
        continue
    
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
    
    # Calcul du tec ; enlève le tec minimal pour le "normer" ; calcul du tec vertical 
    b = data [:,4]
    e = b[(int(borne_inf_time) < a) & (a < int(borne_sup_time))]
    tec = e - min(e)
    vtec = tec * np.cos(x)
    
## =============================================================================
## Plot basique 
## =============================================================================
#    f1.tec(d, tec, vtec, sat, station)

# =============================================================================
# Essai d'animation 
# =============================================================================



# =============================================================================
# Plot traces satellites 
# =============================================================================
    m = f1.plot_nz_anim(lon, lat, sat, station, lon_sat_seism, lat_sat_seism, type_sat)

    
   
        
plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/anim_traces_{1}_{0}.jpeg'.format(station,type_sat))
plt.show()
