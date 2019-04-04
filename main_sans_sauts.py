#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 10:39:41 2019

@author: antoineleblevec
"""
# =============================================================================
# Modules à importer
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
# Constantes
# =============================================================================
Re = 6371032 
H = 350e3
lllat = 31.7; urlat = 43.7; lllon = 130.6; urlon = 145.6
time_inf = 10.5
time_sup = 12.5 
min_ele = 10 
                
# =============================================================================
# Lecture des fichiers 
# =============================================================================
year = str('2016')
day = str('318')
station = 'hast'
type_sat = 'GPS'
name_dir = os.path.join('/Users/antoineleblevec/Desktop/tec-suite-master_2/tec', year, day, station, type_sat)
print (name_dir)
rep = os.path.abspath(os.path.expanduser(name_dir))
files = os.listdir(rep)
print(files)

#Boucle pour parcourir tous les fichiers dans le dossier 
for i in range(0,len(files)) :
    
    #stockage du nom du satellite 
    name_file = files [i]
    x = files[i].split("_")
    sat = x [1]
    print (sat)
    
    #lecture du fichier
    data = np.loadtxt(rep + '/' + name_file)
    
    # calcul de la longitude et lattitude de la station 
    rlon, rlat = f1.lecture_lat_lon_sat(rep, name_file)
    
    # sélection de la période d'observation du séisme
    a = data[:,1]  
    d = a [(time_inf < a) & (a < time_sup)] 
    
    # elevation supérieure au sueil et pour des données dans la fenêtre du séisme  
    ele = data [:,2]    
    elevation = ele [(time_inf < a) & (a < time_sup) & (ele > min_ele)]

    # élimination des satellites inintéressants
    # si ne répond pas aux critères alors itération suivante dans la boucle
    if (len(elevation) < 30) :    
        print("%s bad timing or bad elevation" %sat)
        continue
    
    # azimuth correspondant à la "bonne" elevation
    az = data [:,3]
    az = az [(time_inf < a) & (a < time_sup) & (ele > min_ele)]
    az = np.radians(az)
    
    #convertir l'élévation en radians pour les calculs suivants
    elevation = np.radians(elevation)
    
    # Calcul latitude et longitude du satelite au cours de son passage    
    x = np.arcsin((Re*np.cos(elevation))/(Re + H))
    ksi = np.pi/2 - elevation - x
    lat = np.arcsin(np.sin(rlat)*np.cos(ksi) + np.cos(rlat)*np.sin(ksi)*np.cos(az))
    lon = rlon + np.arcsin(np.sin(ksi)*np.sin(az)/np.cos(lat))
    
#   Lon et lat du sat au moment du séisme 
    indice = np.where(d == 11.03333333333)
    t_seism = d[indice]
    lat_sat_seism = lat[indice]
    lon_sat_seism = lon[indice]
    
#   Calcul du tec ; enlève le tec minimal pour le "normer" ; calcul du tec vertical 
    b = data [:,4]
    e = b[(time_inf < a) & (a < time_sup) & (ele > min_ele)]
#    print (e)
    tec = e - min(e)
#    print (tec) 
    for y in range (0,len(e)-1) : 
        if np.absolute(tec[y+1] - tec[y]) > 10 :
            print ("SAUT DE TEC POUR LE SAT {0} DONNEES IGNOREES !!".format(sat))
            break 
    else : 
        print ("pas de saut de tec pour le sat {0}".format(sat))
        vtec = tec * np.cos(x)
#        print (tec)
#        print ("wesh")
#        print (vtec)

    
## =============================================================================
## Plot du tec 
## =============================================================================
        time_plot = a [(time_inf < a) & (a < time_sup) & (ele > min_ele)]
#    time_plot_f = time_plot[0:len(vtec_f)]
        f1.tec(time_plot, tec, vtec, sat, station, type_sat)

## =============================================================================
## Butterworth filter
## =============================================================================
#        f1.run_vtec(vtec, sat, d)

## =============================================================================
## Polynomial filter 
## =============================================================================
        print (len(vtec))
        window = 27
        order = 4
        vtec_polyn = savgol_filter(vtec, window, order)    
        reduce_tec = vtec - vtec_polyn
        f1.polynomied_tec(time_plot,vtec_polyn,reduce_tec, sat, station, type_sat, window, order)

# =============================================================================
# Plot carte traces satellites 
# =============================================================================
#    f1.plot_nz(lon, lat, sat, station, lon_sat_seism, lat_sat_seism, type_sat)
#    
#plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/traces_{1}_{0}.jpeg'.format(station,type_sat))
#plt.show()
    
## =============================================================================
## Moving average 
## =============================================================================
#    N = int( len(vtec) / 5 )
#    vtec_average = np.convolve(vtec, np.ones((N,))/N, mode='same')
##    print (vtec)
##    print(vtec_average)
#    print (len(vtec_average))
#    print (len(vtec))
#    f1.tec(d,vtec_average,vtec, sat, station)
