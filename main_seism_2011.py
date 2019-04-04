#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 16:35:02 2019

@author: antoineleblevec
"""
# =============================================================================
# Description du code 
# =============================================================================
## Problème :   
## créer un fichier de sortie ?  
## Ajout : 
## choix du temps d'observation 
## carte du japon 
## création d'un module fonction 
## chaque dossier peut être parcouru et les fichiers lu un par un 
## plusieurs fonctions de filtres du signal ont été ajouté 


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


# =============================================================================
# Definition des constantes
# =============================================================================
Re = 6371032 
H = 350e3
lllat = 31.7; urlat = 43.7; lllon = 130.6; urlon = 145.6
                
# =============================================================================
# Interface utilisateur : choix des données à observer
# =============================================================================
#year = input("enter the two last digits of the year __ chosen : ")
#day = input("enter the three digits of the julian day chosen from 0 to 365 : ")
#station = input("enter the four number of the station chosen : ")
year = str('2011')
day = str('070')
station = str ('3009')
name_dir = os.path.join('/Users/antoineleblevec/Desktop/tec-suite-master_2/tec', year, day, station)
#print (" \n Choisir les bornes sup et inf de temps d'observation du satellite :")
#borne_inf = input("borne_inf :")
#borne_sup = input("borne_sup :")
borne_inf = 5
borne_sup = 8 

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
    
    rlon, rlat = f1.lecture_lat_lon_sat(rep, name_data)

    
    a = data[:,1]    
    d = a [(int(borne_inf) < a) & (a < int(borne_sup))]   #temps d'observation du séisme 
    ele = data [:,2]    #sélection des stations répondant aux critères d'observation du seisme 
    elevation = ele [(int(borne_inf) < a) & (a < int(borne_sup))]
    if len(d) == 0 or np.any(elevation < 10) :    
        print("%s didn't observe the seism or didn't have a good elevation" %sat)
        continue
    
    az = data [:,3]
    az = az [(int(borne_inf) < a) & (a < int(borne_sup))]
    az = np.radians(az)
    
    elevation  = np.radians(elevation)
    
    
    # Calcul latitude et longitude du satellite au cours de son passage 
    x = np.arcsin((Re*np.cos(elevation))/(Re + H))	
    ksi = np.pi/2 - elevation - x
    lat = np.arcsin(np.sin(rlat)*np.cos(ksi) + np.cos(rlat)*np.sin(ksi)*np.cos(az))
    lon = rlon + np.arcsin(np.sin(ksi)*np.sin(az)/np.cos(lat))
    
    b = data [:,4]
    e = b[(int(borne_inf) < a) & (a < int(borne_sup))]
    tec = e - min(e)
    vtec = tec * np.cos(x)
    
# =============================================================================
# Plot basique 
# =============================================================================
#    f1.tec(d, vtec, tec, sat)
    
    
# =============================================================================
#   Plot sur carte Japon
# =============================================================================
    f1.plot_japon(lllat, urlat, lllon, urlon, lon, lat, sat)
    
plt.savefig('/Users/antoineleblevec/Desktop/11_results/fig/70/0205/pos_sat/pos.png')

    
# =============================================================================
# Butterworth filter
# =============================================================================
#    f1.run_vtec(vtec, sat, d)
    
       
# =============================================================================
# Moving average 
# =============================================================================
#    N = int( len(vtec) / 5 )
#    vtec_average = np.convolve(vtec, np.ones((N,))/N, mode='same')
##    print (vtec)
##    print(vtec_average)
#    print (len(vtec_average))
#    print (len(vtec))
#    f1.tec(d,vtec_average,vtec, sat)

# =============================================================================
# Polynomial filter 
# =============================================================================
#    vtec_polyn = savgol_filter(vtec, 23, 10)
#    print (vtec)
#    print (vtec_polyn)
#    f1.tec(d,vtec, vtec_polyn, sat)
    
#    reduce_tec = vtec - vtec_polyn
#    f1.tec(d,vtec_polyn,reduce_tec, sat)

# =============================================================================
# Bonuses    
# =============================================================================
#    def f(lat,lon): 
#        return vtec
#    lon , lat = np.meshgrid(np.linspace(2.4,2.5,201),np.linspace(0.55,0.65,201))
#    vtec = f(lat,lon)
#    C = plt.contour(lon, lat, vtec, 20)
#    plt.colorbar()

#plt.savefig('/Users/antoineleblevec/Desktop/11_results/fig/pos_sat_lon.png')
#plt.show()
#plt.close()
        
# =============================================================================
#     Création des fichiers de sortie 
# =============================================================================
#    
#    tsn = data[:,0]
#    tsn = tsn [(int(borne_inf) < a) & (a < int(borne_sup))] #époque d'observation du séisme 
#
#    
#    name_file = "/Users/antoineleblevec/Desktop/11_results/data/{0}_{1}_{2}_{3}_new.dat".format(station, sat, day, year)
#    
#    destination = open(name_file, "w")
#    np.savetxt(name_file, np.c_[tsn,d,elevation, az, tec, vtec], \
#               header = "Columns : tsn, d, elevation, az, slant_tec, vtec", \
#               fmt =  '%1.5f')
#    destination.close()

