#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:14:48 2019

@author: antoineleblevec
"""


""" Compter le nombre de fichiers de chaque dossier.
    Les liens symboliques sont ignorés.
"""

# =============================================================================
# Description du code : Problèmes/Ajouts
# =============================================================================
## Problème :  
## dans la lecture des satellites, car les noms des satellites enregistrés par la station 
## ne s'incrémentent pas de un. 
## => Trouver un moyen de simplement lire tous les fichiers d'un dossier et de stocker leur noms. 
## créer un fichier de sortie ?  
## Ajout : 
## choix du temps d'observation 
## carte du japon 
## création d'un module fonction 


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
import fonction_1 as f1 

# =============================================================================
# Definition des constantes
# =============================================================================
Re = 6371032 
H = 350e3
lllat = 31.7; urlat = 43.7; lllon = 130.6; urlon = 145.6
                
# =============================================================================
# Interface utilisateur : choix des données à observer
# =============================================================================
year = input("enter the two last digits of the year __ chosen : ")
day = input("enter the three digits of the julian day chosen from 0 to 365 : ")
station = input("enter the four number of the station chosen : ")
name_dir = os.path.join('/Users/antoineleblevec/Desktop', year, day, station)
print (" \n Choisir les bornes sup et inf de temps d'observation du satellite :")
borne_inf = input("borne_inf :")
borne_sup = input("borne_sup :")

# =============================================================================
# Compteur du nombres de fichiers
# =============================================================================
counter = f1.Counter(name_dir)
total = 0
for cls in counter.work():
    total += cls.files

# =============================================================================
# Lecture des fichiers 
# =============================================================================
for i in range(2, total + 1) :
    sat = f1.name_sat(i)
    name_data = os.path.join('/Users/antoineleblevec/Desktop', 
                                year, day, station, station + '_' + sat + '_' + day + '_' + year + '.dat')
    data = np.loadtxt(name_data)

    (rlon, rlat) = f1.lecture_lat_lon_sat(name_data)
    
    a = data[:,1]    
    d = a [(int(borne_inf) < a) & (a < int(borne_sup))]   #temps d'observation du séisme 
    ele = data [:,2]    #sélection des stations répondant aux critères d'observation du seisme 
    elevation = ele [(int(borne_inf) < a) & (a < int(borne_sup))]
    if len(d) == 0 or np.any(elevation < 10) :    
        print("The sat %s didn't observe the seism or didn't have a good elevation" %i)
        continue
    
    az = data [:,3]
    az = az [(int(borne_inf) < a) & (a < int(borne_sup))]
    az = np.radians(az)
    
    elevation  = np.radians(elevation)
    
    x = np.arcsin((Re*np.cos(elevation))/(Re + H))	
    ksi = np.pi/2 - elevation - x
    lat = np.arcsin(np.sin(rlat)*np.cos(ksi) + np.cos(rlat)*np.sin(ksi)*np.cos(az))
    lon = rlon + np.arcsin(np.sin(ksi)*np.sin(az)/np.cos(lat))
    
    b = data [:,4]
    e = b[(int(borne_inf) < a) & (a < int(borne_sup))]
    tec = e - min(e)
    vtec = tec * np.cos(x)

    # =============================================================================
    # Plots 
    # =============================================================================

    f1.tec_sat(d, vtec, tec, i, sat)
#    f1.plot_japon(lllat, urlat, lllon, urlon, lon, lat, sat)


    
#    def f(lat,lon): 
#        return vtec
#    lon , lat = np.meshgrid(np.linspace(2.4,2.5,201),np.linspace(0.55,0.65,201))
#    vtec = f(lat,lon)
#    C = plt.contour(lon, lat, vtec, 20)
#    plt.colorbar()

#plt.savefig('/Users/antoineleblevec/Desktop/2011_seism_data/fig/pos_sat_lon.png')
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
#    name_file = "/Users/antoineleblevec/Desktop/2011_seism_data/data/{0}_{1}_{2}_{3}_new.dat".format(station, sat, day, year)
#    
#    destination = open(name_file, "w")
#    np.savetxt(name_file, np.c_[tsn,d,elevation, az, tec, vtec], \
#               header = "Columns : tsn, d, elevation, az, slant_tec, vtec", \
#               fmt =  '%1.5f')
#    destination.close()

