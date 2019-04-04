#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 17:02:55 2019

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
year = str('2016')
day = str('318')
station = 'mqzg'
name_dir = os.path.join('/Users/antoineleblevec/Desktop/tec-suite-master_2/tec', year, day, station)
#print (" \n Choisir les bornes sup et inf de temps d'observation du satellite :")
#borne_inf_time = input("borne_inf_time :")
#borne_sup_time = input("borne_sup_time :")
borne_inf_time = 10.30
borne_sup_time = 12.30 
sueil_elevation = 10 

# =============================================================================
# Lecture des fichiers 
# =============================================================================


# =============================================================================
# pour glonass
# =============================================================================
rep = os.path.abspath(os.path.expanduser(name_dir))
files = os.listdir(rep)
for i in range(0,len(files)) :
    name_data = files [i]
    x = files[i].split("_")
    sat = x [1]
    
    if sat[0] == 'R' : 
        data_r = np.loadtxt(rep + '/' + name_data)
        
        # calcul de la longitude et lattitude de la station 
        rlon_r, rlat_r = f1.lecture_lat_lon_sat(rep, name_data)
        
        a_r = data_r[:,1] 
        # temps d'observation du séisme 
        d_r = a_r [(int(borne_inf_time) < a_r) & (a_r < int(borne_sup_time))] 
        # elevation  
        ele_r = data_r [:,2]
        # elevation au moment du séisme    
        elevation_r = ele_r [(int(borne_inf_time) < a_r) & (a_r < int(borne_sup_time))]
    
        if len(d_r) == 0 or np.any(elevation_r < 10) :    
            print("%s wasn't above the seism or didn't have the requested elevation" %sat)
            continue
        
        az_r = data_r [:,3]
        az_r = az_r [(int(borne_inf_time) < a_r) & (a_r < int(borne_sup_time))]
        az_r = np.radians(az_r)
        
        elevation_r  = np.radians(elevation_r)
        
        # Calcul latitude et longitude du satelite au cours de son passage    
        x_r = np.arcsin((Re*np.cos(elevation_r))/(Re + H))	
        ksi_r = np.pi/2 - elevation_r - x_r
        lat_r = np.arcsin(np.sin(rlat_r)*np.cos(ksi_r) + np.cos(rlat_r)*np.sin(ksi_r)*np.cos(az_r))
        lon_r = rlon_r + np.arcsin(np.sin(ksi_r)*np.sin(az_r)/np.cos(lat_r))
        
        # Trouver la lon et lat du sat au moment du séisme 
        indice_r = np.where(d_r == 11.2)
        lat_sat_seism_r = lat_r[indice_r]
        lon_sat_seism_r = lon_r[indice_r]
        
        # Calcul du tec ; enlève le tec minimal pour le "normer" ; calcul du tec vertical 
        b_r = data_r [:,4]
        e_r = b_r[(int(borne_inf_time) < a_r) & (a_r < int(borne_sup_time))]
        tec_r = e_r - min(e_r)
        vtec_r = tec_r * np.cos(x_r)
        
    ## =============================================================================
    ## Plot basique 
    ## =============================================================================
    #    f1.tec(d, tec, vtec, sat, station)
    
    # =============================================================================
    # Plot traces satellites 
    # =============================================================================
        f1.plot_nz_r(lon_r, lat_r, sat, station, lon_sat_seism_r, lat_sat_seism_r)

    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/traces_sat_{0}.png'.format(station))
    plt.show()
    
    
    
## =============================================================================
## Butterworth filter
## =============================================================================
##    f1.run_vtec(vtec, sat, d)
#    
#       
## =============================================================================
## Moving average 
## =============================================================================
##    N = int( len(vtec) / 5 )
##    vtec_average = np.convolve(vtec, np.ones((N,))/N, mode='same')
###    print (vtec)
###    print(vtec_average)
##    print (len(vtec_average))
##    print (len(vtec))
##    f1.tec(d,vtec_average,vtec, sat)
#
## =============================================================================
## Polynomial filter 
## =============================================================================
##    vtec_polyn = savgol_filter(vtec, 23, 10)
##    print (vtec)
##    print (vtec_polyn)
##    f1.tec(d,vtec, vtec_polyn, sat)
#    
##    reduce_tec = vtec - vtec_polyn
##    f1.tec(d,vtec_polyn,reduce_tec, sat)
#
## =============================================================================
## Bonuses    
## =============================================================================
##    def f(lat,lon): 
##        return vtec
##    lon , lat = np.meshgrid(np.linspace(2.4,2.5,201),np.linspace(0.55,0.65,201))
##    vtec = f(lat,lon)
##    C = plt.contour(lon, lat, vtec, 20)
##    plt.colorbar()
#
##plt.savefig('/Users/antoineleblevec/Desktop/11_results/fig/pos_sat_lon.png')
##plt.show()
##plt.close()
#        
## =============================================================================
##     Création des fichiers de sortie 
## =============================================================================
##    
##    tsn = data[:,0]
##    tsn = tsn [(int(borne_inf_time) < a) & (a < int(borne_sup_time))] #époque d'observation du séisme 
##
##    
##    name_file = "/Users/antoineleblevec/Desktop/11_results/data/{0}_{1}_{2}_{3}_new.dat".format(station, sat, day, year)
##    
##    destination = open(name_file, "w")
##    np.savetxt(name_file, np.c_[tsn,d,elevation, az, tec, vtec], \
##               header = "Columns : tsn, d, elevation, az, slant_tec, vtec", \
##               fmt =  '%1.5f')
##    destination.close()
#
