#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 10:57:39 2019

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
                
# =============================================================================
# Lecture des fichiers 
# =============================================================================
name_dir = os.path.join('/Users/antoineleblevec/Desktop/G20')
print (name_dir)
rep = os.path.abspath(os.path.expanduser(name_dir))
files = os.listdir(rep)

#Trie des fichiers en fonction de la latitude  
files.sort()
print(files)
#print (len(files))
#name_file_2 = []
#for j in range(0,len(files)):         
#    name_file_1 = files [j]
#    x_1 = files[j].split("_")
#    station_1 = x_1[0] 
##    print (station_1)
#    data_1 = np.loadtxt(rep + '/' + name_file_1)
#    rlon, rlat = f1.lecture_lat_lon_sat(rep, name_file_1)
##    print (rlon, rlat)
#    name_file_2.append( '{0}'.format(-rlat) + '_' + name_file_1) 
##print (name_file_2)
#name_file_2.sort()
##print (name_file_2)
#
#
#Boucle pour parcourir tous les fichiers dans le dossier 
reduce_tec_stations = []
for i in range(0,len(files)) :    
#    #stockage du nom du satellite 
    name_file = files [i]
#    print (name_file)
    x = name_file.split("_")
#    print(x)
    station = x [0]
    print (station)

#    
#    #lecture du fichier
    data = np.loadtxt(rep + '/' + name_file)
#    
#    # calcul de la longitude et lattitude de la station 
    rlon, rlat = f1.lecture_lat_lon_sat(rep, name_file)
#    print (rlon, rlat)
    
#    for z in range (1320,1393,1) :
    #    
    #    # sélection de la période d'observation du séisme
    a = data[:,1]  
    d = a [(time_inf < a) & (a < time_sup)]
    
    epoque = data[:,0]
    ep = epoque[(epoque_inf < epoque) & (epoque < epoque_sup)]
    print (ep)
    
#    for x in ep :
#        
#        print (x)
   
# elevation supérieure au sueil et pour des données dans la fenêtre du séisme  
    ele = data [:,2]  
    ele = ele [(epoque_inf < epoque) & (epoque < epoque_sup)]
    print (ele)
    #
        # élimination des satellites inintéressants
        # si ne répond pas aux critères alors itération suivante dans la boucle
#        if (len(elevation) < 30) :    
#            print("{0} bad timing or bad elevation".format(station))
#            continue
#    #    
#    #    # azimuth correspondant à la "bonne" elevation
#        az = data [:,3]
#        az = az [(time_inf < a) & (a < time_sup) & (ele > min_ele)]
#        az = np.radians(az)
#    #    
#    #    #convertir l'élévation en radians pour les calculs suivants
#        elevation = np.radians(elevation)
#    #    
#    #    # Calcul latitude et longitude du satelite au cours de son passage    
#        x = np.arcsin((Re*np.cos(elevation))/(Re + H))
#        ksi = np.pi/2 - elevation - x
#        lat = np.arcsin(np.sin(rlat)*np.cos(ksi) + np.cos(rlat)*np.sin(ksi)*np.cos(az))
#        lon = rlon + np.arcsin(np.sin(ksi)*np.sin(az)/np.cos(lat))
#    #    print (lon,lat)
#    #    
#    ##   Lon et lat du sat au moment du séisme 
#    #    indice = np.where(d == 11.03333333333)
#    #    t_seism = d[indice]
#    #    lat_sat_seism = lat[indice]
#    #    lon_sat_seism = lon[indice]
#    #    
#    ##   Calcul du tec ; enlève le tec minimal pour le "normer" ; calcul du tec vertical 
#        b = data [:,4]
#        e = b[(time_inf < a) & (a < time_sup) & (ele > min_ele)]
#        tec = e - min(e)
#        for y in range (0,len(e)-1) : 
#            if np.absolute(tec[y+1] - tec[y]) > 10 :
#    #            print ("SAUT DE TEC POUR LA STATION {0} DONNEES IGNOREES !!".format(station))
#                break 
#        else : 
#    #        print ("pas de saut de tec pour la station {0}".format(station))
#            vtec = tec * np.cos(x)  
#    #    
#    #### =============================================================================
#    #### Plot du tec 
#    #### =============================================================================
#            time_plot = a [(time_inf < a) & (a < time_sup) & (ele > min_ele)]
#    #### =============================================================================
#    #### Polynomial filter 
#    #### =============================================================================
#    ##        print (len(vtec))
#            window = 27
#            order = 4
#            vtec_polyn = savgol_filter(vtec, window, order)    
#            reduce_tec = vtec - vtec_polyn
#    #        print (reduce_tec)
#            
#            
#    
#        # remonte le prochain vtec de 0.5 à chaque fois 
#    #    reduce_tec = reduce_tec + 0.5*i
#    #    f1.polynomied_tec_1sat_stations(time_plot, reduce_tec, station)
#    #    plt.gcf()
#    #lines = plt.gca().get_lines()        
#    #labelLines(lines, xvals = [11.9]*13, align = False, fontsize = 8)
#    #plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/1sat_stations_bis_tec_polyn_filtered.png')
#    #plt.show()
