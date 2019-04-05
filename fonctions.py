#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 11:34:21 2019

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
import fonction_1 as f1
from scipy.signal import butter, lfilter
import matplotlib.animation as animation
from labellines import labelLine, labelLines


# =============================================================================
# Varion TEC and VTEC from RINEX fonction
# =============================================================================
def tec(d, courbe_1, courbe_2, sat, station, type_sat) : 
    plt.xlabel('Temps IUT (h)')
    plt.ylabel('TEC (TECU)')
    plt.title('TEC and VTEC measured by {2}: {0} and recorded by : {1}'.format(sat,station,type_sat))
    plt.plot(d, courbe_1, label = "{0}".format('tec'), color = 'green')
    plt.plot(d, courbe_2, label = "{0}".format('vtec'), color = 'blue')
    plt.axvline(x = 11.033, color = 'red', label = "seism", linestyle = ':')
    delta_tec = plt.legend()
    plt.gcf()
    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/{0}_{1}_tec.png'.format(sat,station))
    plt.show() 
    return delta_tec

# =============================================================================
# Varion TEC and VTEC from RINEX fonction for uniq_data
# =============================================================================
def tec_uniq(d, courbe_1, courbe_2) : 
    plt.xlabel('Temps IUT (h)')
    plt.ylabel('TEC (TECU)')
    plt.title('TEC and VTEC')
    plt.plot(d, courbe_1, label = "tec", color = 'green')
    plt.plot(d, courbe_2, label = "vtec", color = 'blue')
    plt.axvline(x = 11.033, color = 'red', label = "seism", linestyle = ':')
    delta_tec = plt.legend()
    plt.gcf()
    plt.show() 
    return delta_tec

# =============================================================================
# Plot TEC d'un satellite pour plusieurs stations 
# =============================================================================
def tec_1sat_stations(d, courbe_1, courbe_2,station) : 
    plt.xlabel('Temps IUT (h)')
    plt.ylabel('TEC (TECU)')
    plt.title('TEC and VTEC measured by {0} and recorded by G20'.format(station))
    plt.plot(d, courbe_1, label = "{0}".format('tec'), color = 'green')
    plt.plot(d, courbe_2, label = "{0}".format('vtec'), color = 'blue')
    plt.axvline(x = 11.033, color = 'red', label = "seism", linestyle = ':')
    delta_tec = plt.legend()
    plt.gcf()
    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/1sat_stations_tec.png')
    plt.show() 
    return delta_tec


# =============================================================================
# Fonction plot filtered TEC 
# =============================================================================
def polynomied_tec(time_plot, courbe_1, courbe_2, sat, station, type_sat, window, order) : 
    # plot du tec polynomial et du vtec - vtecpolynomial 
    plt.xlabel('Temps IUT (h)')
    plt.ylabel('TEC (TECU)')
    plt.title('TEC and VTEC measured by {2}: {0} and recorded by : {1}'.format(sat,station, type_sat))
    plt.plot(time_plot, courbe_1, label = "vtec_polyn", color = 'red')
    plt.plot(time_plot, courbe_2, label = "vtec_reduced", color = 'blue')
    # ligne verticale représentant le séisme 
    plt.axvline(x = 11.033, color = 'k', label = "seism", linestyle = ':')
    delta_tec = plt.legend()
    plt.gcf()
    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/{0}_{1}_tec_polyn_filtered_compar.png'.format(sat,station))
    plt.show()
    plt.plot(time_plot, courbe_2, label = "vtec_filtered", color = 'b')
    plt.axis([10.5, 12, -0.25, 0.25])
    plt.xlabel('Temps IUT (h)')
    plt.ylabel('TEC (TECU)')
    plt.title('Polynomied VTEC: (window {3}, {4}th order) measured by {2}: {0} recorded by : {1}'.format(sat,station, type_sat,window,order))
    plt.axvline(x = 11.033, color = 'k', label = "seism", linestyle = ':')
    delta_tec = plt.legend()
    plt.gcf()
    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/{0}_{1}_tec_polyn_filtered.png'.format(sat,station))
    plt.show()
    return delta_tec

# =============================================================================
# Fonction plot filtered TEC 
# =============================================================================
def polynomied_tec_1sat_stations(time_plot, courbe_2, station) : 

    plt.plot(time_plot, courbe_2, label = "{0}".format(station))
#    labelLines(plt.gca().get_lines(),align = True)
    plt.axis([10.5, 12, -0.5, 20])
    plt.xlabel('Temps IUT (h)')
    plt.ylabel('TEC (TECU)')
    plt.title('Polynomied VTEC of R07 for stations with latitudes smaller than the epicentre')
    delta_tec = plt.axvline(x = 11.033, color = 'k', linestyle = ':')
#    delta_tec = plt.legend(bbox_to_anchor=(1, 1), loc='upper left', 
#                   borderaxespad=0., fancybox = True, shadow = True)
#    plt.gcf()
#    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/1sat_stations_bis_tec_polyn_filtered.png')
#    plt.show()
    return delta_tec

# =============================================================================
# Calcul longitude et lattitude de la station
# =============================================================================
def lecture_lat_lon_sat(rep, name_data): 
        with open(rep + '/' + name_data,'r') as fich:
            for j in range(7):
                ligne = fich.readline()
        liste_mots = ligne.split(":")
        mot=liste_mots[1]
        liste_bis_mots = mot.split(",")
        rlon = np.radians(float(liste_bis_mots[0])) #longitude de la station
        rlat = np.radians(float(liste_bis_mots[1])) #lattitude de la station
        return(rlon, rlat)
 
# =============================================================================
# Calcul longitude et lattitude de la station test un fichier
# =============================================================================
def lecture_lat_lon_sat_uniq(files): 
        with open(files,'r') as fich:
            for j in range(7):
                ligne = fich.readline()
        liste_mots = ligne.split(":")
        mot=liste_mots[1]
        liste_bis_mots = mot.split(",")
        rlon = np.radians(float(liste_bis_mots[0])) #longitude de la station
        rlat = np.radians(float(liste_bis_mots[1])) #lattitude de la station
        return(rlon, rlat)    

# =============================================================================
# lecture longitude
# =============================================================================
def lecture_lon(files): 
    with open(files,'r') as fich:
        for j in range(7):
            ligne = fich.readline()
        liste_mots = ligne.split(":")
        mot=liste_mots[1]
        liste_bis_mots = mot.split(",")
        rlon = np.radians(float(liste_bis_mots[0]))
        return rlon 

# =============================================================================
# lecture latitude
# =============================================================================
def lecture_lat(files): 
    with open(files,'r') as fich:
        for j in range(7):
            ligne = fich.readline()
        liste_mots = ligne.split(":")
        mot=liste_mots[1]
        liste_bis_mots = mot.split(",")
        rlat = np.radians(float(liste_bis_mots[1]))
        return rlat 

# =============================================================================
# Filtre passe bande butterworth
# =============================================================================
def butter_bandpass(lowcut, highcut, fs, order = 6):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order = 6):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


def run_vtec(vtec, sat, d):
    # Sample rate and desired cutoff frequencies (in Hz).
    fs = 0.008
    lowcut = 0.0007
    highcut = 0.003
    order = 6
    y = butter_bandpass_filter(vtec, lowcut, highcut, fs, order)
    #plot
    plt.plot(d, y, label = 'Filtered vtec')
    plt.plot(d, vtec, label = 'Vtec', color = 'g')
    plt.xlabel('Temps UT (h)')
    plt.grid(True)
    plt.axis('tight')
    plt.title('Variation of TEC and filtered TEC during seism observed by satelite : %s ' %sat)
    plt.axvline(x = 11.033, color = 'red', label = 'seism', linestyle = ':')
    plt.legend(loc='upper left')
    plt.gcf()
    plt.show()
    plt.plot(d, y, label = 'Filtered vtec')
    plt.title('Buttered {0}th order TEC observed by : {1} '.format(order,sat))
    plt.grid(True)
    plt.axvline(x = 11.033, color = 'red', label = "seism", linestyle = ':')
    plt.xlabel('Temps UT (h)')
    plt.legend(loc='upper left')
    plt.gcf()
    plt.show()


# =============================================================================
# Création du fond de carte de NZ
# =============================================================================
def basic_nz_map() : 
    
    # =============================================================================
    # Lecture du fichier en deux fois pour problème de type
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
    # Plot du fond de carte de nouvelle-zélande
    # =============================================================================
    m = Basemap(width=1500000,height=1500000,projection='lcc',
                    resolution='c',lat_1=-80.,lat_2=-20,lat_0=-40,lon_0=176.)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    parallels = np.arange(-60, -30, 5.)
    m.drawparallels(parallels, labels = [False,True,True,False], dashes = [1,5], linewidth = 0.5)
    meridians = np.arange(160., 190., 5.)
    m.drawmeridians(meridians, labels = [True, False, False, True], dashes = [1,5], linewidth = 0.5 )
    m.bluemarble()
    # =============================================================================
    # Plot des stations sur le fond de carte  
    # =============================================================================
    for j in range (len(lon)) : 
        xpt,ypt = m(lon[j],lat[j])
        m.plot(xpt,ypt,'wo', ms = 1) 
        plt.text(xpt+10000,ypt+10000,'%s' %station[j], color = 'white', fontsize = 5)
        expt, eypt = m(elon,elat)
        m.plot(expt,eypt,'r*', ms = 7)

    return m 

# =============================================================================
#Fonction plot Nouvelle-Zélande 
# =============================================================================
def plot_nz(lon, lat, sat, station, lon_sat_seism, lat_sat_seism, type_sat) : 
 
        plt.title('Trace des satellites {1} au cours du séisme vu par {0}'.format(station, type_sat))
        m = basic_nz_map()
        x, y = m(np.degrees(lon), np.degrees(lat))
        m.plot(x, y, label = "{0}".format(sat))
        u, v = m(np.degrees(lon_sat_seism), np.degrees(lat_sat_seism))
        m.plot(u, v, 'w*', ms = 5)
        plt.legend(bbox_to_anchor=(1.15, 1), loc='upper left', 
                   borderaxespad=0., fancybox = True, shadow = True)
        return m
# =============================================================================
# Fonction minimum
# =============================================================================
def minimum(x) : 
    return (x.min())

def maximum(x) : 
    return(x.max())



# =============================================================================
#  Animation 
# =============================================================================
# =============================================================================
# Calcul longitude et lattitude de la station pour module animation
# =============================================================================
def lecture_lat_lon_sat_anim(name_dir, files): 
        with open(name_dir + '/' + files[14],'r') as fich:
            for j in range(7):
                ligne = fich.readline()
        liste_mots = ligne.split(":")
        mot=liste_mots[1]
        liste_bis_mots = mot.split(",")
        rlon = np.radians(float(liste_bis_mots[0])) #longitude de la station
        rlat = np.radians(float(liste_bis_mots[1])) #lattitude de la station
        return(rlon, rlat)
# =============================================================================
#Fonction plot Nouvelle-Zélande pour l'animation
# =============================================================================
def plot_nz_anim(lon, lat, station, lon_sat_seism, lat_sat_seism, type_sat) : 
 
        
        bmap = basic_nz_map()
        plt.title('Trace du R07')
        x, y = bmap(np.degrees(lon), np.degrees(lat))
        
#        m.plot(x, y)
#        u, v = m(np.degrees(lon_sat_seism), np.degrees(lat_sat_seism))
#        m.plot(u, v, 'w*', ms = 5)
#        plt.legend(bbox_to_anchor=(1.3, 1), loc='upper left', 
#                   borderaxespad=0., fancybox = True, shadow = True)
        
        return bmap
# =============================================================================
# Création du fond de carte de NZ pour l'animation
# =============================================================================
def basic_nz_map_anim() : 
    
    # =============================================================================
    # Lecture du fichier en deux fois pour problème de type
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
    # Plot du fond de carte de nouvelle-zélande
    # =============================================================================
    m = Basemap(width=1500000,height=1500000,projection='lcc',
                    resolution='c',lat_1=-80.,lat_2=-20,lat_0=-40,lon_0=176.)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    parallels = np.arange(-60, -30, 5.)
    m.drawparallels(parallels, labels = [False,True,True,False], dashes = [1,5], linewidth = 0.5)
    meridians = np.arange(160., 190., 5.)
    m.drawmeridians(meridians, labels = [True, False, False, True], dashes = [1,5], linewidth = 0.5 )
    m.bluemarble()
    # =============================================================================
    # Plot des stations sur le fond de carte  
    # =============================================================================
    for j in range (len(lon)) : 
        xpt,ypt = m(lon[j],lat[j])
        m.plot(xpt,ypt,'wo', ms = 1) 
        plt.text(xpt+10000,ypt+10000,'%s' %station[j], color = 'white', fontsize = 5)
        expt, eypt = m(elon,elat)
        m.plot(expt,eypt,'r*', ms = 7)

    return m 
# =============================================================================
# Varion TEC and VTEC from RINEX fonction for animation
# =============================================================================
def tec_anim(d, courbe_1, courbe_2) : 
    plt.xlabel('Temps IUT (h)')
    plt.ylabel('TEC (TECU)')
    plt.title('TEC and VTEC measured by Glonass R21')
    plt.plot(d, courbe_1, label = "tec", color = 'green')
    plt.plot(d, courbe_2, label = "vtec", color = 'blue')
    plt.axvline(x = 11.033, color = 'red', label = "seism", linestyle = ':')
    delta_tec = plt.legend()
#    plt.gcf()
#    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/anim_tec.png')
#    plt.show() 
    return delta_tec








#
## =============================================================================
## Fonctions pour étudier le Japon 
## =============================================================================
#
#
## =============================================================================
## Création de la carte du  Japon 
## =============================================================================
#def basic_japon_map(lllat = 31.7, urlat = 43.7, lllon = 130.6, urlon = 145.6) : 
#    m = Basemap(projection='stere', 
#                lon_0 = (urlon + lllon) / 2, 
#                lat_0 = (urlat + lllat) / 2, 
#                llcrnrlat = lllat, urcrnrlat = urlat, 
#                llcrnrlon = lllon, urcrnrlon = urlon, 
#                resolution = 'c')
#    m.drawcoastlines()
#    m.drawstates()
#    m.drawcountries()
#    m.bluemarble()
#    return m 
#
## =============================================================================
##    Def fonction plot Japon  
## =============================================================================
#def plot_japon(lllat, urlat, lllon, urlon, lon, lat, sat) : 
#
#    axes = plt.gca()
#    plt.xlabel('Longitude')
#    axes.set_xlim(30, 45)
#    plt.ylabel('Latitude')
#    plt.title('Mouvement des satellites au cours du séisme')
#    m = basic_japon_map(lllat = lllat, urlat = urlat, 
#                        lllon = lllon, urlon = urlon)
#    x, y = m(np.degrees(lon),np.degrees(lat))
#    jap = m.plot(x, y, label = "%s" %sat)
#    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', 
#               borderaxespad=0., fancybox = True, shadow = True)
#    return jap 