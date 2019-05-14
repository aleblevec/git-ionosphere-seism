#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 11:34:21 2019

@author: antoineleblevec

ajouter :     for j in range (len(col_names_station)):
        xpt,ypt = m(np.degrees(lo_station[j]),np.degrees(la_station[j]))
        m.plot(xpt,ypt,'wo', ms = 1) 
        plt.text(xpt+10000,ypt+10000,'%s' %col_names_station[j], color = 'white', fontsize = 5)

    => pour plotter les stations sur la carte 

"""

# =============================================================================
# Modules 
# =============================================================================
import numpy as np
#import os.path
#from os import scandir
import matplotlib.pyplot as plt
#import csv
#from os import chdir 
#from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap
#from scipy.signal import butter, lfilter
#import matplotlib.animation as animation
#from labellines import labelLine, labelLines
import pandas as pd

# =============================================================================
# DataFrames from RINEX
# =============================================================================
def read(files, epoque1, epoque2):
    df = pd.read_csv(
            files,
            skiprows=11,
            delim_whitespace=True,
            header=None,
            names=["tsn", "hour", "el", "az", "tec", "tec.p1p2", "validity"]
            )
    
    # column of interest
    df = df[["tsn", "el", "az", "tec"]] 
    
    df = df[epoque1 < df["tsn"]]
    return df[df["tsn"] < epoque2]

# =============================================================================
# DataFrames from Sound_speed
# =============================================================================
def read_sound(files):
    df = pd.read_csv(
            files,
            delim_whitespace=True,
            header=None,
            names=['a','b','H','c','Cs','d','e']
            )
    
    df = df[['H','Cs']]
    
    return df
# =============================================================================
# Latitude et longitude station from RINEX
# =============================================================================
def lecture_lat_lon(files): 
    rlat = [] ; rlon = []
    with open(files,'r') as fich:
        for j in range(7):
            ligne = fich.readline()
        liste_mots = ligne.split(":")
        mot=liste_mots[1]
        liste_bis_mots = mot.split(",")
        rlat = np.radians(float(liste_bis_mots[1]))
        rlon = np.radians(float(liste_bis_mots[0]))
        return rlat,rlon 

## =============================================================================
## Tohoku Functions
## =============================================================================
## Basemap Japan ; rajouter lon_station, lat_station et station dans les paramètres 
## pour voir les stations sur la carte et décommenter l'itération
def japan_map(lllat,urlat,lllon,urlon,elon,elat,lon_sip_max,lat_sip_max,tda,
              H,v_i_s,v_s_s,sat) : 
    fig = plt.figure()
    m = Basemap(projection='stere', 
                lon_0 = (urlon + lllon) / 2, 
                lat_0 = (urlat + lllat) / 2, 
                llcrnrlat = lllat, urcrnrlat = urlat, 
                llcrnrlon = lllon, urcrnrlon = urlon, 
                resolution = 'c')
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.bluemarble()
    expt, eypt = m(elon,elat)
    m.plot(expt,eypt,'r*', ms = 7)
    x, y = m(np.degrees(lon_sip_max), np.degrees(lat_sip_max))
    m.hexbin(x,
             y,
             C=tda,
             reduce_C_function=np.mean,
             gridsize=30, 
             cmap="plasma")
    plt.title('Max VTEC Hion: {0}m ; {1}<seuil<{2} ; sat: {3}'.format(H,v_i_s,v_s_s,sat))
    cbaxes = fig.add_axes([0.90,0.15,0.01,0.7]) 
    plt.colorbar(cax=cbaxes)
    m.colorbar()
    plt.gcf()
    plt.show()
    
    return m 

def tda(station,df_vtec,v_i_s,v_s_s,df_lon,df_lat): 
    a = df_vtec['tec_{0}'.format(station)]
    for i in range(len(a)-1): 
        if abs(a.iloc[i+1]-a.iloc[i]) > v_i_s and abs(a.iloc[i+1]-a.iloc[i]) < v_s_s :
            break
    # lon_tda = df_lon['el_{0}'.format(station)].iloc[i]
    # lat_tda = df_lat['el_{0}'.format(station)].iloc[i]
    # vtec_tda = a.iloc[i]
    return i
    # , lon_tda, lat_tda, vtec_tda



# détection de la première onde, en prenant les indices à partir de 0. 
def detec_first_onde(station,df_vtec,v_i_s,v_s_s,df_lon,df_lat): 
    a = df_vtec['tec_{0}'.format(station)]
    for i in range(len(a)-1): 
        if abs(a.iloc[i+1]-a.iloc[i]) > v_i_s and abs(a.iloc[i+1]-a.iloc[i]) < v_s_s :
            break
    lon_tda = df_lon['el_{0}'.format(station)].iloc[i]
    lat_tda = df_lat['el_{0}'.format(station)].iloc[i]
    vtec_tda = a.iloc[i]
    return i, lon_tda, lat_tda, vtec_tda

# max vtec avec longitude et latitude correspondante, 100s après la détection de la première onde
def lon_lat_tecmax(station,df_vtec,p_p,e1,df_lon,df_lat,v_i_s,v_s_s):
    a = df_vtec['tec_{0}'.format(station)]
    a = a[detec_first_onde(station,df_vtec,v_i_s,v_s_s) : detec_first_onde(station,df_vtec,v_i_s,v_s_s)+p_p]
    indice_tec = a.idxmax() - (e1+1)
    lon_max = df_lon['el_{0}'.format(station)].iloc[indice_tec]
    lat_max = df_lat['el_{0}'.format(station)].iloc[indice_tec]
    return lon_max, lat_max, a.max()

## =============================================================================
## Butterworth
## =============================================================================
#def butter_bandpass(lowcut, highcut, fs, order = 6):
#    nyq = 0.5 * fs
#    low = lowcut / nyq
#    high = highcut / nyq
#    b, a = butter(order, [low, high], btype='band')
#    return b, a
#
#def butter_bandpass_filter(data, lowcut, highcut, fs, order = 6):
#    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
#    y = lfilter(b, a, data)
#    return y
#
#
#def run_vtec(vtec, sat, d):
#    # Sample rate and desired cutoff frequencies (in Hz).
#    fs = 0.008
#    lowcut = 0.0007
#    highcut = 0.003
#    order = 6
#    y = butter_bandpass_filter(vtec, lowcut, highcut, fs, order)
#    #plot
#    plt.plot(d, y, label = 'Filtered vtec')
#    plt.plot(d, vtec, label = 'Vtec', color = 'g')
#    plt.xlabel('Temps UT (h)')
#    plt.grid(True)
#    plt.axis('tight')
#    plt.title('Variation of TEC and filtered TEC during seism observed by satelite : %s ' %sat)
#    plt.axvline(x = 11.033, color = 'red', label = 'seism', linestyle = ':')
#    plt.legend(loc='upper left')
#    plt.gcf()
#    plt.show()
#    plt.plot(d, y, label = 'Filtered vtec')
#    plt.title('Buttered {0}th order TEC observed by : {1} '.format(order,sat))
#    plt.grid(True)
#    plt.axvline(x = 11.033, color = 'red', label = "seism", linestyle = ':')
#    plt.xlabel('Temps UT (h)')
#    plt.legend(loc='upper left')
#    plt.gcf()
#    plt.show()

## =============================================================================
## Hodochrones 
## =============================================================================
#def hodor(lon, lat, sat): 
# 
#        plt.title('Trace et mesure du TEC par les satellites {0} au cours du séisme '.format(sat))
#        m = basic_nz_map()
#        x, y = m(lon, lat)
#        m.plot(x, y, label = "{0}".format(sat))
##        u, v = m(lon_seism, lat_seism)
##        m.plot(u, v, 'w*', ms = 5)
#        plt.legend(bbox_to_anchor=(1.15, 1), loc='upper left', 
#                   borderaxespad=0., fancybox = True, shadow = True)
#        return m

## =============================================================================
## Basemap NZ
## =============================================================================
#def basic_nz_map() : 
#    
#    # =============================================================================
#    # Lecture du fichier en deux fois pour problème de type
#    # =============================================================================
#    datat = np.loadtxt('/Users/antoineleblevec/Desktop/2016_seism/coord_stations.txt', usecols = (1,2,3,4,5))
#    data = open('/Users/antoineleblevec/Desktop/2016_seism/coord_stations.txt', "r")
#    # =============================================================================
#    # Stockage des noms des stations dans : station 
#    # =============================================================================
#    contenu = data.readlines()
#    station = []
#    for j in range (36) : 
#        list = contenu[j].split()
#        station.append(list [0])
#    # =============================================================================
#    # Stockage des latitudes et longitudes des stations dans lon et lat
#    # =============================================================================
#    lon = []
#    lat = []
#    for i in datat : 
#        lon.append(i[4])
#        lat.append(i[3])
#    # =============================================================================
#    # Epicentre longitude latitude 
#    # =============================================================================
##    elon = 173.054
##    elat = -42.737
#    # =============================================================================
#    # Plot du fond de carte de nouvelle-zélande
#    # =============================================================================
#    m = Basemap(width=1500000,height=1500000,projection='lcc',
#                    resolution='c',lat_1=-80.,lat_2=-20,lat_0=-40,lon_0=176.)
#    m.drawcoastlines()
#    m.drawstates()
#    m.drawcountries()
#    parallels = np.arange(-60, -30, 5.)
#    m.drawparallels(parallels, labels = [False,True,True,False], dashes = [1,5], linewidth = 0.5)
#    meridians = np.arange(160., 190., 5.)
#    m.drawmeridians(meridians, labels = [True, False, False, True], dashes = [1,5], linewidth = 0.5 )
#    m.bluemarble()
#    # =============================================================================
#    # Plot des stations sur le fond de carte  
#    # =============================================================================
##    for j in range (len(lon)) : 
##        xpt,ypt = m(lon[j],lat[j])
##        m.plot(xpt,ypt,'wo', ms = 1) 
##        plt.text(xpt+10000,ypt+10000,'%s' %station[j], color = 'white', fontsize = 5)
##        expt, eypt = m(elon,elat)
##        m.plot(expt,eypt,'r*', ms = 7)
#
#    return m 
#
## =============================================================================
## Plot TEC on map NZ 
## =============================================================================
#def plot_nz(lon, lat, sat, station, lon_sat_seism, lat_sat_seism, type_sat): 
# 
#        plt.title('Trace des satellites {1} au cours du séisme vu par {0}'.format(station, type_sat))
#        m = basic_nz_map()
#        x, y = m(np.degrees(lon), np.degrees(lat))
#        m.plot(x, y, label = "{0}".format(sat))
#        u, v = m(np.degrees(lon_sat_seism), np.degrees(lat_sat_seism))
#        m.plot(u, v, 'w*', ms = 5)
#        plt.legend(bbox_to_anchor=(1.15, 1), loc='upper left', 
#                   borderaxespad=0., fancybox = True, shadow = True)
#        return m

## =============================================================================
## Graphes TECs
## =============================================================================
## Varion TEC and VTEC from RINEX fonction
#def tec(d, courbe_1, courbe_2, sat, station, type_sat) : 
#    plt.xlabel('Temps IUT (h)')
#    plt.ylabel('TEC (TECU)')
#    plt.title('TEC and VTEC measured by {2}: {0} and recorded by : {1}'.format(sat,station,type_sat))
#    plt.plot(d, courbe_1, label = "{0}".format('tec'), color = 'green')
#    plt.plot(d, courbe_2, label = "{0}".format('vtec'), color = 'blue')
#    plt.axvline(x = 11.033, color = 'red', label = "seism", linestyle = ':')
#    delta_tec = plt.legend()
#    plt.gcf()
#    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/{0}_{1}_tec.png'.format(sat,station))
#    plt.show() 
#    return delta_tec
#
## Varion TEC and VTEC from RINEX fonction for uniq_data
#def tec_uniq(d, courbe_1, courbe_2) : 
#    plt.xlabel('Temps IUT (h)')
#    plt.ylabel('TEC (TECU)')
#    plt.title('TEC and VTEC')
#    plt.plot(d, courbe_1, label = "tec", color = 'green')
#    plt.plot(d, courbe_2, label = "vtec", color = 'blue')
#    plt.axvline(x = 11.033, color = 'red', label = "seism", linestyle = ':')
#    delta_tec = plt.legend()
#    plt.gcf()
#    plt.show() 
#    return delta_tec
#
## Plot TEC d'un satellite pour plusieurs stations 
#def tec_1sat_stations(d, courbe_1, courbe_2,station) : 
#    plt.xlabel('Temps IUT (h)')
#    plt.ylabel('TEC (TECU)')
#    plt.title('TEC and VTEC measured by {0} and recorded by G20'.format(station))
#    plt.plot(d, courbe_1, label = "{0}".format('tec'), color = 'green')
#    plt.plot(d, courbe_2, label = "{0}".format('vtec'), color = 'blue')
#    plt.axvline(x = 11.033, color = 'red', label = "seism", linestyle = ':')
#    delta_tec = plt.legend()
#    plt.gcf()
#    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/1sat_stations_tec.png')
#    plt.show() 
#    return delta_tec
#
## Plot Filtered TEC
#def polynomied_tec(time_plot, courbe_1, courbe_2, sat, station, type_sat, window, order) : 
#    # plot du tec polynomial et du vtec - vtecpolynomial 
#    plt.xlabel('Temps IUT (h)')
#    plt.ylabel('TEC (TECU)')
#    plt.title('TEC and VTEC measured by {2}: {0} and recorded by : {1}'.format(sat,station, type_sat))
#    plt.plot(time_plot, courbe_1, label = "vtec_polyn", color = 'red')
#    plt.plot(time_plot, courbe_2, label = "vtec_reduced", color = 'blue')
#    # ligne verticale représentant le séisme 
#    plt.axvline(x = 11.033, color = 'k', label = "seism", linestyle = ':')
#    delta_tec = plt.legend()
#    plt.gcf()
#    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/{0}_{1}_tec_polyn_filtered_compar.png'.format(sat,station))
#    plt.show()
#    plt.plot(time_plot, courbe_2, label = "vtec_filtered", color = 'b')
#    plt.axis([10.5, 12, -0.25, 0.25])
#    plt.xlabel('Temps IUT (h)')
#    plt.ylabel('TEC (TECU)')
#    plt.title('Polynomied VTEC: (window {3}, {4}th order) measured by {2}: {0} recorded by : {1}'.format(sat,station, type_sat,window,order))
#    plt.axvline(x = 11.033, color = 'k', label = "seism", linestyle = ':')
#    delta_tec = plt.legend()
#    plt.gcf()
#    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/{0}_{1}_tec_polyn_filtered.png'.format(sat,station))
#    plt.show()
#    return delta_tec
#
#def polynomied_tec_1sat_stations(time_plot, courbe_2, station) : 
#
#    plt.plot(time_plot, courbe_2, label = "{0}".format(station))
##    labelLines(plt.gca().get_lines(),align = True)
#    plt.axis([10.5, 12, -0.5, 20])
#    plt.xlabel('Temps IUT (h)')
#    plt.ylabel('TEC (TECU)')
#    plt.title('Polynomied VTEC of R07 for stations with latitudes smaller than the epicentre')
#    delta_tec = plt.axvline(x = 11.033, color = 'k', linestyle = ':')
##    delta_tec = plt.legend(bbox_to_anchor=(1, 1), loc='upper left', 
##                   borderaxespad=0., fancybox = True, shadow = True)
##    plt.gcf()
##    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/1sat_stations_bis_tec_polyn_filtered.png')
##    plt.show()
#    return delta_tec





## =============================================================================
##  Animation 
## =============================================================================
## =============================================================================
## Calcul longitude et lattitude de la station pour module animation
## =============================================================================
#def lecture_lat_lon_sat_anim(name_dir, files): 
#        with open(name_dir + '/' + files[14],'r') as fich:
#            for j in range(7):
#                ligne = fich.readline()
#        liste_mots = ligne.split(":")
#        mot=liste_mots[1]
#        liste_bis_mots = mot.split(",")
#        rlon = np.radians(float(liste_bis_mots[0])) #longitude de la station
#        rlat = np.radians(float(liste_bis_mots[1])) #lattitude de la station
#        return(rlon, rlat)
## =============================================================================
##Fonction plot Nouvelle-Zélande pour l'animation
## =============================================================================
#def plot_nz_anim(lon, lat, station, lon_sat_seism, lat_sat_seism, type_sat) : 
# 
#        
#        bmap = basic_nz_map()
#        plt.title('Trace du R07')
#        x, y = bmap(np.degrees(lon), np.degrees(lat))
#        
##        m.plot(x, y)
##        u, v = m(np.degrees(lon_sat_seism), np.degrees(lat_sat_seism))
##        m.plot(u, v, 'w*', ms = 5)
##        plt.legend(bbox_to_anchor=(1.3, 1), loc='upper left', 
##                   borderaxespad=0., fancybox = True, shadow = True)
#        
#        return bmap
## =============================================================================
## Création du fond de carte de NZ pour l'animation
## =============================================================================
#def basic_nz_map_anim() : 
#    
#    # =============================================================================
#    # Lecture du fichier en deux fois pour problème de type
#    # =============================================================================
#    datat = np.loadtxt('/Users/antoineleblevec/Desktop/2016_seism/coord_stations.txt', usecols = (1,2,3,4,5))
#    data = open('/Users/antoineleblevec/Desktop/2016_seism/coord_stations.txt', "r")
#    # =============================================================================
#    # Stockage des noms des stations dans : station 
#    # =============================================================================
#    contenu = data.readlines()
#    station = []
#    for j in range (36) : 
#        list = contenu[j].split()
#        station.append(list [0])
#    # =============================================================================
#    # Stockage des latitudes et longitudes des stations dans lon et lat
#    # =============================================================================
#    lon = []
#    lat = []
#    for i in datat : 
#        lon.append(i[4])
#        lat.append(i[3])
#    # =============================================================================
#    # Epicentre longitude latitude 
#    # =============================================================================
#    elon = 173.054
#    elat = -42.737
#    # =============================================================================
#    # Plot du fond de carte de nouvelle-zélande
#    # =============================================================================
#    m = Basemap(width=1500000,height=1500000,projection='lcc',
#                    resolution='c',lat_1=-80.,lat_2=-20,lat_0=-40,lon_0=176.)
#    m.drawcoastlines()
#    m.drawstates()
#    m.drawcountries()
#    parallels = np.arange(-60, -30, 5.)
#    m.drawparallels(parallels, labels = [False,True,True,False], dashes = [1,5], linewidth = 0.5)
#    meridians = np.arange(160., 190., 5.)
#    m.drawmeridians(meridians, labels = [True, False, False, True], dashes = [1,5], linewidth = 0.5 )
#    m.bluemarble()
#    # =============================================================================
#    # Plot des stations sur le fond de carte  
#    # =============================================================================
#    for j in range (len(lon)) : 
#        xpt,ypt = m(lon[j],lat[j])
#        m.plot(xpt,ypt,'wo', ms = 1) 
#        plt.text(xpt+10000,ypt+10000,'%s' %station[j], color = 'white', fontsize = 5)
#        expt, eypt = m(elon,elat)
#        m.plot(expt,eypt,'r*', ms = 7)
#
#    return m 


# =============================================================================
# Varion TEC and VTEC from RINEX fonction for animation
# =============================================================================
#def tec_anim(d, courbe_1, courbe_2) : 
#    plt.xlabel('Temps IUT (h)')
#    plt.ylabel('TEC (TECU)')
#    plt.title('TEC and VTEC measured by Glonass R21')
#    plt.plot(d, courbe_1, label = "tec", color = 'green')
#    plt.plot(d, courbe_2, label = "vtec", color = 'blue')
#    plt.axvline(x = 11.033, color = 'red', label = "seism", linestyle = ':')
#    delta_tec = plt.legend()
##    plt.gcf()
##    plt.savefig('/Users/antoineleblevec/Desktop/2016_seism/anim_tec.png')
##    plt.show() 
#    return delta_tec
#
#def lecture_lat_lon_sat_anim(name_dir, files): 
#        with open(name_dir + '/' + files[14],'r') as fich:
#            for j in range(7):
#                ligne = fich.readline()
#        liste_mots = ligne.split(":")
#        mot=liste_mots[1]
#        liste_bis_mots = mot.split(",")
#        rlon = np.radians(float(liste_bis_mots[0])) #longitude de la station
#        rlat = np.radians(float(liste_bis_mots[1])) #lattitude de la station
#        return(rlon, rlat)
