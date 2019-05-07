#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Apr 17 15:18:58 2019

@author: antoineleblevec

Changer de répertoire en fonction du set de données étudié


Pour séisme de Tohoku, le séisme a eu lieu à l'époque 20 783 ; donc prendre entre 21 183 et 21 883 
Pour séisme de Kaikura, le séisme a eu lieu à l'époque  39 776 ; donc prendre entre 39 076 et 43 199

tos : time of seism
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fonctions as f1
import sys
import timeit 
#import math
from functools import reduce
#from statistics import mean 

start = timeit.default_timer()

# =============================================================================
# Définition des constantes et variables utilisées
# =============================================================================
Re = 6371032 
H = 250e3
# latitude et longitude pour la carte désirée 
lllat = 33.7; urlat = 43.7; lllon = 133.6; urlon = 150.6
# epicentre du séisme 
elon = 142 ; elat = 38 
# period of observation : à changer en fonction du séisme étudié 
propagation_period = 100 
tos = 20783
epoque1 = 21183
epoque2 = 21883
# seuil de detection
valeur_inf_seuil = 0.031
valeur_sup_seuil = 0.033
# satelite choisi
sat = 'G26'
# Listes et dataframes
station = []; df = []; lon_station = []; lat_station = [];
lon_sip_max = []; lat_sip_max = []; vtec=[]; saq = [];
df_el = pd.DataFrame(); df_az = pd.DataFrame(); df_tec = pd.DataFrame()

# =============================================================================
# Accès et lecture des données
# =============================================================================
# chemin vers données
directory = os.path.join('/Users/antoineleblevec/Desktop/tohoku_1s/tohoku_1s_{0}'.format(sat))
rep = os.path.abspath(os.path.expanduser(directory))
files = os.listdir(rep) 
files.sort()
os.chdir(directory)

# lecture des données et stockage des latitudes et longitudes des stations présentes
for file in files:
    name = file.split('_')
    da_fr = f1.read(file,epoque1,epoque2)
    if len(da_fr) < 10 :
        continue 
    else :
        df.append(da_fr)
        station.append(name[0])
        lon_station.append(f1.lecture_lon(file))
        lat_station.append(f1.lecture_lat(file))

df_lat_station = pd.DataFrame([lat_station])
df_lon_station = pd.DataFrame([lon_station])
df_lat_station.columns = station
df_lon_station.columns = station

if len(df) == 0 : 
        sys.exit("mauvais sat")

# merging de tous les dataframes et création des listes de longitude et latitude des stations
for i,stat in enumerate(station):
    df[i] = df[i].set_index("tsn")
    df[i].columns = ['el_{0}'.format(stat), 'az_{0}'.format(stat), 'tec_{0}'.format(stat)]

dftot = reduce(lambda x, y: pd.merge(x,y,on="tsn",how='outer'),df)

# extraction des trois tableaux permettant de calculer le tec
for j in station: 
    df_el = df_el.append(dftot['el_{0}'.format(j)])
    df_az = df_az.append(dftot['az_{0}'.format(j)])
    df_tec = df_tec.append(dftot['tec_{0}'.format(j)])

# =============================================================================
# Calcul du TEC à partir des dataframes créés
# =============================================================================
df_el = np.radians(df_el).T
df_az = np.radians(df_az).T
df_tec = df_tec.T
df_tec = df_tec.replace(0, np.NaN)
df_tec = df_tec.dropna(axis=1)
df_tec = df_tec - df_tec.iloc[0]
col_names_tec = list(df_tec.columns.values)
col_names_el = [x.replace('tec','el') for x in col_names_tec]
df_el = df_el[col_names_el]
df_x = np.arcsin((Re * np.cos(df_el)) / (Re + H))
df_x = df_x[col_names_el]
col_names_az = [x.replace('el','az') for x in col_names_el]
df_az = df_az[col_names_az] 
col_names_station = [x[3:] for x in col_names_az]
df_lat_station = df_lat_station[col_names_station]
df_lon_station = df_lon_station[col_names_station]
df_ksi = (np.pi / 2) - (df_el.add((df_x), fill_value=0))
tplat = pd.np.multiply(np.sin(df_ksi),np.cos(df_az))
tplon = pd.np.multiply(np.sin(df_ksi),np.sin(df_az))
df_lat = np.arcsin(pd.np.multiply(np.cos(df_ksi),np.sin(df_lat_station)) + 
                   pd.np.multiply(tplat,np.cos(df_lat_station)))
tpplon = np.arcsin(pd.np.divide(tplon,np.cos(df_lat)))
df_lon = pd.np.add(tpplon,df_lon_station)
df_x = df_x.dropna()
df_vtec = pd.np.multiply(df_tec,np.cos(df_x))

# =============================================================================
# Fonctions pour analyse des données
# =============================================================================
# détection de l'indice, commencant à 0, de la première apparition de l'onde sismique à partir du VTEC
def devtec(station): 
    a = df_vtec['tec_{0}'.format(station)]
    for i in range(len(a)-1): 
        if abs(a.iloc[i+1]-a.iloc[i]) > valeur_inf_seuil and abs(a.iloc[i+1]-a.iloc[i]) < valeur_sup_seuil :
            break
    return i

# retourne l'indice tsn du maximum de vtec 100s après la détection de l'onde sismique
def ivtecmax(station):
    a = df_vtec['tec_{0}'.format(station)]
    a = a[devtec(station):devtec(station) + propagation_period]
    return a.idxmax()

# retourne la valeur du maximum de vtec
def devtecmax(station):
    a = df_vtec['tec_{0}'.format(station)]
    a = a[devtec(station):devtec(station) + propagation_period]
    return a.max()

# retourne la longitude et latitude correspondant à ce maximum de vtec
def delon(station): 
    return df_lon['el_{0}'.format(station)].iloc[ivtecmax(station) - (epoque1 + 1)]
def delat(station):
    return df_lat['el_{0}'.format(station)].iloc[ivtecmax(station) - (epoque1 + 1)]

# =============================================================================
# Plot des données
# =============================================================================
for k in col_names_station : 
    tod = devtec(k)+ (epoque1 - tos + 1)  
    if tod < 700: 
        saq.append(tod)
        vtec.append(devtecmax(k))
        lon_sip_max.append(delon(k))
        lat_sip_max.append(delat(k))

if not saq: 
    print ("pas de variation assez grande dans le tec pour détecter une onde")

# plot map of max vtec
if saq : 
    fig = plt.figure()
    m = f1.basic_japan_map(lllat, urlat, lllon, urlon, elon, elat) 
    x, y = m(np.degrees(lon_sip_max), np.degrees(lat_sip_max))
    m.hexbin(x,
             y,
             C=vtec,
             reduce_C_function=np.mean,
             gridsize=20, 
             cmap="plasma")
    plt.title('Max VTEC with Hion:{0}m ; {1}<seuil<{2} ; sat:{3}'.format(H, valeur_inf_seuil, valeur_sup_seuil, sat))
    cbaxes = fig.add_axes([0.90, 0.1, 0.01, 0.8]) 
    cb = plt.colorbar(cax=cbaxes)
    m.colorbar()
    plt.gcf()
    plt.show()
#
## paramètres généraux des données 
#df_param = pd.DataFrame({
#                'saq' : saq,
##                'GPS Site' : col_names_station,
#                'VTEC max': vtec,
#                'Lon of SIP max' : np.degrees(lon_sip_max), 
#                'Lat of SIP max' : np.degrees(lat_sip_max)
#                })
#print(df_param.describe())

# =============================================================================
# Graphe des vtec de toutes les stations pour le satellite
# =============================================================================
#df_vtec.plot()
#plt.show()

# =============================================================================
# À utiliser pour ploter les stations ; rajouter dans fonctions.py lon_station, lat_station, station
# =============================================================================
#lon_station = np.degrees(lon_station)
#lat_station = np.degrees(lat_station)
#fig = plt.figure()
#m = f1.basic_japan_map(lllat, urlat, lllon, urlon, elon, elat, lon_station, lat_station, station)
#plt.show()

stop = timeit.default_timer()
print ('Time: ', stop-start)