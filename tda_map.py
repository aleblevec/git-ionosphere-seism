#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 14:41:45 2019
@author: antoineleblevec

Pour séisme de Tohoku, le séisme a eu lieu à l'époque 20 783 ; donc prendre entre 21 183 et 21 883 
Pour séisme de Kaikura, le séisme a eu lieu à l'époque  39 776 ; donc prendre entre 39 076 et 43 199
East Cost of Honshu / foreshock de Tohoku : séisme a lieu à l'époque 9920 : e1 = 10200 et e2 = 10900 ; coordonnées de l'épicentre: 38.435 et 142.842
les ondes sont censées arrivées à environ 420s après le séisme. 

tos : time of seism
p_p : propagation_period temps pendant lequel on recherche le maximum de TEC après la première détection de l'onde acoustisque (?)
pour l'instant l'algorithme supprime toutes les stations qui ont enregistré un 0 ou un NA
e1 : epoque1
e2 : epoque2 
v_i_s: valeur_inf_seuil
v_s_s : valeur_sup_seuil 

col_names_station,lo_station,la_station : à utiliser dans l'appel de japan_map pour plotter les stations 
lo_station = df_lon_station.iloc[0]
la_station = df_lat_station.iloc[0]
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
H = 187e3
# latitude et longitude pour la carte désirée 
lllat = 31.7; urlat = 45.7; lllon = 131.6; urlon = 152.6
# epicentre du séisme 
elon = 142 ; elat = 38 
# period of observation : à changer en fonction du séisme étudié 
p_p = 100 
tos = 9920
e1 = 10200
e2 = 10900
# seuil de detection
v_i_s = 0.03
v_s_s = 0.04    
# satelite choisi
sat = 'G10'
# Listes et dataframes
station = []; df = []; lon_station = []; lat_station = [];
lon_sip_max = []; lat_sip_max = []; vtec=[]; saq = [];
df_el = pd.DataFrame(); df_az = pd.DataFrame(); df_tec = pd.DataFrame()
vtec_sip_tda = []; lon_sip_tda = []; lat_sip_tda = [] ; tda = []

# =============================================================================
# Accès, lecture des données et création des DataFrames
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
    da_fr = f1.read(file,e1,e2)
    if len(da_fr) < 10 :  # données incomplètes pour cette station
        continue 
    else :
        df.append(da_fr)
        station.append(name[0])
        lon_station.append(f1.lecture_lat_lon(file)[1])
        lat_station.append(f1.lecture_lat_lon(file)[0])

df_lat_station = pd.DataFrame([lat_station])
df_lon_station = pd.DataFrame([lon_station])
df_lat_station.columns = station
df_lon_station.columns = station

if len(df) == 0 :   # satellite inutilisable
        sys.exit("aucunes données correspondent")

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
# Calcul du TEC
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
# Plot des données
# =============================================================================
for k in col_names_station : 
    tod = f1.tda(k,df_vtec,v_i_s,v_s_s)) + (e1-tos+1)
    tda.append(tod)

print (df_vtec.iloc[0])




stop = timeit.default_timer()
print ('Time: ', stop-start)











