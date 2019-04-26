#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 15:18:58 2019

@author: antoineleblevec

Changer de répertoire en fonction du set de données étudié


Pour séisme de Tohoku, le séisme a eu lieu à l'époque 20 783 ; donc prendre entre 21 183 et 21 883 
Pour séisme de Kaikura, le séisme a eu lieu à l'époque  39 776 ; donc prendre entre 39 076 et 43 199

pmt : period of researching maximum tec
tos : time of seism
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fonctions as f1
import sys
import timeit 
from functools import reduce
from statistics import mean 

start = timeit.default_timer()

# =============================================================================
# Constantes
# =============================================================================
Re = 6371032 
# assumed 
H = 250e3
# pour basemap
lllat = 33.7; urlat = 43.7; lllon = 133.6; urlon = 150.6
# epicentre
elon = 142 ; elat = 38 
pmt = 100 
# period of observation
tos = 20783
epoque1 = 21183
epoque2 = 21883
valeur_inf_seuil = 0.031
valeur_sup_seuil = 0.033
sat = 'G12'

station = []; station1 = []; df = []; lon_station = []; lat_station = []; lat = []; lon= [];
lon_sip_max = []; lat_sip_max = []; tec = []; vtec=[]; saq = [];
df_el = pd.DataFrame(); df_az = pd.DataFrame(); df_tec = pd.DataFrame()

directory = os.path.join('/Users/antoineleblevec/Desktop/tohoku_1s/tohoku_1s_{0}'.format(sat))
rep = os.path.abspath(os.path.expanduser(directory))
files = os.listdir(rep) 
files.sort()
os.chdir(directory)

for i in range(len(files)):
    name = files[i].split('_')
    a = f1.read(files[i],epoque1,epoque2)
    if len(a) < 10 :  
        continue 
    else :
        df.append(f1.read(files[i],epoque1,epoque2))
        station.append(name[0])
    if len(df) == 0 : 
        sys.exit("mauvais sat")
        
for i in range (len(df)):    
    df[i] = df[i].set_index("tsn")
    df[i].columns = ['el_{0}'.format(station[i]), 'az_{0}'.format(station[i]), 'tec_{0}'.format(station[i])]
    lon_station.append(f1.lecture_lon(files[i]))
    lat_station.append(f1.lecture_lat(files[i]))
    
dftot = reduce(lambda x, y: pd.merge(x, y, on = "tsn"), df)

for j in range(len(df)):
    df_el = df_el.append(dftot['el_{0}'.format(station[j])])
    df_az = df_az.append(dftot['az_{0}'.format(station[j])])
    df_tec = df_tec.append(dftot['tec_{0}'.format(station[j])])
#    
df_lat_station = pd.DataFrame([lat_station])
df_lon_station = pd.DataFrame([lon_station])
df_el = np.radians(df_el).T
df_az = np.radians(df_az).T
df_tec = df_tec.T
df_tec = df_tec - df_tec.iloc[0]
df_x = np.arcsin((Re * np.cos(df_el)) / (Re + H))
df_ksi = (np.pi / 2) - (df_el.add((df_x), fill_value=0))
tplat = pd.np.multiply(np.sin(df_ksi),np.cos(df_az))
tplon = pd.np.multiply(np.sin(df_ksi),np.sin(df_az))
df_lat = np.arcsin(pd.np.multiply(np.cos(df_ksi),np.sin(df_lat_station)) + 
                   pd.np.multiply(tplat,np.cos(df_lat_station)))
tpplon = np.arcsin(pd.np.divide(tplon,np.cos(df_lat)))
df_lon = pd.np.add(tpplon,df_lon_station)
df_vtec = pd.np.multiply(df_tec,np.cos(df_x))
#
## détection de l'indice, en partant de 0, de la première apparition de l'onde en regardant le VTEC
def devtec(station): 
    a = df_vtec['tec_{0}'.format(station)]
    for i in range(len(a)-1): 
        if abs(a.iloc[i+1]-a.iloc[i]) > valeur_inf_seuil and abs(a.iloc[i+1]-a.iloc[i]) < valeur_sup_seuil :
            break
    return i
#
### retourne le vrai indice du maximum de vtec à partir du moment de l'apparition 
### de l'onde jusqu'à 100s après 
def ivtecmax(station):
    a = df_vtec['tec_{0}'.format(station)]
    a = a[devtec(station):devtec(station) + pmt]
    return a.idxmax()
#
##retourne la valeur du maximum de vtec
def devtecmax(station):
    a = df_vtec['tec_{0}'.format(station)]
    a = a[devtec(station):devtec(station) + pmt]
    return a.max()

def delon(station): 
    return df_lon['el_{0}'.format(station)].iloc[ivtecmax(station) - (epoque1 + 1)]

def delat(station):
    return df_lat['el_{0}'.format(station)].iloc[ivtecmax(station) - (epoque1 + 1)]

## =============================================================================
## Plot des maximum de VTEC 
# =============================================================================
for k in station : 
    a = devtec(k)+ (epoque1 - tos + 1)
    if a < 700: 
        saq.append(devtec(k)+ (epoque1 - tos + 1))
        vtec.append(devtecmax(k))
        lon_sip_max.append(delon(k))
        lat_sip_max.append(delat(k))
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

#df_param = pd.DataFrame({
#                'saq' : saq,
#                'GPS Site' : station1,
#                'VTEC max': vtec,
#                'Lon of SIP max' : np.degrees(lon_sip_max), 
#                'Lat of SIP max' : np.degrees(lat_sip_max)
#                })
#print(df_param.describe())


#
## =============================================================================
## Plot des vtec de toutes les stations pour le satellite
## =============================================================================
##df_vtec.plot()
##plt.show()


### =============================================================================
### À utiliser pour ploter les stations ; rajouter dans fonctions.py lon_station, lat_station, station
### =============================================================================
##lon_station = np.degrees(lon_station)
##lat_station = np.degrees(lat_station)
##fig = plt.figure()
##m = f1.basic_japan_map(lllat, urlat, lllon, urlon, elon, elat, lon_station, lat_station, station)
##plt.show()
#
#stop = timeit.default_timer()
#print ('Time: ', stop-start)