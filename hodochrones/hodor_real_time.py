#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 13:13:13 2019

@author: antoineleblevec
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fonctions as f1
from functools import reduce
#from statistics import mean 

Re = 6371032 
H = 250e3
lllat = 33.7; urlat = 43.7; lllon = 133.6; urlon = 150.6
elon = 142 
elat = 38 
epoque1 = 21183
epoque2 = 21883
station = []; df = []; lon_station = []; lat_station = []; lat = []; lon= [];
lon_sip_max = []; lat_sip_max = []; tec = []; vtec=[]; saq = []
df_el = pd.DataFrame(); df_az = pd.DataFrame(); df_tec = pd.DataFrame()

directory = os.path.join('/Users/antoineleblevec/Desktop/G26')
rep = os.path.abspath(os.path.expanduser(directory))
files = os.listdir(rep) 
files.sort()
os.chdir(directory)

for i in range(len(files)):
    name = files[i].split('_')
    station.append(name[0])
    df.append(f1.read(files[i],epoque1,epoque2))
    df[i] = df[i].set_index("tsn")
    df[i].columns = ['el_{0}'.format(station[i]), 'az_{0}'.format(station[i]), 'tec_{0}'.format(station[i])]
    lon_station.append(f1.lecture_lon(files[i]))
    lat_station.append(f1.lecture_lat(files[i]))
dftot = reduce(lambda x, y: pd.merge(x, y, on = "tsn"), df)
df_lat_station = pd.DataFrame([lat_station])
df_lon_station = pd.DataFrame([lon_station])
for j in range(len(files)):
    df_el = df_el.append(dftot['el_{0}'.format(station[j])])
    df_az = df_az.append(dftot['az_{0}'.format(station[j])])
    df_tec = df_tec.append(dftot['tec_{0}'.format(station[j])])
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

## détection de l'indice, en partant de 0, de la première apparition de l'onde en regardant le VTEC
def detec_tid(station): 
    a = df_vtec['tec_{0}'.format(station)]
    for i in range(len(a)): 
        if abs(a.iloc[i+1]-a.iloc[i]) > 0.030 : 
            break
    return i

def devtec(station, i): 
    return df_vtec['tec_{0}'.format(station)].iloc[i]

def delon(station, i): 
    return df_lon['el_{0}'.format(station)].iloc[i]

def delat(station, i):
    return df_lat['el_{0}'.format(station)].iloc[i]
#

for i in range (100):  
    fig = plt.figure()
    for k in station :        
        vtec.append(devtec(k, detec_tid(k)+i))
        lon_sip_max.append(delon(k, detec_tid(k)+i))
        lat_sip_max.append(delat(k, detec_tid(k)+i))

    m = f1.basic_japan_map(lllat, urlat, lllon, urlon, elon, elat) 
    x, y = m(np.degrees(lon_sip_max), np.degrees(lat_sip_max))
    m.hexbin(x,
             y,
             C=vtec,
             reduce_C_function=np.mean,
             gridsize=20, 
             cmap="plasma")
    plt.title('Max VTEC with Hion:{0} m {1}s après earthquake'.format(H,i+400))
    cbaxes = fig.add_axes([0.90, 0.1, 0.01, 0.8]) 
    cb = plt.colorbar(cax=cbaxes)
    m.colorbar()
    plt.gcf()
    fig.savefig(f"/Users/antoineleblevec/Desktop/frames/realtime_{i:04d}.png", 
                frameon=False, pad_inches=0)
    cbaxes.clear()


# =============================================================================
# à utiliser si points aberrants
# =============================================================================
#for i,a in enumerate(saq) : 
#    if abs(mean(saq)-a) > 30 : 
#        saq.pop(i)
#        lon_sip_max.pop(i)
#        lat_sip_max.pop(i)
#        station.pop(i)
#        vtec.pop(i)
