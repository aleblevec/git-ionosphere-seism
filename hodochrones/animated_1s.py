#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 11:46:22 2019

@author: antoineleblevec
"""
import os
import pandas as pd
 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time 

from scipy.signal import savgol_filter
from functools import reduce

import fonctions as f1  

# =============================================================================
# Constantes 
# =============================================================================
start_time = time.time()
Re = 6371032 
H = 350e3
lllat = 31.7; urlat = 43.7; lllon = 130.6; urlon = 145.6
window = 811
order = 4
epoque1 = 39600
epoque2 = 43199

# =============================================================================
# Lecture du répertoire ; construction d'une liste avec nom des stations
# =============================================================================
directory = os.path.join('/Users/antoineleblevec/Desktop/G20_1s')
rep = os.path.abspath(os.path.expanduser(directory))
files = os.listdir(rep) 
files.sort()
os.chdir(directory)
stations = []
for a in range (len(files)): 
    name = files[a].split('_')
    stations.append(name[0])
# =============================================================================
# Liste de DataFrame où chaque élément correspond aux données d'une station
# =============================================================================
df = []
for i in range(len(files)) : 
    df.append(f1.read(files[i],epoque1,epoque2))
    df[i] = df[i].set_index("tsn")
    df[i].columns = ['el_{0}'.format(stations[i]), 'az_{0}'.format(stations[i]), 'tec_{0}'.format(stations[i])]
dftot = reduce(lambda x, y: pd.merge(x, y, on = "tsn"), df)
#
###### =============================================================================
###### Lon et Lat des Stations 
###### =============================================================================
rlon = []
rlat = []
for k in range(len(files)) :
    rlon.append(f1.lecture_lon(files[k]))
    rlat.append(f1.lecture_lat(files[k]))
df_rlat = pd.DataFrame([rlat])
df_rlon = pd.DataFrame([rlon])
##
###### =============================================================================
###### Calcul du tec réduit, de la latitude et de la longitude 
###### =============================================================================
df_el = pd.DataFrame()
df_az = pd.DataFrame()
df_tec = pd.DataFrame()
for j in range(len(files)):
    df_el = df_el.append(dftot['el_{0}'.format(stations[j])])
    df_az = df_az.append(dftot['az_{0}'.format(stations[j])])
    df_tec = df_tec.append(dftot['tec_{0}'.format(stations[j])])
#
df_el = np.radians(df_el)
df_az = np.radians(df_az)
#
df_ele = df_el.T
df_azi = df_az.T
df_tecr = df_tec.T
#
df_x = np.arcsin((Re * np.cos(df_ele)) / (Re + H))
df_ksi = (np.pi / 2) - (df_ele.add((df_x), fill_value=0))
#
##### étapes intermédiaire de calcul
tplat = pd.np.multiply(np.sin(df_ksi),np.cos(df_azi))
tplon = pd.np.multiply(np.sin(df_ksi),np.sin(df_azi))
#####
df_lat = np.arcsin(pd.np.multiply(np.cos(df_ksi),np.sin(df_rlat)) + 
                   pd.np.multiply(tplat,np.cos(df_rlat)))
##### étape intermédiaire de calcul
tpplon = np.arcsin(pd.np.divide(tplon,np.cos(df_lat)))
#####
df_lon = pd.np.add(tpplon,df_rlon)
#

# soustrait le tec des mesures de tec prises respectivement par chaque station
df_tecr_reduit = df_tecr - df_tecr.apply(f1.minimum,axis=0)
df_vtec = pd.np.multiply(df_tecr_reduit,np.cos(df_x))
df_vtecf = savgol_filter(df_vtec,window,order,axis=0)
df_vtecr = df_vtec - df_vtecf

# =============================================================================
# Plot des tecs pour voir à ce que je dois m'attendre dans les snapshots 
# =============================================================================
#data1 = pd.DataFrame()
#for i in range(10) : 
#    data1 = data1.append(df_vtecr['tec_{0}'.format(stations[i])])
#data1 = data1.T
#data1.plot(subplots=True, sharex=True, figsize=(20,20))
#plt.show()
#
#print(min(df_vtecr.min()))
### min without v47b : -0.0213447610165849
### max withour v47b : 0.019251844721856037
####### =============================================================================
####### Plot
####### =============================================================================
for i in range (1,360):
    a = i * 10
    fig = plt.figure()
    lon = np.degrees(df_lon.iloc[a:a+1].values.flatten())
    lat = np.degrees(df_lat.iloc[a:a+1].values.flatten())
    tec = df_vtecr.iloc[a:a+1].values.flatten()

    m = f1.basic_nz_map()
    x, y = m(lon, lat)
    m.hexbin(x,
             y,
             C=tec,
             reduce_C_function=np.mean,
             gridsize=20, 
             cmap="viridis", 
             vmin=-0.23, vmax=0.32)

#positionnement de la colorbar
    plt.title('tec at epoque {0}'.format(epoque1+a))
    cbaxes = fig.add_axes([0.90, 0.1, 0.01, 0.8]) 
    cb = plt.colorbar(cax=cbaxes)
    m.colorbar()
    plt.gcf()
    fig.savefig(f"/Users/antoineleblevec/Desktop/frames/frame_{i:04d}.png", 
                frameon=False, pad_inches=0)
    cbaxes.clear()

print("--- %s seconds ---" % (time.time() - start_time))
