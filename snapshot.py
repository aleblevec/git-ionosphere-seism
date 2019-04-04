#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:08:51 2019

@author: antoineleblevec

amélioration du read_file pour rajouter la carte sur le hexbin.plot 
"""

import os
import pandas as pd
from functools import reduce 
import numpy as np

from scipy.signal import savgol_filter
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show

import fonction_1 as f1  

# =============================================================================
# Constantes 
# =============================================================================
Re = 6371032 
H = 350e3
lllat = 31.7; urlat = 43.7; lllon = 130.6; urlon = 145.6
window = 27
order = 4
stations = ('auck','bluf','chti','dnvk','dund','gisb','gldb','haas','hamt','hast','hikb',
            'hoki','kaik','ktia','lexa','lkta','maho','mavl','meth','mqzg','mtjo','nlsn',
            'nply','pygr','taup','trng','vgmt','waim','wang','wark','west','wgtn','whkt',
            'whng','wrpa')

# =============================================================================
# Lecture du dossier comprenant les différents fichiers 
# =============================================================================
directory = os.path.join('/Users/antoineleblevec/Desktop/G20')
rep = os.path.abspath(os.path.expanduser(directory))
files = os.listdir(rep) 
files.sort()
os.chdir(directory)

# =============================================================================
# Fonction de lecture de mon dataframe 
# =============================================================================
def read(files):
    df = pd.read_csv(
            files,
            skiprows=11,
            delim_whitespace=True,
            header=None,
            names=["tsn", "hour", "el", "az", "tec", "tec.p1p2", "validity"]
            )
    
    df = df[["tsn", "el", "az", "tec"]]
    
    df = df[1320 < df["tsn"]]
    return df[df["tsn"] < 1393]

# =============================================================================
# Concaténation des fichiers de chaque station dans un dataframe en fonction de 
# l'index correspondant aux époques 
# =============================================================================
df = []

for i in range(len(files)) : 
    df.append(read(files[i]))
for j in range (len(df)) : 
    df[j] = df[j].set_index("tsn")
dftot = reduce(lambda x, y: pd.merge(x, y, on = "tsn"), df)

# =============================================================================
# Lecture de la lon et lat de chaque stations
# =============================================================================
rlon = []
rlat = []
for k in range(len(files)) :
    rlon.append(f1.lecture_lon(files[k]))
    rlat.append(f1.lecture_lat(files[k]))
df_rlat = pd.DataFrame([rlat])
df_rlon = pd.DataFrame([rlon])

# =============================================================================
# Calcul du tec réduit, de la latitude et de la longitude 
# =============================================================================
df_el = dftot.el
df_elx = dftot.el_x
df_ely = dftot.el_y
df_ele = df_elx.add(df_ely, fill_value=0)
df_ele['el'] = pd.Series(df_el)
df_ele = df_ele.rename(index=str, columns={"el_x":"ele", "el_y":"ele", "el":"ele"})
df_ele = np.radians(df_ele)

df_az = dftot.az
df_azx = dftot.az_x
df_azy = dftot.az_y
df_azi = df_azx.add(df_azy, fill_value=0)
df_azi['az'] = pd.Series(df_az)
df_azi = df_azi.rename(index=str, columns={"az_x":"azi", "az_y":"azi", "az":"azi"})
df_azi = np.radians(df_azi)

df_x = np.arcsin((Re * np.cos(df_ele)) / (Re + H))
df_x = df_x.rename(index=str, columns={"x" : "ele"})

df_ksi = (np.pi / 2) - (df_ele.add((df_x), fill_value=0))
df_ksi = df_ksi.rename(index=str, columns={'ele':'ksi'})

# étapes de calcul nécessaire à cause des opérations entre dataframe
tplat = pd.np.multiply(np.sin(df_ksi),np.cos(df_azi))
tplon = pd.np.multiply(np.sin(df_ksi),np.sin(df_azi))

df_lat = np.arcsin(pd.np.multiply(np.cos(df_ksi),np.sin(df_rlat)) + 
                   pd.np.multiply(tplat,np.cos(df_rlat)))
df_lat = df_lat.rename(index=str, columns={'ksi':'lat'})

tpplon = np.arcsin(pd.np.divide(tplon,np.cos(df_lat)))
df_lon = pd.np.add(tpplon,df_rlon)
df_lon = df_lon.rename(index=str, columns={'ksi':'lon'})

df_tec = dftot.tec
df_tecx = dftot.tec_x
df_tecy = dftot.tec_y
df_tecr = df_tecx.add(df_tecy, fill_value = 0)
df_tecr['tec'] = pd.Series(df_tec)
df_tecr = df_tecr.rename(index=str, columns={"tec_x" : "tecr", "tec_y" : "tecr", "tec" : "tecr"})

df_tecr = df_tecr - df_tecr.apply(f1.minimum,axis=0)

df_vtec = pd.np.multiply(df_tecr,np.cos(df_x))
df_vtecf = savgol_filter(df_vtec,window,order)
df_vtecr = df_vtec - df_vtecf

# =============================================================================
# création de tableau par époques 
# =============================================================================
# index commence à 0 
df_lon.reset_index(drop = True, inplace = True)
df_lat.reset_index(drop = True, inplace = True)
df_vtecr.reset_index(drop = True, inplace = True)

# rename columns en fonction des stations
df_lon.columns = stations
df_lat.columns = stations
df_vtecr.columns = stations
#print (df_vtecr)

# valeurs tec : -1.5 -> 0.8

for i in range (1, 8):
    a = i * 10
#    print (a)
    lon = np.degrees(df_lon.iloc[a:a+1].values.flatten())
    lat = np.degrees(df_lat.iloc[a:a+1].values.flatten())
    tec = df_vtecr.iloc[a:a+1].values.flatten()
    
    data = {'lon' : lon,
            'lat' : lat,
            'tec': tec}
    frame = pd.DataFrame(data)
    
    plt.subplot(2,4,i)
    m = f1.basic_nz_map()
    x, y = m(lon, lat)
    m.hexbin(x,
             y,
             C=tec,
             reduce_C_function=np.mean,
             gridsize=20, 
             cmap="viridis", 
             vmin=-1.5, vmax=0.8)
    plt.title('tec at epoque {0}'.format(a))

cb = m.colorbar(pad='20%', label='tec')
#cb.set_ticks([-1.5, -1, -0.5, 0, 0.5, 1])
plt.show()








