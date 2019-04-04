#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 20:54:47 2019

@author: antoineleblevec
"""

import os
import pandas as pd
from functools import reduce 
import numpy as np
import fonction_1 as f1  
from scipy.signal import savgol_filter
from matplotlib.font_manager import FontProperties
from mpl_toolkits.basemap import Basemap
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show

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
# Lecture du dossier avec tous les fichiers 
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
    
    df = df[["tsn",  "el", "az", "tec"]]
    
    df = df[1320 < df["tsn"]]
    return df[df["tsn"] < 1393]


# =============================================================================
# Concaténation de toutes les stations dans un grand dataframe
# =============================================================================
df = []
dfh = []

for i in range(len(files)) : 
    df.append(read(files[i]))
for j in range (len(df)) : 
    df[j] = df[j].set_index("tsn")
dftot = reduce(lambda x, y: pd.merge(x, y, on = "tsn"), df)
#print (dftot)

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
df_rlon.columns = ['ksi']*35
#print(df_rlon)

# =============================================================================
# Calcul du tec réduit, de la latitude et de la longitude 
# =============================================================================
df_el = dftot.el
df_elx = dftot.el_x
df_ely = dftot.el_y
df_ele = df_elx.add(df_ely, fill_value = 0)
df_ele['el'] = pd.Series(df_el)
df_ele = df_ele.rename(index=str, columns={"el_x" : "ele", "el_y" : "ele", "el" : "ele"})
df_ele = np.radians(df_ele)

df_az = dftot.az
df_azx = dftot.az_x
df_azy = dftot.az_y
df_azi = df_azx.add(df_azy, fill_value = 0)
df_azi['az'] = pd.Series(df_az)
df_azi = df_azi.rename(index=str, columns={"az_x" : "azi", "az_y" : "azi", "az" : "azi"})
df_azi = np.radians(df_azi)

df_x = np.arcsin((Re*np.cos(df_ele))/(Re + H))
df_x = df_x.rename(index=str, columns={"ele" : "x"})

df_x = df_x.rename(index=str, columns={"x" : "ele"})
df_ksi = np.pi/2 - (df_ele.add((df_x), fill_value = 0))
df_ksi = df_ksi.rename(index=str, columns={'ele':'ksi'})

# étapes de calcul pour la lat et lon
tplat = pd.np.multiply(np.sin(df_ksi),np.cos(df_azi))
tplon = pd.np.multiply(np.sin(df_ksi),np.sin(df_azi))

#lat = np.arcsin(np.sin(rlat)*np.cos(ksi) + np.cos(rlat)*np.sin(ksi)*np.cos(az))
#lon = rlon + np.arcsin(np.sin(ksi)*np.sin(az)/np.cos(lat))

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

def minimum(x) : 
    return (x.min())

df_tecr = df_tecr - df_tecr.apply(minimum,axis=0)

df_vtec = pd.np.multiply(df_tecr,np.cos(df_x))
#print (df_vtec)
df_vtecf = savgol_filter(df_vtec,window,order)
#print (df_vtecf)
df_vtecr = df_vtec - df_vtecf
#print (df_vtecr)

# =============================================================================
# création de tableau par époques 
# =============================================================================

#print (df_lon.loc['1321'])
df_lon.reset_index(drop = True, inplace = True)
df_lat.reset_index(drop = True, inplace = True)
df_vtecr.reset_index(drop = True, inplace = True)

#rename columns
df_lon.columns = stations
df_lat.columns = stations
df_vtecr.columns = stations

#print (type(df_vtecr.iloc[:1]))
#print (type(df_lat.iloc[:1]))

#df0 = pd.concat([df_lon.iloc[:1],df_lat.iloc[:1],df_vtecr.iloc[:1]], axis = 0, keys =["lon","lat","tec"])
#data = pd.concat([df_lon.iloc[:1].values.flatten(),df_lat.iloc[:1].values.flatten(),df_vtecr.iloc[:1].values.flatten()])
data = {'lon' : df_lon.iloc[:1].values.flatten(),
        'lat' :  df_lat.iloc[:1].values.flatten(),
        'tec': df_vtecr.iloc[:1].values.flatten()}
frame = pd.DataFrame(data)
#df0 = pd([df_lon.iloc[:1].values.flatten(),df_lat.iloc[:1].values.flatten(),df_vtecr.iloc[:1].values.flatten()],
#                  index =["lon","lat","tec"])

#print (df0)
#df0_transposed = df0.T
#print(df0_transposed)
#print(df0_transposed.dtypes)
#
#m = f1.basic_nz_map()
ax = frame.plot.hexbin(x='lon',
                       y='lat',
                       C='tec',
                       reduce_C_function=np.mean,
                       gridsize=20, 
                       cmap="viridis")
#
#ax = m.hexbin(x=frame['lon'],
#              y=frame['lat'],
#              C=frame['tec'],
#              reduce_C_function=np.mean,
#              gridsize=20)
#np.plot(ax)









