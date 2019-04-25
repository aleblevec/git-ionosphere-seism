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
from functools import reduce
#from statistics import mean 

Re = 6371032 
H = 250e3
#lllat = 33.7; urlat = 43.7; lllon = 133.6; urlon = 150.6
#elon = 142 
#elat = 38 
pmt = 100 
tos = 39776
epoque1 = 40400
epoque2 = 41100
station = []; df = []; lon_station = []; lat_station = []; lat = []; lon= [];
lon_sip_max = []; lat_sip_max = []; tec = []; vtec=[]; saq = []
df_el = pd.DataFrame(); df_az = pd.DataFrame(); df_tec = pd.DataFrame()

directory = os.path.join('/Users/antoineleblevec/Desktop/kaikura_G20_1s_bis')
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



# détection de l'indice, en partant de 0, de la première apparition de l'onde en regardant le VTEC
def devtec(station): 
    a = df_vtec['tec_{0}'.format(station)]
    for i in range(len(a)-1): 
        if abs(a.iloc[i+1]-a.iloc[i]) > 0.025:
            break
    return i


## retourne le vrai indice du maximum de vtec à partir du moment de l'apparition 
## de l'onde jusqu'à 100s après 
def ivtecmax(station):
    a = df_vtec['tec_{0}'.format(station)]
    a = a[devtec(station):devtec(station)+ pmt]
    return a.idxmax()


##retourne la valeur du maximum de vtec
def devtecmax(station):
    a = df_vtec['tec_{0}'.format(station)]
    a = a[devtec(station):devtec(station)+ pmt]
    return a.max()

def delon(station): 
    return df_lon['el_{0}'.format(station)].iloc[ivtecmax(station) - (epoque1 + 1)]

def delat(station):
    return df_lat['el_{0}'.format(station)].iloc[ivtecmax(station) - (epoque1 + 1)]

for k in station : 
    saq.append(devtec(k)+ (epoque1 - tos + 1))
    vtec.append(devtecmax(k))
    lon_sip_max.append(delon(k))
    lat_sip_max.append(delat(k))

df_param = pd.DataFrame({
                        'saq' : saq,
                        'GPS Site' : station,
                        'VTEC max': vtec,
                        'Lon of SIP max' : np.degrees(lon_sip_max), 
                        'Lat of SIP max' : np.degrees(lat_sip_max)
                        })
print(df_param.describe())
#
#fig = plt.figure()
#m = f1.basic_japan_map(lllat, urlat, lllon, urlon, elon, elat) 
#x, y = m(np.degrees(lon_sip_max), np.degrees(lat_sip_max))
#m.hexbin(x,
#         y,
#         C=vtec,
#         reduce_C_function=np.mean,
#         gridsize=20, 
#         cmap="plasma")
#plt.title('Max VTEC with Hion:{0} m avec seuil à 0.03'.format(H))
#cbaxes = fig.add_axes([0.90, 0.1, 0.01, 0.8]) 
#cb = plt.colorbar(cax=cbaxes)
#m.colorbar()
#plt.gcf()
#plt.show()


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