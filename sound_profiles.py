#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:25:57 2019

@author: antoineleblevec
"""
import os 
import fonctions as f1
import pandas as pd
from functools import reduce
import matplotlib.pyplot as plt

#import pandas as pd 

directory = os.path.join('/Users/antoineleblevec/Desktop/sound_speed')
rep = os.path.abspath(os.path.expanduser(directory))
files = os.listdir(rep) 
files.sort()
os.chdir(directory)

df_50 = f1.read_sound('modele_tohF50_out.txt')
df_100 = f1.read_sound('modele_tohF100_out.txt')
df_150 = f1.read_sound('modele_tohF150_out.txt')
df_200 = f1.read_sound('modele_tohF200_out.txt')
df_250 = f1.read_sound('modele_tohF250_out.txt')
df_300 = f1.read_sound('modele_tohF300_out.txt')
df_350 = f1.read_sound('modele_tohF350_out.txt')

df_50['H'] = df_50['H'] - 6371
df_100['H'] = df_100['H'] - 6371
df_150['H'] = df_150['H'] - 6371
df_200['H'] = df_200['H'] - 6371
df_250['H'] = df_250['H'] - 6371
df_300['H'] = df_300['H'] - 6371
df_350['H'] = df_350['H'] - 6371

df = [df_50, df_100, df_150, df_200, df_250, df_300, df_350]
dftot = reduce(lambda x, y: pd.merge(x,y,on="H",how='outer'),df)
dftot.columns = ['H','Cs_50','Cs_100','Cs_150', 'Cs_200', 'Cs_250', 'Cs_300', 'Cs_350']
#dftot = dftot.set_index('H')
Cs = ['Cs_50','Cs_100','Cs_150', 'Cs_200', 'Cs_250', 'Cs_300', 'Cs_350']
fig,ax = plt.subplots()
for i in Cs:
    dftot.plot(
            ax=ax,
            y='H',
            x=i)
    ax.legend(Cs)
    ax.set(xlabel="Cs", ylabel="Hauteur")
plt.show()











