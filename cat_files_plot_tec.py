#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 13:50:05 2019

@author: antoineleblevec
"""
import numpy as np 
import matplotlib.pyplot as plt
from glob import glob

##Lecture d'un fichier concaténé mais mauvaix choix 

 
fnames = glob('2011/068/0047/0047_G*_068_11.dat')
fnames.sort()
arrays = [np.loadtxt(f) for f in fnames]
final_array = np.concatenate(arrays)








#Temps de mesure de tous les satellites 
a = final_array[:,1]     #temps d'enregistrement des données de tous les satellites
d = a [(2 < a) & (a < 4)]   #temps d'enregistrement des satellites au moment du seisme

#TEC calculé par tous les satellites à tous les moments 
b = final_array [:,4]
#c = b - min(final_array[:,4]) # peut être pas un bon plan de soustraire le min pour des mesures différentes 

e = []
f = []
tec = [ ]

########## Passage à un autre satellite ########
for i in range(0,len(d)-1) : #je parcours toute ma liste de données des temps au moment du séisme   
    if abs(d[i+1] - d[i]) < 0.0084 : # tant que l'écart n'est pas supérieure à cette valeur alors c'est toujours 
                                     #le même satellite et je stocke  le tec et le temps dans deux listes 
        f.append(d[i])
        e.append(b[i])
        i = i + 1
    else :                           # quand l'ecart est supérieur alors on passe à un enregistrement fait 
                                     # par un autre satellite ; je plot donc mes premières données correspondant à
                                     # mon premier satellite je calcul le tec en soustrayant le min de ma liste
        
        f.append(d[i])
        e.append(b[i])
        tec = e - min(e)             
        print (tec)
        #print(f)
        plt.xlabel('Temps IUT (h)')
        plt.ylabel('TEC (TECU)')
        plt.title('Variation of TEC during seism for sat:')
        plt.plot(f,tec)
        plt.show() 
        e = []
        f = []
        tec = []
        i = i + 1 

    
        
        


