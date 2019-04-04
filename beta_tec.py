#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 09:42:24 2019

@author: antoineleblevec
"""

""" Compter le nombre de fichiers de chaque dossier.
    Les liens symboliques sont ignorés.
"""

## Problème :  dans la lecture des satellites, car les noms des satellites enregistrés par la station 
## ne s'incrémentent pas de un. 
## Trouver un moyen de simplement lire tous les fichiers d'un dossier et de stocker leur noms. 
## Ajout du choix du temps d'observation  



############### modules ##########################
import numpy as np
import os.path
from os import scandir
import matplotlib.pyplot as plt


############Definition des constantes #######################
Re = 6371032 
H = 350e3

########################################### Class pour compter le nombre de fichiers dans un dossier #######################
class Counter(object):
    def __init__(self, path):
        if not path:
            raise ValueError('A girl needs a path.')
        self.path = path
        self.files = 0

    def work(self):
        for entry in scandir(self.path):
            if entry.is_dir() and not entry.is_symlink():
                path = os.path.join(self.path, entry.name)
                counter = Counter(path)
                yield from counter.work()
            else:
                self.files += 1

       
        yield self

    def __str__(self):
       return '{} {}'.format(self.path, self.files)
                
#######################################################################                


########Interface utilisateur#####################
year = input("enter the two last digits of the year __ chosen : ")
day = input("enter the three digits of the julian day chosen from 0 to 365 : ")
station = input("enter the four number of the station chosen : ")
name_dir = os.path.join('/Users/antoineleblevec/Desktop', year, day, station)
print (" \n Choisir les bornes sup et inf de temps d'observation du satellite :")
borne_inf = input("borne_inf :")
borne_sup = input("borne_sup :")


 ####### Nombre de fichier de données #########
counter = Counter(name_dir)
total = 0
for cls in counter.work():
    total += cls.files
#print(name_dir, total)

############# Lecture du fichier de data du sat choisi #########
for i in range(2,total + 1) : 
    if i < 10 : 
        sat = 'G0' + str(i)
        #print (sat)
    if i >= 10 : 
        sat = 'G' + str(i)
    name_data = os.path.join('/Users/antoineleblevec/Desktop', year, day, station, station + '_' + sat + '_' + day + '_' + year + '.dat')
    #print(name_data)
    
    
    ####### Lecture du fichier ###########
    data = np.loadtxt(name_data)
        
    ##########Calcul of TEC################# 
    b = data [:,4]
 
    ###########Sélection du moment du séisme ######################
    a = data[:,1]    #liste correspondant aux heures 
    d = a [(int(borne_inf) < a) & (a < int(borne_sup))]   #heures encadrant le séisme 
    ele = data [:,2]    ##### selection en fonction de l'elevation 
    elevation = ele [(int(borne_inf) < a) & (a < int(borne_sup))]
    
    if len(d) == 0 or np.any(elevation < 10) :    
        print("The sat %s didn't observe the seism or didn't have a good elevation" %i)
        continue
    
    
    ############Sélection des données du TEC au moment du séisme ######################
    e = b[(int(borne_inf) < a) & (a < int(borne_sup))]
    #print (e)
    #print(min(e))
    tec = e - min(e)
    elevation  = np.radians(elevation)
    x = np.arcsin((Re*np.cos(elevation))/(Re + H))
    vtec = tec * np.cos(x)
    #print(tec)
 
    
    ############ Plot du TEC au moment du séisme ###################
    plt.xlabel('Temps IUT (h)')
    plt.ylabel('TEC (TECU)')
    plt.title('Variation of TEC during seism observed by satellite : G%s ' %i)
    plt.plot(d, vtec, label = "vTEC")
    plt.plot(d, tec, label = "Slant TEC")
    plt.legend()
    plt.show()
    
    
######### Création des colonnes de sortie pour la période du séisme #########
    tsn = data[:,0]
    tsn = tsn [(int(borne_inf) < a) & (a < int(borne_sup))]


########## Création d'un fichier de sortie ######################

    f = open("{0}_{1}_{2}_{3}_new.dat".format(station, sat, day, year), "w")
    np.savetxt("{0}_{1}_{2}_{3}_new.dat".format(station, sat, day, year), (tsn,d))
    