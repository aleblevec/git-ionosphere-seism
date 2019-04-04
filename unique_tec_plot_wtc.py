#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:00:53 2019

@author: antoineleblevec
"""

import numpy as np 
import matplotlib.pyplot as plt
import os.path 


############Definition des constantes #######################
Re = 6371032 
H = 350e3

while "Bon choix du satellite" : 
    ########Interface utilisateur#####################
    year = input("enter the two last digits of the year __ chosen : ")
    day = input("enter the three digits of the julian day chosen from 0 to 365 : ")
    station = input("enter the four number of the station chosen : ")
    sat = input("enter the name of the satellite G__ :")
    print (" \n Choisir les bornes sup et inf de temps d'observation du satellite :")
    borne_inf = input("borne_inf :")
    borne_sup = input("borne_sup :")
    name = os.path.join('/Users/antoineleblevec/Desktop', year, day, station, station + '_' + sat + '_' + day + '_' + year + '.dat')
    print (name)
    
    data = np.loadtxt(name)
     
    ##########Calcul of TEC################# 
    b = data [:,4]
    #c = b - min(data[:,4]) #il ne faut en fait pas le soustraire ici, car il soustrait tous les TEC même quand 
    #le sattelite n'est plus au moment du séisme
    #print(np.shape(c))
    
    ###########Sélection du moment du séisme ######################
    a = data[:,1]    #liste correspondant aux heures 
    d = a [(int(borne_inf) < a) & (a < int(borne_sup))]
#    d = a [(4 < a) & (a < 6)]   #heures encadrant le séisme 
    #print(d)
    if len(d) == 0  :    ###### sortie de boucle si le sat n'étais pas au dessus du séisme au bon moment
        print("This sat didn't observe the seism")
        print("Launch the program with another satellite")
        break ; 
    #print(np.shape(d))
    
    ############Sélection des données du TEC au moment du séisme ######################
    e = b[(int(borne_inf) < a) & (a < int(borne_sup))]
#    e = b[(4 < a) & (a < 6)]
    ele = data [:,2]    ##### selection en fonction de l'elevation 
    elevation = ele [(int(borne_inf) < a) & (a < int(borne_sup))]
#    elevation = ele [(4 < a) & (a < 6)]
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
    plt.title('Variation of TEC during seism observed by satellite : %s ' %sat)
    plt.plot(d, vtec, label = "vTEC")
    plt.plot(d, tec, label = "Slant TEC")
    plt.legend()
    plt.show()

###################Affichage########################################
    break ; 
