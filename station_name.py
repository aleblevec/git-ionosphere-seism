#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 14:37:22 2019

@author: antoineleblevec
"""

data = open('/Users/antoineleblevec/Desktop/2016_seism/coord_stations.txt', "r")
contenu = data.readlines()
station = []
for j in range (36) : 
    list = contenu[j].split()
    station.append(list [0])
print (station)
