#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 10:03:13 2019

@author: antoineleblevec
"""
import os
from os import mkdir 
import shutil

year = '2011' 
jday = '070' 
seism = 'tohoku'
sate = []

directory = os.path.join('/Users/antoineleblevec/Desktop/tec-suite-master_2/tec/{0}/{1}'.format(year,jday))
rep = os.path.abspath(os.path.expanduser(directory))
stations = os.listdir(rep) 
stations.sort()
os.chdir(directory)

for k in stations: 
    directory1 = os.path.join('/Users/antoineleblevec/Desktop/tec-suite-master_2/tec/{0}/{1}/{2}'.format(year,jday,k))
    rep1 = os.path.abspath(os.path.expanduser(directory1))
    sat_station = os.listdir(rep1) 
    sat_station.sort()
    os.chdir(directory1)
    for i in sat_station : 
        name = i.split('_')
        if name[1] not in sate :
            sate.append(name[1])
            mkdir ('/Users/antoineleblevec/Desktop/tohoku_1s/tohoku_1s_{0}'.format(name[1]))
            shutil.copy(i, '/Users/antoineleblevec/Desktop/tohoku_1s/tohoku_1s_{0}'.format(name[1]))
        if name[1] in sate : 
            shutil.copy(i, '/Users/antoineleblevec/Desktop/tohoku_1s/{0}_1s_{1}'.format(seism,name[1]))
            


        
    