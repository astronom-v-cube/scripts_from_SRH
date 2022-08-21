#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 02:43:37 2022

@author: sergey_lesovoi
"""

import base2uvw_1224
import base2uvw_612
import base2uvw_36
import numpy as NP
import pylab as PL

hourAngle = NP.deg2rad(0)
declination = 0.903338787600965#NP.deg2rad(53)

antennas36 = NP.linspace(1,129,129,dtype='int')
antennas612 = NP.linspace(1,192,192,dtype='int')
antennas1224 = NP.linspace(1,207,207,dtype='int')

uv36 = []
lambda3 = 3e8/3e9
for antA in range(129):
    for antB in range(129 - antA):
        uv36.append(base2uvw_36.base2uvw(hourAngle,declination,antennas36[antA],antennas36[antB + antA]))

uv612 = []
lambda6 = 3e8/6e9
for antA in range(192):
    for antB in range(192 - antA):
        uv612.append(base2uvw_612.base2uvw(hourAngle,declination,antennas612[antA],antennas612[antB + antA]))

uv1224 = []
lambda12 = 3e8/12e9
for antA in range(207):
    for antB in range(207 - antA):
        uv1224.append(base2uvw_1224.base2uvw(hourAngle,declination,antennas1224[antA],antennas1224[antB + antA]))

uv36 = NP.array(uv36)
uv612 = NP.array(uv612)
uv1224 = NP.array(uv1224)

PL.figure(figsize=(10,10))
PL.scatter(uv36[:,0],uv36[:,1],s=.1,color='red',label='3-6 GHz')
PL.scatter(uv612[:,0],uv612[:,1],s=.1,color='green',label='6-12 GHz')
PL.scatter(uv1224[:,0],uv1224[:,1],s=.1,color='blue',label='12-24 GHz')
PL.xlim(-1000,1000)
PL.ylim(-1000,1000)
PL.xlabel('west-east')
PL.ylabel('south-north')
PL.legend(markerscale=9)

PL.figure(figsize=(10,10))
PL.scatter(uv36[:,0] / lambda3,     uv36[:,1] / lambda3,s=.1,color='red',label='3 GHz')
PL.scatter(uv612[:,0] / lambda6,    uv612[:,1] / lambda6,s=.1,color='green',label='6 GHz')
PL.scatter(uv1224[:,0] / lambda12,  uv1224[:,1] / lambda12,s=.1,color='blue',label='12 GHz')
PL.xlim(-16000,16000)
PL.ylim(-16000,16000)
PL.xlabel('U')
PL.ylabel('V')
PL.legend(markerscale=9)
