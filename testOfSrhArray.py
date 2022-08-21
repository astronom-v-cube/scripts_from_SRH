#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  2 09:44:52 2018

@author: sergey
"""

import srhArray
import numpy as NP
import pylab as PL

N = 100
hourAngle = NP.linspace(NP.deg2rad(-60),NP.deg2rad(60),N)
declination = NP.deg2rad(0)

uvwEW = NP.zeros((3,N))
uvwNS = NP.zeros((3,N))
uvwES = NP.zeros((3,N))
uvwWS = NP.zeros((3,N))

SRH = srhArray.SrhArray()

for i in range(N):
    uvwEW[:,i] = SRH.baseline2uvw(hourAngle[i],declination,2,1)
    uvwNS[:,i] = SRH.baseline2uvw(hourAngle[i],declination,192,191)
    uvwES[:,i] = SRH.baseline2uvw(hourAngle[i],declination,192,65)
    uvwWS[:,i] = SRH.baseline2uvw(hourAngle[i],declination,192,64)
    
PL.plot(NP.rad2deg(hourAngle),uvwEW[0,:])
PL.plot(NP.rad2deg(hourAngle),uvwNS[1,:])
PL.plot(NP.rad2deg(hourAngle),NP.sqrt(uvwES[0,:]**2 + uvwES[1,:]**2))
PL.plot(NP.rad2deg(hourAngle),NP.sqrt(uvwWS[0,:]**2 + uvwWS[1,:]**2))
PL.plot([NP.rad2deg(hourAngle[0]),NP.rad2deg(hourAngle[-1])],[SRH.dishDiameter,SRH.dishDiameter])


