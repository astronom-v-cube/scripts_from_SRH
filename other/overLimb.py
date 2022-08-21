#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 15:22:47 2021

@author: sergey_lesovoi
"""

import numpy as NP
import pylab as PL

N = 720
H = 50
XY0 = 256
R0 = 190
overLimb = NP.zeros((16,N,H))

alpha = 2*NP.pi*NP.linspace(0,N,N)/N
for ff in range(16):
    print(ff)
    phaseEdit.onFrequencyChannelChanged(ff)
#    phaseEdit.onCenterDisk()
    for radius in range(H):
        indX = XY0 + NP.array(NP.ceil((radius + R0)*NP.cos(alpha)),dtype='int')
        indY = XY0 + NP.array(NP.ceil((radius + R0)*NP.sin(alpha)),dtype='int')
        overLimb[ff,:,radius] = phaseEdit.lcpData[indY,indX]

PL.figure()
PL.imshow(NP.concatenate(overLimb,axis=1).T,origin='lower')
