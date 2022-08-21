#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 09:07:54 2020

@author: svlesovoi
"""
from astropy.io import fits
import pylab as PL
import numpy as NP
from scipy.interpolate import interp1d

rainbowColors = PL.get_cmap('rainbow')

fd = fits.open('srh_cp_20200529.fits')
frequencies = fd[1].data['frequencies']
freqMean = frequencies.mean()
t0_ind = 2750
dt_ind = 300
flareI = fd[2].data['I'][:,t0_ind:t0_ind+dt_ind]
flareV = fd[2].data['V'][:,t0_ind:t0_ind+dt_ind]
flareTime = fd[2].data['time'][:,t0_ind:t0_ind+dt_ind]
for i in range(32):
    flareI[i] -= flareI[i,0]
    flareI[i] *= (frequencies[i] / freqMean)**2#(i/31 + 1)
    flareV[i] *= (frequencies[i] / freqMean)**2#(i/31 + 1)
#for i in range(32):
#    flareI[i] = flareI[i] / flareI[i,flareI.shape[1] - 140:flareI.shape[1] - 1].mean()

intFlareTime = NP.linspace(0,1,dt_ind*32)*(flareTime[-1,-1] - flareTime[0,0]) + flareTime[0,0]
intFlareI = NP.zeros((32,dt_ind*32))
intFlareV = NP.zeros((32,dt_ind*32))
for i in range(32):
    intFlareI[i] = NP.interp(intFlareTime,flareTime[i],flareI[i])
    intFlareV[i] = NP.interp(intFlareTime,flareTime[i],flareV[i])
    
PL.figure()
for i in range(32): 
#    PL.plot(intFlareTime,intFlareI[i]*(i/31 + 1)+0.0*i,color=rainbowColors((32-i)*8))
    fIntFlareI = interp1d(flareTime[i],flareI[i],kind='cubic')
    fIntFlareV = interp1d(flareTime[i],flareV[i],kind='cubic')
    intFlareTime = NP.linspace(flareTime[i,0],flareTime[i,-1],dt_ind*32)
    PL.plot(intFlareTime,fIntFlareI(intFlareTime),color=rainbowColors((32-i)*8))
    PL.plot(intFlareTime,fIntFlareV(intFlareTime),color=rainbowColors((32-i)*8))

#PL.imshow(intFlareI,aspect=100,origin='lower')