#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 01:37:52 2021

@author: svlesovoi
"""

import numpy as NP
import pylab as PL
import srh36Utils
from astropy.io import fits
from scipy.stats import linregress

iFitses2800 = srh36Utils.findFits('SRH36_temp_20210815_1/','*I*2800*.fit')
vFitses2800 = srh36Utils.findFits('SRH36_temp_20210815_1/','*V*2800*.fit')
iFitses3100 = srh36Utils.findFits('SRH36_temp_20210815_1/','*I*3100*.fit')
vFitses3100 = srh36Utils.findFits('SRH36_temp_20210815_1/','*V*3100*.fit')
iFitses3400 = srh36Utils.findFits('SRH36_temp_20210815_1/','*I*3400*.fit')
vFitses3400 = srh36Utils.findFits('SRH36_temp_20210815_1/','*V*3400*.fit')

iImages2800 = []
vImages2800 = []
iImages3100 = []
vImages3100 = []
iImages3400 = []
vImages3400 = []
iT0 = 1e5
vT0 = 1e4
for fitsName in iFitses2800:
    fitsHandle = fits.open(fitsName)
    iImages2800.append(fitsHandle[0].data)
for fitsName in vFitses2800:
    fitsHandle = fits.open(fitsName)
    vImages2800.append(fitsHandle[0].data)
for fitsName in iFitses3100:
    fitsHandle = fits.open(fitsName)
    iImages3100.append(fitsHandle[0].data)
for fitsName in vFitses3100:
    fitsHandle = fits.open(fitsName)
    vImages3100.append(fitsHandle[0].data)
for fitsName in iFitses3400:
    fitsHandle = fits.open(fitsName)
    iImages3400.append(fitsHandle[0].data)
for fitsName in vFitses3400:
    fitsHandle = fits.open(fitsName)
    vImages3400.append(fitsHandle[0].data)
    
iMeanImage2800 = NP.array(iImages2800).mean(axis=0)
vMeanImage2800 = NP.array(vImages2800).mean(axis=0)
iMeanImage3100 = NP.array(iImages3100).mean(axis=0)
vMeanImage3100 = NP.array(vImages3100).mean(axis=0)
iMeanImage3400 = NP.array(iImages3400).mean(axis=0)
vMeanImage3400 = NP.array(vImages3400).mean(axis=0)

#PL.figure()
#levels = NP.linspace(iT0,NP.max(iMeanImage2800),5)
#PL.imshow(iMeanImage2800,origin='lower',cmap='hot',vmin=0,vmax=iT0)
#PL.contour(iMeanImage2800,levels=levels,cmap='hot',linewidths=0.5)
#
#PL.figure()
#PL.imshow(iMeanImage3100,origin='lower',cmap='hot',vmin=0,vmax=iT0)
#PL.contour(iMeanImage3100,levels=levels,cmap='hot',linewidths=0.5)
#
#PL.figure()
#PL.imshow(iMeanImage3400,origin='lower',cmap='hot',vmin=0,vmax=iT0)
#PL.contour(iMeanImage3400,levels=levels,cmap='hot',linewidths=0.5)

vT0 = 3e3
maxV = 1.*NP.max(NP.abs(vMeanImage3400))
nLevels = NP.linspace(-maxV,-vT0,5)
pLevels = NP.linspace(vT0,maxV,5)

PL.figure()
PL.imshow(iMeanImage2800,origin='lower',cmap='hot',vmin=0,vmax=iT0)
PL.contour(vMeanImage2800,levels=nLevels,colors='blue',linewidths=0.5)
PL.contour(vMeanImage2800,levels=pLevels,colors='red',linewidths=0.5)

PL.figure()
PL.imshow(iMeanImage3100,origin='lower',cmap='hot',vmin=0,vmax=iT0)
PL.contour(vMeanImage3100,levels=nLevels,colors='blue',linewidths=0.5)
PL.contour(vMeanImage3100,levels=pLevels,colors='red',linewidths=0.5)

PL.figure()
PL.imshow(iMeanImage3400,origin='lower',cmap='hot',vmin=0,vmax=iT0)
PL.contour(vMeanImage3400,levels=nLevels,colors='blue',linewidths=0.5)
PL.contour(vMeanImage3400,levels=pLevels,colors='red',linewidths=0.5)

