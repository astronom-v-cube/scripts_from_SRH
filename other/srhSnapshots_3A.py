#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 03:28:47 2022

@author: sergey_lesovoi
"""
from astropy.io import fits
from ZirinTb import ZirinTb
from matplotlib.ticker import MultipleLocator
import pylab as PL
import numpy as NP

Zi = ZirinTb()

fitsSnapshot0306 = fits.open('images_maria/36_fits/36_0000_I.fits')
fitsSnapshot0612 = fits.open('images_maria/612_fits/612_000_I.fits')
fitsSnapshot1224 = fits.open('images_maria/1224_fits/1224_000_I.fits')

skyY0 = 10
skyY1 = 45
skyX0 = 110
skyX1 = 350

sunY0 = 90
sunY1 = 130
sunX0 = 190
sunX1 = 320

snapshot0306 = fitsSnapshot0306[0].data
snapshot0612 = fitsSnapshot0612[0].data
snapshot1224 = fitsSnapshot1224[0].data

sky0306 = snapshot0306[skyY0:skyY1,skyX0:skyX1].mean()
sky0612 = snapshot0612[skyY0:skyY1,skyX0:skyX1].mean()
sky1224 = snapshot1224[skyY0:skyY1,skyX0:skyX1].mean()

sun0306 = snapshot0306[sunY0:sunY1,sunX0:sunX1].mean()
sun0612 = snapshot0612[sunY0:sunY1,sunX0:sunX1].mean()
sun1224 = snapshot1224[sunY0:sunY1,sunX0:sunX1].mean()

tSnapshot0306 = (snapshot0306 - sky0306)/(sun0306 - sky0306) * Zi.getTbAtFrequency(3.0) * 1e3
tSnapshot0612 = (snapshot0612 - sky0612)/(sun0612 - sky0612) * Zi.getTbAtFrequency(6.0) * 1e3
tSnapshot1224 = (snapshot1224 - sky1224)/(sun1224 - sky1224) * Zi.getTbAtFrequency(12.0) * 1e3

tHist0306 = NP.histogram(tSnapshot0306, bins=10000)
tHist0612 = NP.histogram(tSnapshot0612, bins=10000)
tHist1224 = NP.histogram(tSnapshot1224, bins=10000)

logtSnapshot0306 = NP.nan_to_num(NP.log10(tSnapshot0306),nan=0.0,copy=True)
logtSnapshot0612 = NP.nan_to_num(NP.log10(tSnapshot0612),nan=0.0,copy=True)
logtSnapshot1224 = NP.nan_to_num(NP.log10(tSnapshot1224),nan=0.0,copy=True)

fig = PL.figure(figsize=(12,8))
fig.suptitle('SRH snapshots 20220425')
pl = fig.subplots(nrows=2,ncols=3)
pl[0,0].imshow(logtSnapshot0306,origin='lower',cmap='hot',vmin=3,vmax=6)
pl[0,1].imshow(logtSnapshot0612,origin='lower',cmap='hot',vmin=3,vmax=6)
pl[0,2].imshow(logtSnapshot1224,origin='lower',cmap='hot',vmin=3,vmax=6)
pl[0,0].axis('off')
pl[0,1].axis('off')
pl[0,2].axis('off')

pl[1,0].plot(tHist0306[1][1:]*1e-3,tHist0306[0],label='03-06 GHz')
pl[1,1].plot(tHist0612[1][1:]*1e-3,tHist0612[0],label='06-12 GHz')
pl[1,2].plot(tHist1224[1][1:]*1e-3,tHist1224[0],label='12-24 GHz')
pl[1,0].set_xlim(-40,80)
pl[1,1].set_xlim(-40,80)
pl[1,2].set_xlim(-40,80)
for cc in range(3):
    pl[1,cc].set_ylim(0,3000)
    pl[1,cc].grid()
    pl[1,cc].legend()
    pl[1,cc].set_xlabel('1E3 K')
    
