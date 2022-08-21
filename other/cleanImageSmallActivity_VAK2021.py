#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 07:28:28 2021

@author: svlesovoi
"""

import numpy as NP
import pylab as PL
import srh36Utils
from astropy.io import fits

iFitses = srh36Utils.findFits('VAK2021_smallActivity/','*I*.fit')
iFitses.sort()

iImages = []
iT0 = 1e5
vT0 = 1e4

for fitsName in iFitses:
    fitsHandle = fits.open(fitsName)
    iImages.append(fitsHandle[0].data)
    fitsHandle.close()

iImages = NP.array(iImages)

nCols = 5
fig, pl = PL.subplots(nrows=1,ncols=nCols,figsize=(12,4))
fig.subplots_adjust(left=0.01,right=0.99,top=0.99,bottom=0.01,hspace=0.01,wspace=0.01)
for col in range(nCols):
    pl[col].axis('off')
    meanI = iImages[3*col:3*col+3].mean(axis=0)
    pl[col].imshow(meanI,origin='lower',vmin=2e4,vmax=1e5,cmap='hot')
    pl[col].contour(meanI,origin='lower',levels=[1e5,1.5e5,2.0e5,2.5e5,3.e5],colors='black',linewidths=0.5)

