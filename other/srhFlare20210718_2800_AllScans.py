#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 00:39:04 2021

@author: svlesovoi
"""

import numpy as NP
import pylab as PL
import srh36Utils
from astropy.io import fits
from scipy.stats import linregress
from skimage.transform import warp, AffineTransform
import scipy
from ZirinTb import ZirinTb

frequency = 2800
frames = 20
#iFitses = srh36Utils.findFits('cleanMaps/20210703/%d'%frequency,'*I*.fit')
#vFitses = srh36Utils.findFits('cleanMaps/20210703/%d'%frequency,'*V*.fit')
iFitses = srh36Utils.findFits('cleanMaps/20210703/snapshot_0718/%d'%frequency,'*I*.fit')
vFitses = srh36Utils.findFits('cleanMaps/20210703/snapshot_0718/%d'%frequency,'*V*.fit')

iFitses.sort()
vFitses.sort()
iFitses = iFitses
vFitses = vFitses

iImages = []
vImages = []
time_0718_2800 = []
iT0 = 1e5
vT0 = 1e4

for fitsName in iFitses:
    fitsHandle = fits.open(fitsName)
    iImages.append(fitsHandle[0].data)
    tObs = fitsHandle[0].header['T-OBS'].split('T')[1].split('.')
    tHMS = tObs[0].split(':')
    hh = int(tHMS[0])
    mm = int(tHMS[1])
    ss = int(tHMS[2])
    time_0718_2800.append(hh*3600 + mm*60 + ss + int(tObs[1])*1e-6)
    fitsHandle.close()
for fitsName in vFitses:
    fitsHandle = fits.open(fitsName)
    vImages.append(fitsHandle[0].data)
    fitsHandle.close()

iImages = NP.array(iImages)
vImages = NP.array(vImages)
rImages = iImages + vImages
lImages = iImages - vImages

#iHist = NP.histogram(iImages[30],10000)
#PL.figure()
#PL.xlim(-1e5,1e5)
#PL.plot(iHist[1][1:],iHist[0])
#peaks, _ = scipy.signal.find_peaks(iHist[0],prominence = 200)
#skySunInd = peaks
#PL.plot(iHist[1][peaks+1],iHist[0][peaks],'.',color='green')
#    
#sunStair = 1
#skyLevel = 0
#if(skySunInd.shape[0] > 1):
#    sunStair = iHist[1][skySunInd[1]] - iHist[1][skySunInd[0]]
#    skyLevel = iHist[1][skySunInd[0]]
#for i in range(len(iFitses)): 
#    iImages[i] -= skyLevel
#    iImages[i] /= sunStair
#    iImages[i] *= qSunTb
#    vImages[i] /= sunStair
#    vImages[i] *= qSunTb
#

iHist = []
for i in range(len(iFitses)):
    iI = iImages[i].copy()
    iHist.append(NP.histogram(iImages[i],5000))

skySunInd = []
PL.figure()
for i in range(len(iFitses)): 
    PL.xlim(-1e5,1e5)
    PL.plot(iHist[i][1][1:],iHist[i][0] + i*10000)
    peaks, _ = scipy.signal.find_peaks(iHist[i][0],prominence = 400)
    ind1 = peaks[0]
    skySunInd.append(peaks)
    PL.plot(iHist[i][1][peaks+1],iHist[i][0][peaks] + i*10000,'.',color='green')
    
sunStair = []
skyLevel = []
for i in range(len(iFitses)): 
    if(skySunInd[i].shape[0] > 1):
        sunStair.append(iHist[i][1][skySunInd[i][1]] - iHist[i][1][skySunInd[i][0]])
        skyLevel.append(iHist[i][1][skySunInd[i][0]])

sunStair = NP.array(sunStair)
skyLevel = NP.array(skyLevel)

ZirinQSunTb = ZirinTb()
qSunTb = ZirinQSunTb.getTbAtFrequncy(frequency/1e3)*1e3

#for i in range(len(iFitses)): 
#    iImages[i] -= skyLevel[i]
#    iImages[i] /= sunStair[i]
#    iImages[i] *= qSunTb
#    vImages[i] /= sunStair[i]
#    vImages[i] *= qSunTb

#iSpectrum = NP.zeros((20,512,512,3))
#for i in range(20):
#    iSpectrum[i,:,:,0] = NP.abs(iImages2800[i]/iImages2800[i].max())**.5
#    iSpectrum[i,:,:,1] = NP.abs(iImages3100[i]/iImages3100[i].max())**.5
#    iSpectrum[i,:,:,2] = NP.abs(iImages3900[i]/iImages3900[i].max())**.5
    
iImages2800 = iImages

#for i in range(20,frames):
#    iImages3400[i] = NP.roll(iImages[i],(-358+382,-325+337),axis=(0,1))
    
nCols=10
nRows=frames//nCols
fig, pl = PL.subplots(nrows=nRows,ncols=nCols,figsize=(12,4))
fig.subplots_adjust(left=0.01,right=0.99,top=0.95,bottom=0.01,hspace=0.01,wspace=0.01)
fig.suptitle('%d'%frequency)
for row in range(nRows):
    for col in range(nCols):
        pl[row,col].axis('off')
        pl[row,col].imshow(iImages2800[row*nCols + col],origin='lower',cmap='hot',vmin=0,vmax=5e5)

flux_0718_2800 = []
Tbmax_0718_2800 = []
for i in range(frames):
    Tbmax_0718_2800.append(iImages2800[i,325,420])
    flux_0718_2800.append(2*iImages2800[i,303:350,402:438].sum()*scipy.constants.k/(scipy.constants.c/3.4e9)**2 * NP.deg2rad(4.911/3600)**2 / 1e-22)
