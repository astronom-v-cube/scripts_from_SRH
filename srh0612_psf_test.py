#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 03:53:37 2022

@author: sergey_lesovoi
"""
import numpy as NP
import pylab as PL
import scipy.optimize as OPT

def fitGauss(x, A, B, C):
    return  A * NP.exp(-(x - B)**2 / C)
    
deltaT = phaseEdit.srhFits.freqTime[0,1] - phaseEdit.srhFits.freqTime[0,0]
deltaH = 15*deltaT

# rcpArr = phaseEdit.srhFits.ampRcp[:,:,:].copy()
# lcpArr = phaseEdit.srhFits.ampLcp[:,:,:].copy()
rcpArr = phaseEdit.srhFits.ampRcp[:,50:270,:].copy()
lcpArr = phaseEdit.srhFits.ampLcp[:,50:270,:].copy()

# rcpArr[:,:,16:31] = 0
# lcpArr[:,:,16:31] = 0
# rcpArr[:,0,16:31] = 1
# lcpArr[:,0,16:31] = 1

antennaNumber = rcpArr.shape[2]
freqNumber = rcpArr.shape[0]

for ant in range(antennaNumber):
    for freq in range(freqNumber):
        curMin = rcpArr[freq,:,ant].min()
        curMax = rcpArr[freq,:,ant].max()
        rcpArr[freq,:,ant] = (rcpArr[freq,:,ant] - curMin)/(curMax - curMin)
        curMin = lcpArr[freq,:,ant].min()
        curMax = lcpArr[freq,:,ant].max()
        lcpArr[freq,:,ant] = (lcpArr[freq,:,ant] - curMin)/(curMax - curMin)

antMaxInd = NP.zeros((antennaNumber,freqNumber))
psfCenter = NP.zeros((antennaNumber,freqNumber))
psfWidth = NP.zeros((antennaNumber,freqNumber))
for ant in range(antennaNumber):
    for freq in range(freqNumber):
        antMaxInd[ant,freq] = NP.argmax(rcpArr[freq,:,ant])
        try:
            fitPar = OPT.curve_fit(fitGauss,NP.linspace(0,219,220),rcpArr[freq,:,ant],p0=[1,100,50])
            psfCenter[ant,freq] = fitPar[0][1]
            if (psfCenter[ant,freq] > 5000):
                psfCenter[ant,freq] = 100
            psfWidth[ant,freq] = 2.355*NP.sqrt(fitPar[0][2])
            if (psfWidth[ant,freq] > 5000):
                psfWidth[ant,freq] = 50
        except:
            pass

PL.figure()
PL.title('SRH 6-12 20220314, offsets of the PSF Gaussian fit')
PL.gca().set_xticks(NP.arange(0,192,16) - 1)
PL.plot(deltaH/3600*psfCenter[:,1], '.',label='%d MHz'%(phaseEdit.srhFits.freqList[1]/1e3))
PL.plot(deltaH/3600*psfCenter[:,14],'.',label='%d MHz'%(phaseEdit.srhFits.freqList[14]/1e3))
PL.legend()
PL.grid()
PL.xlabel('antenna index')
PL.ylabel('degree')
PL.ylim(0,3)
PL.xlim(0,192)

PL.figure()
PL.title('SRH 6-12 20220314, FWHM of the PSF Gaussian fit')
PL.gca().set_xticks(NP.arange(0,192,16) - 1)
PL.plot(deltaH/3600*psfWidth[:,1], '.',label='%d MHz'%(phaseEdit.srhFits.freqList[1]/1e3))
PL.plot(deltaH/3600*psfWidth[:,14],'.',label='%d MHz'%(phaseEdit.srhFits.freqList[14]/1e3))
PL.legend()
PL.grid()
PL.xlabel('antenna index')
PL.ylabel('degree')
PL.ylim(0,3)
PL.xlim(0,192)

xlabel = ["{0:.1f}".format(phaseEdit.srhFits.freqList[f] * 1e-6) for f in range(freqNumber)]
PL.figure(figsize=(8,6))
PL.title('SRH 6-12 20220314, offsets and FWHM of the PSF Gaussian fit')
PL.gca().set_xticks(NP.arange(0,16,1))
PL.gca().set_xticklabels(xlabel)
PL.plot(deltaH/3600*psfCenter[:,:].mean(axis=0),'-',label='offset')
PL.plot(deltaH/3600*psfWidth[:,:].mean(axis=0),'-',label='FWHM')
PL.legend()
PL.grid()
PL.xlabel('GHz')
PL.ylabel('degree')
PL.ylim(0,3)
PL.xlim(0,15)