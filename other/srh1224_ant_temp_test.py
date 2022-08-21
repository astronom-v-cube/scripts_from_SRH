#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 03:53:37 2022

@author: sergey_lesovoi
"""
import numpy as NP
import pylab as PL
from ZirinTb import ZirinTb

def fitGauss(x, A, B, C):
    return  A * NP.exp(-(x - B)**2 / C)
    
zirinObj = ZirinTb()

deltaT = phaseEdit.srhFits.freqTime[0,1] - phaseEdit.srhFits.freqTime[0,0]
deltaH = 15*deltaT

rcpArr = phaseEdit.srhFits.ampRcp[:,:,:].copy()
lcpArr = phaseEdit.srhFits.ampLcp[:,:,:].copy()

sunInd = [0,50]
skyInd = [113,125]
forestInd = [145,165]

rcpArr[:,:,16:31] = 0
lcpArr[:,:,16:31] = 0
rcpArr[:,145:165,16:31] = 1
lcpArr[:,145:165,16:31] = 1

antennaNumber = rcpArr.shape[2]
freqNumber = rcpArr.shape[0]

freqListGHz = phaseEdit.srhFits.freqList/1e6
freqFlux = zirinObj.getSfuAtFrequency(freqListGHz)

for ant in range(antennaNumber):
    for freq in range(freqNumber):
        curMin = rcpArr[freq,skyInd[0]:skyInd[1],ant].mean()
        curMax = rcpArr[freq,sunInd[0]:sunInd[1],ant].mean()
        rcpArr[freq,:,ant] = (rcpArr[freq,:,ant] - curMin)/(curMax - curMin)
        forestMean = rcpArr[freq,forestInd[0]:forestInd[1],ant].mean()
        rcpArr[freq,:,ant] *= 300/forestMean
        curMin = lcpArr[freq,skyInd[0]:skyInd[1],ant].mean()
        curMax = lcpArr[freq,sunInd[0]:sunInd[1],ant].mean()
        lcpArr[freq,:,ant] = (lcpArr[freq,:,ant] - curMin)/(curMax - curMin)
        forestMean = lcpArr[freq,forestInd[0]:forestInd[1],ant].mean()
        lcpArr[freq,:,ant] *= 300/forestMean

rcpAntT = rcpArr[:,0:50,:].mean(axis=1)
lcpAntT = lcpArr[:,0:50,:].mean(axis=1)
rcpAntT[:,16:31] = 1000
lcpAntT[:,16:31] = 1000

expectedAntT1 = zirinObj.getTbAtFrequency(freqListGHz[1])*(0.5**2/1.9**2)
expectedAntT12 = zirinObj.getTbAtFrequency(freqListGHz[12])*(0.5**2/1.4**2)

PL.figure(figsize=(8,6))
PL.title('SRH 20220315, antenna temperature')
PL.gca().set_xticks(NP.arange(0,207,16) - 1)
PL.plot(rcpAntT[1],'.',label = '%.1f GHz, mean %.1f, expected %.1f'%(freqListGHz[1],rcpAntT[1].mean(),expectedAntT1*1e3))
PL.plot(rcpAntT[12],'.',label = '%.1f GHz, mean %.1f, expected %.1f'%(freqListGHz[12],rcpAntT[12].mean(),expectedAntT12*1e3))
PL.legend()
PL.grid()
PL.xlabel('antenna index')
PL.ylabel('K')
PL.ylim(0,5000)
PL.xlim(0,207)
