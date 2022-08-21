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

rcpWeVis = NP.abs(phaseEdit.srhFits.visRcp[:,:,11730:11867].copy())
lcpWeVis = NP.abs(phaseEdit.srhFits.visLcp[:,:,11730:11867].copy())
rcpSVis = NP.abs(phaseEdit.srhFits.visRcp[:,:,9452:9518].copy())
lcpSVis = NP.abs(phaseEdit.srhFits.visLcp[:,:,9452:9518].copy())

sunInd = [0,35]
skyInd = [170,205]

rcpArr[:,:,16:31] = 0
lcpArr[:,:,16:31] = 0
rcpArr[:,0,16:31] = 1
lcpArr[:,0,16:31] = 1

antennaNumber = rcpArr.shape[2]
freqNumber = rcpArr.shape[0]

freqListGHz = phaseEdit.srhFits.freqList/1e6
freqFlux = zirinObj.getSfuAtFrequency(freqListGHz)

for ant in range(antennaNumber):
    for freq in range(freqNumber):
        curMin = rcpArr[freq,skyInd[0]:skyInd[1],ant].mean()
        curMax = rcpArr[freq,sunInd[0]:sunInd[1],ant].mean()
        rcpArr[freq,:,ant] = (rcpArr[freq,:,ant] - curMin)/(curMax - curMin)
        curMin = lcpArr[freq,skyInd[0]:skyInd[1],ant].mean()
        curMax = lcpArr[freq,sunInd[0]:sunInd[1],ant].mean()
        lcpArr[freq,:,ant] = (lcpArr[freq,:,ant] - curMin)/(curMax - curMin)

antSnr = NP.zeros((antennaNumber,freqNumber))
for ant in range(antennaNumber):
    for freq in range(freqNumber):
        antSnr[ant,freq] = NP.std(rcpArr[freq,sunInd[0]:sunInd[1],ant])
        if (antSnr[ant,freq] > 0.2):
            antSnr[ant,freq] = 0.01

PL.figure(figsize=(8,6))
PL.title('SRH 20220313, antenna sensitivity')
PL.gca().set_xticks(NP.arange(0,207,16) - 1)
PL.plot(freqFlux[1]*antSnr[:,1],'.',label='%.1f GHz, %.1f s.f.u.'%(freqListGHz[1],freqFlux[1]*antSnr[32:,1].mean(axis=0)))
PL.plot(freqFlux[14]*antSnr[:,14],'.',label='%.1f GHz, %.1f s.f.u.'%(freqListGHz[14],freqFlux[14]*antSnr[32:,14].mean(axis=0)))
PL.legend()
PL.grid()
PL.xlabel('antenna index')
PL.ylabel('s.f.u.')
PL.ylim(0,50)
PL.xlim(0,207)

xlabel = ["{0:.1f}".format(phaseEdit.srhFits.freqList[f] * 1e-6) for f in range(freqNumber)]
fig = PL.figure(figsize=(8,6))
pl = fig.subplots()
fig.suptitle('SRH 20220313, antenna sensitivity mean versus frequency')
pl.set_xticks(NP.arange(0,16,1))
pl.set_xticklabels(xlabel)
pl.plot(antSnr[32:].mean(axis=0)*freqFlux,'-',label='s.f.u.',color='red')
pl.legend(loc='upper left')
pl.grid()
pl.set_xlabel('GHz')
pl.set_ylabel('s.f.u.')
pl.set_xlim(0,15)
pl.set_ylim(0,50)
pl_t = pl.twinx()
pl_t.plot(antSnr[32:].mean(axis=0)*100,color='blue',label='%')
pl_t.set_ylim(0,5)
pl_t.set_ylabel('%')
pl_t.legend(loc='upper right')


visSunInd = [105,145]
visSkyInd = [160,200]

ewVisNumber = 46
rcpVisEwSnr = NP.zeros((ewVisNumber,freqNumber))
lcpVisEwSnr = NP.zeros((ewVisNumber,freqNumber))
for vis in range(ewVisNumber):
    for freq in range(freqNumber):
        visSun = NP.mean(rcpWeVis[freq,visSunInd[0]:visSunInd[1],vis + 47])
        visSky = NP.std(rcpWeVis[freq,visSkyInd[0]:visSkyInd[1],vis + 47])
        rcpVisEwSnr[vis,freq] = visSun / visSky
        visSun = NP.mean(lcpWeVis[freq,visSunInd[0]:visSunInd[1],vis + 47])
        visSky = NP.std(lcpWeVis[freq,visSkyInd[0]:visSkyInd[1],vis + 47])
        lcpVisEwSnr[vis,freq] = visSun / visSky

sVisNumber = 22
rcpVisSSnr = NP.zeros((sVisNumber,freqNumber))
lcpVisSSnr = NP.zeros((sVisNumber,freqNumber))
for vis in range(sVisNumber):
    for freq in range(freqNumber):
        visSun = NP.mean(rcpSVis[freq,visSunInd[0]:visSunInd[1],vis])
        visSky = NP.std(rcpSVis[freq,visSkyInd[0]:visSkyInd[1],vis])
        rcpVisSSnr[vis,freq] = visSun / visSky
        visSun = NP.mean(lcpSVis[freq,visSunInd[0]:visSunInd[1],vis])
        visSky = NP.std(lcpSVis[freq,visSkyInd[0]:visSkyInd[1],vis])
        lcpVisSSnr[vis,freq] = visSun / visSky

fig = PL.figure(figsize=(8,6))
pl = fig.subplots()
fig.suptitle('SRH 20220313, shortest baselines SNR versus frequency')
pl.set_xticks(NP.arange(0,16,1))
pl.set_xticklabels(xlabel)
pl.plot(0.5*(rcpVisSSnr.mean(axis=0) + lcpVisSSnr.mean(axis=0)),label='South 1b',color='blue')
pl.plot(0.5*(rcpVisEwSnr.mean(axis=0) + lcpVisEwSnr.mean(axis=0)),label='West-East 1b',color='red')
pl.legend(loc='upper right')
pl.grid()
pl.set_xlabel('GHz')
pl.set_ylabel('(mean of the Sun) / (std of the sky)')
