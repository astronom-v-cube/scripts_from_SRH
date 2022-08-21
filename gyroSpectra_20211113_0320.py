#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 11:49:23 2021

@author: sergey_lesovoi
"""

import numpy as NP
import pylab as PL
import ZirinTb

def freqGHzFormat(f, pos):
    return '%3.2f'%(f * 0.4 + 5.8)

Zi = ZirinTb.ZirinTb()

lcpImages = NP.zeros((16,512,512))
rcpImages = NP.zeros((16,512,512))

for ff in range(16):
    phaseEdit.onFrequencyChannelChanged(ff)
#    phaseEdit.onCenterDisk()
    lcpImages[ff] = phaseEdit.lcpData
    skyLevel = lcpImages[ff,470:500,200:300].mean()
    sunLevel = lcpImages[ff,390:420,220:270].mean()
    qSunTb = 0.5*Zi.getTbAtFrequency(phaseEdit.srhFits.freqList[ff]*1e-6)*1e3
    lcpImages[ff] = (lcpImages[ff] - skyLevel)/(sunLevel - skyLevel) * qSunTb

    rcpImages[ff] = phaseEdit.rcpData
    skyLevel = rcpImages[ff,470:500,200:300].mean()
    sunLevel = rcpImages[ff,390:420,220:270].mean()
    qSunTb = 0.5*Zi.getTbAtFrequency(phaseEdit.srhFits.freqList[ff]*1e-6)*1e3
    rcpImages[ff] = (rcpImages[ff] - skyLevel)/(sunLevel - skyLevel) * qSunTb

fig = PL.figure()
pl = fig.subplots(4,16)

lcpMaxT1 = []
rcpMaxT1 = []
lcpMaxT2 = []
rcpMaxT2 = []

srcNwLcpT = NP.zeros((16,5,5))
srcNwRcpT = NP.zeros((16,5,5))
# srcNW = [300,420]
# srcSE = [140,280]
srcNW = [300,410]
srcSE = [140,240]
dXY = 30

for ii in range(16):
    pl[0,ii].imshow(lcpImages[ii,srcNW[0]:srcNW[0]+dXY,srcNW[1]:srcNW[1]+dXY],origin='lower',vmin=0,vmax=6e4)
    pl[0,ii].axis('off')
    pl[1,ii].imshow(rcpImages[ii,srcNW[0]:srcNW[0]+dXY,srcNW[1]:srcNW[1]+dXY],origin='lower',vmin=0,vmax=6e4)
    pl[1,ii].axis('off')
    pl[2,ii].imshow(lcpImages[ii,srcSE[0]:srcSE[0]+dXY,srcSE[1]:srcSE[1]+dXY],origin='lower',vmin=0,vmax=6e4)
    pl[2,ii].axis('off')
    pl[3,ii].imshow(rcpImages[ii,srcSE[0]:srcSE[0]+dXY,srcSE[1]:srcSE[1]+dXY],origin='lower',vmin=0,vmax=6e4)
    pl[3,ii].axis('off')

    lcpMaxT1.append(lcpImages[ii,srcNW[0]:srcNW[0]+dXY,srcNW[1]:srcNW[1]+dXY].max())
    rcpMaxT1.append(rcpImages[ii,srcNW[0]:srcNW[0]+dXY,srcNW[1]:srcNW[1]+dXY].max())
    lcpMaxT2.append(lcpImages[ii,srcSE[0]:srcSE[0]+dXY,srcSE[1]:srcSE[1]+dXY].max())
    rcpMaxT2.append(rcpImages[ii,srcSE[0]:srcSE[0]+dXY,srcSE[1]:srcSE[1]+dXY].max())
    
    lcpCell = lcpImages[ii,srcNW[0]:srcNW[0]+dXY,srcNW[1]:srcNW[1]+dXY]
    rcpCell = rcpImages[ii,srcNW[0]:srcNW[0]+dXY,srcNW[1]:srcNW[1]+dXY]
    for k in range(5):
        for j in range(5):
            maxInd = NP.unravel_index(NP.argmax(lcpCell),lcpCell.shape)
            srcNwLcpT[ii,k,j] = lcpCell[maxInd[0] + k - 2, maxInd[1] + j - 2]
            maxInd = NP.unravel_index(NP.argmax(rcpCell),rcpCell.shape)
            srcNwRcpT[ii,k,j] = rcpCell[maxInd[0] + k - 2, maxInd[1] + j - 2]
            
    
fig = PL.figure(figsize=(12,8))
pl = fig.subplots(1)
pl.xaxis.set_major_formatter(PL.FuncFormatter(freqGHzFormat))
pl.xaxis.set_major_locator(PL.MultipleLocator(1.0))

fig.suptitle(phaseEdit.srhFits.dateObs)
pl.plot(rcpMaxT1,'--',color='red', label='RCP NW')
pl.plot(lcpMaxT1,'--',color='blue', label='LCP NW')
pl.plot(rcpMaxT2,color='red', label='RCP SE')
pl.plot(lcpMaxT2,color='blue', label='LCP SE')
pl.set_xlabel('GHz')
pl.set_ylabel('Tb [K]')
pl.legend()
pl.grid()

Tmin = 3e3
Tmax = 2e4
fig = PL.figure(figsize=(30,5))
pl = fig.subplots(2,8)
PL.tight_layout()
for f in NP.arange(5,13,1):       
    pl[0,f-5].imshow(lcpImages[f],origin='lower',vmin=Tmin,vmax=Tmax,cmap='gnuplot')
    pl[0,f-5].axis('off')
    pl[1,f-5].imshow(rcpImages[f],origin='lower',vmin=Tmin,vmax=Tmax,cmap='gnuplot')
    pl[1,f-5].axis('off')
    