#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 14 15:13:31 2021

@author: sergey_lesovoi
"""
import numpy as NP
import pylab as PL
import ZirinTb
from scipy.signal import find_peaks, peak_widths

Zi = ZirinTb.ZirinTb()

lcpImages = NP.zeros((16,512,512))
rcpImages = NP.zeros((16,512,512))

for ff in range(16):
    print(ff)
    phaseEdit.onFrequencyChannelChanged(ff)
    phaseEdit.onCenterDisk()
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

iSpectrum = NP.zeros((512,512,3))
iSpectrum[:,:,0] = (lcpImages+rcpImages)[0:5].mean(axis=0)
iSpectrum[:,:,1] = (lcpImages+rcpImages)[5:10].mean(axis=0)
iSpectrum[:,:,2] = (lcpImages+rcpImages)[10:16].mean(axis=0)
iSpectrum /= iSpectrum.max()

vSpectrum = NP.zeros((512,512,3))
vSpectrum[:,:,0] = (lcpImages-rcpImages)[0:5].mean(axis=0)
vSpectrum[:,:,1] = (lcpImages-rcpImages)[5:10].mean(axis=0)
vSpectrum[:,:,2] = (lcpImages-rcpImages)[10:16].mean(axis=0)
vSpectrum /= vSpectrum.max()
vSpectrum /= 2
vSpectrum += .5

fig = PL.figure(figsize=(10,10))
fig.suptitle(phaseEdit.srhFits.dateObs)
pl = fig.subplots(2,2)
PL.tight_layout()
pl[0,0].imshow(iSpectrum[:,:,0],origin='lower',cmap='Reds',vmin=0.0,vmax=0.6)
pl[0,0].axis('off')
pl[0,1].imshow(iSpectrum[:,:,1],origin='lower',cmap='Greens',vmin=0.0,vmax=0.6)
pl[0,1].axis('off')
pl[1,0].imshow(iSpectrum[:,:,2],origin='lower',cmap='Blues',vmin=0.0,vmax=0.6)
pl[1,0].axis('off')
pl[1,1].imshow(iSpectrum,origin='lower',vmin=0.0,vmax=0.3)
pl[1,1].axis('off')


fig = PL.figure(figsize=(10,10))
fig.suptitle(phaseEdit.srhFits.dateObs)
pl = fig.subplots(2,2)
PL.tight_layout()
pl[0,0].imshow(vSpectrum[:,:,0],origin='lower',cmap='Reds',vmin=0.0,vmax=0.6)
pl[0,0].axis('off')
pl[0,1].imshow(vSpectrum[:,:,1],origin='lower',cmap='Greens',vmin=0.0,vmax=0.6)
pl[0,1].axis('off')
pl[1,0].imshow(vSpectrum[:,:,2],origin='lower',cmap='Blues',vmin=0.0,vmax=0.6)
pl[1,0].axis('off')
pl[1,1].imshow(vSpectrum,origin='lower',vmin=0.0,vmax=1.0)
pl[1,1].axis('off')

iHist = []
iImages = 0.5*(lcpImages + rcpImages)
for f in range(16):
    iHist.append(NP.histogram(iImages[f],bins=1000,range=(-1e4,4e4)))

iHist = NP.array(iHist)

iWidths = []
for f in range(16):
    peaks, _ = find_peaks(iHist[f,0],prominence=1000)
    widths = peak_widths(iHist[f,0],peaks)
    iWidths.append(widths[0]*5e4/1e3)

outFile = open('srh_sens_' + phaseEdit.srhFits.dateObs + '.txt','w')
outFile.write('Histograms peaks of I images in frequency range\n')
outFile.write('Frequency [GHz] Sky [K] Sun [K]\n')
for f in range(16):
    if (len(iWidths[f]) > 1):
        outFile.write('%.1f %.1f %.1f\n'%(phaseEdit.srhFits.freqList[f]*1e-6,iWidths[f][0],iWidths[f][1]))
outFile.close()

