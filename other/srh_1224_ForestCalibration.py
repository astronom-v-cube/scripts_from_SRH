#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 08:16:33 2022

@author: sergey_lesovoi
"""

import numpy as NP
import pylab as PL
from MoonTb import MoonTb
from scipy.stats import linregress
import matplotlib

mt = MoonTb()
zi = ZirinTb()
freqs = phaseEdit.srhFits.freqList * 1e-6


# sun0 = 0
# sun1 = 10
# sun0 = 0
# sun1 = 10
sun0 = 80
sun1 = 90

# sky0 = 70
# sky1 = 80
# sky0 = 95
# sky1 = 105
sky0 = 100
sky1 = 110

#moon0 = 298
#moon1 = 300
# moon0 = 285
# moon1 = 305
moon0 = 397
moon1 = 399

weVis0 = 11775
weVis1 = 11822

sVis0 = 9452
sVis1 = 9474

#FWHM at f_min 2.2 degree, f_max 1.5 degree, sun size 0.5 degree
psfCoef = NP.linspace(19,9,16)/1e3

rcpAmp = phaseEdit.srhFits.ampRcp[:,:,:]
lcpAmp = phaseEdit.srhFits.ampLcp[:,:,:]
iAmp = 0.5*(rcpAmp + lcpAmp)

shortestSVis = NP.abs(phaseEdit.srhFits.visLcp[:,:,sVis0:sVis1]).mean(axis=2) + NP.abs(phaseEdit.srhFits.visRcp[:,:,sVis0:sVis1]).mean(axis=2)
shortestWeVis = NP.abs(phaseEdit.srhFits.visLcp[:,:,weVis0:weVis1]).mean(axis=2) + NP.abs(phaseEdit.srhFits.visRcp[:,:,weVis0:weVis1]).mean(axis=2)

antSunStairs = iAmp[:,sun0:sun1,:].mean(axis=(0,1))-iAmp[:,sky0:sky1,:].mean(axis=(0,1))
antSunStairMean = antSunStairs.mean()

invalidAntsInds = NP.where(antSunStairs < 0.3*antSunStairMean)

validAnts = NP.delete(antSunStairs,invalidAntsInds)

validRcpAmp = NP.delete(rcpAmp,invalidAntsInds,axis=2)
validLcpAmp = NP.delete(lcpAmp,invalidAntsInds,axis=2)
validIAmp = 0.5*(validRcpAmp + validLcpAmp)

PL.figure()
sunZirinFlux = []
sunFlux = []
moonFlux = []
for ff in range(16):
    sunSkyMoon = validIAmp[ff,:,:].mean(axis=1)
    sunSkyMoon -= sunSkyMoon.min()
    sunSkyMoon /= sunSkyMoon[moon0:moon1].mean()
    sunStair = sunSkyMoon[sun0:sun1].mean() - sunSkyMoon[sky0:sky1].mean()
    sunFlux.append(sunStair*mt.getSfuAtFrequency(freqs[ff]))
    sunZirinFlux.append(zi.getSfuAtFrequency(freqs[ff]))
    moonFlux.append(mt.getSfuAtFrequency(freqs[ff]))

freqLabels = NP.linspace(12,24,16)
PL.plot(freqLabels,sunFlux,'.',label='SRH sun')
PL.plot(freqLabels,sunZirinFlux,label='Zirirn Sun')
PL.plot(freqLabels,moonFlux,label='Daywitt Moon')
PL.grid()
PL.legend()
PL.xlabel('GHz')
PL.ylabel('s.f.u.')
PL.title('SRH 20220331 the flux density, calibrated by the Moon observation')

PL.figure()
sunZirinTb = []
sunTb = []
moonTb = []
for ff in range(16):
    sunSkyMoon = validIAmp[ff,:,:].mean(axis=1)
    sunSkyMoon -= sunSkyMoon.min()
    sunSkyMoon /= sunSkyMoon[moon0:moon1].mean()
    sunStair = sunSkyMoon[sun0:sun1].mean() - sunSkyMoon[sky0:sky1].mean()
    sunTb.append(sunStair*mt.getTbAtFrequency(freqs[ff])/psfCoef[ff])
    sunZirinTb.append(zi.getTbAtFrequency(freqs[ff])/psfCoef[ff])
    moonTb.append(mt.getTbAtFrequency(freqs[ff])/psfCoef[ff])

freqLabels = NP.linspace(12,24,16)
PL.plot(freqLabels,sunTb,'.',label='SRH sun')
PL.plot(freqLabels,sunZirinTb,label='Zirirn Sun')
PL.plot(freqLabels,moonTb,label='Daywitt Moon')
PL.grid()
PL.legend()
PL.xlabel('GHz')
PL.ylabel('K')
PL.title('SRH 20220331 the Ta, calibrated by the Moon observation')

c_list = matplotlib.colors.LinearSegmentedColormap.from_list(PL.cm.datad['gist_rainbow'], colors=['r','g','b'], N = freqs.shape[0]);
PL.figure(figsize=(8,9))
for ff in range(16):
    sunSkyMoon = validIAmp[ff,:,:].mean(axis=1)
    # slope, intercept, A, B, C = linregress(NP.linspace(0,49,50), sunSkyMoon[0:50])
    # sunSkyMoon -= NP.linspace(0,sunSkyMoon.shape[0]-1,sunSkyMoon.shape[0])*slope + intercept
    sunSkyMoon /= sunSkyMoon[moon0:moon1].mean() - sunSkyMoon.min()
    sunSkyMoon *= (mt.getTbAtFrequency(freqs[ff])/psfCoef[ff])
    PL.plot(sunSkyMoon,label='%.1f GHz'%freqs[ff],color=c_list(ff))

arrow_style = {
    "head_width": 5.,
    "head_length": 100.,
    "color":"green"
}
PL.arrow(x=65,y=1000,dx=0,dy=1000,**arrow_style)
PL.text(50,800,'sun->sky')
PL.arrow(x=300,y=1000,dx=0,dy=1000,**arrow_style)
PL.text(280,800,'sky->moon')
PL.grid()
PL.legend()
PL.xlabel('sample number')
PL.ylabel('K')
PL.title('SRH %s the Ta, calibrated by the Moon observation'%phaseEdit.srhFits.dateObs)
PL.xlim(-100,500)
PL.ylim(0,8000)

forestSkyT0 = 90
forestSkyT1 = 100
forestT0 = 130
forestT1 = 140
moonSkyT0 = 280
moonSkyT1 = 290
moonT0 = 298
moonT1 = 302

forestTb = []
for ff in range(freqs.shape[0]):
    forestStair = iAmp[ff,forestT0:forestT1,:].mean() - iAmp[ff,forestSkyT0:forestSkyT1,:].mean()
    moonStair = iAmp0331[ff,moonT0:moonT1,:].mean() - iAmp0331[ff,moonSkyT0:moonSkyT1,:].mean()
    forestTb.append(moonStair/moonStair.max()*mt.getTbAtFrequency(freqs[ff])*1e3)

PL.figure()
PL.plot(freqLabels, forestTb, '.', label='Forest Az -30, Alt 3')
PL.xlabel('GHz')
PL.ylabel('K')
PL.ylim(200,220)
PL.grid()
PL.title('Forest (Az -40, Alt 3) brightness temperature')
