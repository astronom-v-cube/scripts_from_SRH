#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 07:20:37 2022

@author: sergey_lesovoi
"""

from astropy.io import fits
from ZirinTb import ZirinTb
from matplotlib.ticker import MultipleLocator
import pylab as PL
import numpy as NP

def freq0306format(f, pos):
    if (f > 0):
        return '%2.1f' % (f / 128 * 0.2 + 2.6)
    else:
        return ''

def freq0612format(f, pos):
    if (f > 0):
        return '%2.1f' % (f / 192 * 0.4 + 5.4)
    else:
        return ''

def freq1224format(f, pos):
    if (f > 0):
        return '%2.1f' % (f / 207 * 0.8 + 11.2)
    else:
        return ''
    
fitsAmps = fits.open('srh_amp_snr_20220425.fits')
antNumber0306 = int(fitsAmps[0].header['ANT_0306'])
antNumber0612 = int(fitsAmps[0].header['ANT_0612'])
antNumber1224 = int(fitsAmps[0].header['ANT_1224'])

freqs0306 = fitsAmps[1].data['freqs0306']*1e-6
freqs0612 = fitsAmps[1].data['freqs0612']*1e-6
freqs1224 = fitsAmps[1].data['freqs1224']*1e-6

time0306 = fitsAmps[2].data['time0306']
time0612 = fitsAmps[2].data['time0612']
time1224 = fitsAmps[2].data['time1224']


lcpAmp0306 = fitsAmps[3].data['lcpAmp0306'].reshape((freqs0306.shape[0],time0306.shape[1],antNumber0306))
rcpAmp0306 = fitsAmps[3].data['rcpAmp0306'].reshape((freqs0306.shape[0],time0306.shape[1],antNumber0306))
lcpAmp0612 = fitsAmps[4].data['lcpAmp0612'].reshape((freqs0612.shape[0],time0306.shape[1],antNumber0612))
rcpAmp0612 = fitsAmps[4].data['rcpAmp0612'].reshape((freqs0612.shape[0],time0306.shape[1],antNumber0612))
lcpAmp1224 = fitsAmps[5].data['lcpAmp1224'].reshape((freqs0612.shape[0],time0306.shape[1],antNumber1224))
rcpAmp1224 = fitsAmps[5].data['rcpAmp1224'].reshape((freqs0612.shape[0],time0306.shape[1],antNumber1224))

Zi = ZirinTb()
sunInd0 = 0
sunInd1 = 20
skyInd0 = 60
skyInd1 = 80

Tref = 290

#beamVsFreqFactor = NP.linspace(0.1,0.4,freqs0306.shape[0])*.5
beamVsFreqFactor = NP.linspace(0.0625,0.25,freqs0306.shape[0])*.7

for freq in range(freqs0306.shape[0]):
    for ant in range(antNumber0306):
        sunLevel = lcpAmp0306[freq,sunInd0:sunInd1,ant].mean()
        skyLevel = lcpAmp0306[freq,skyInd0:skyInd1,ant].mean()
        lcpAmp0306[freq,:,ant] *= Zi.getTbAtFrequency(freqs0306[freq]) / (sunLevel - skyLevel) * 1e3 * beamVsFreqFactor[freq]
        rcpAmp0306[freq,:,ant] *= Zi.getTbAtFrequency(freqs0306[freq]) / (sunLevel - skyLevel) * 1e3 * beamVsFreqFactor[freq]
        
for freq in range(freqs0612.shape[0]):
    for ant in range(antNumber0612):
        sunLevel = lcpAmp0612[freq,sunInd0:sunInd1,ant].mean()
        skyLevel = lcpAmp0612[freq,skyInd0:skyInd1,ant].mean()
        lcpAmp0612[freq,:,ant] *= Zi.getTbAtFrequency(freqs0612[freq]) / (sunLevel - skyLevel) * 1e3 * beamVsFreqFactor[freq]
        rcpAmp0612[freq,:,ant] *= Zi.getTbAtFrequency(freqs0612[freq]) / (sunLevel - skyLevel) * 1e3 * beamVsFreqFactor[freq]
        
for freq in range(freqs1224.shape[0]):
    for ant in range(antNumber1224):
        sunLevel = lcpAmp1224[freq,sunInd0:sunInd1,ant].mean()
        skyLevel = lcpAmp1224[freq,skyInd0:skyInd1,ant].mean()
        lcpAmp1224[freq,:,ant] *= Zi.getTbAtFrequency(freqs1224[freq]) / (sunLevel - skyLevel) * 1e3 * beamVsFreqFactor[freq]
        rcpAmp1224[freq,:,ant] *= Zi.getTbAtFrequency(freqs1224[freq]) / (sunLevel - skyLevel) * 1e3 * beamVsFreqFactor[freq]

flag0612 = [45, 104, 112, 189]
lcpAmp0612[:,:,flag0612] = 0

flag1224 = [1, 77, 99, 120, 134, 140, 194]
lcpAmp1224[:,:,flag1224] = 0
rcpAmp1224[:,:,flag1224] = 0

full0306 = NP.zeros(freqs0306.shape[0] * antNumber0306)
sefd0306 = NP.zeros(freqs0306.shape[0] * antNumber0306)
for freq in range(freqs0306.shape[0]):
    full0306[freq*antNumber0306:(freq+1)*antNumber0306] = (lcpAmp0306[freq,sunInd0:sunInd1,:].mean(axis=0))
    sefd0306[freq*antNumber0306:(freq+1)*antNumber0306] = (lcpAmp0306[freq,skyInd0:skyInd1,:].mean(axis=0))

full0612 = NP.zeros(freqs0612.shape[0] * antNumber0612)
sefd0612 = NP.zeros(freqs0612.shape[0] * antNumber0612)
for freq in range(freqs0612.shape[0]):
    full0612[freq*antNumber0612:(freq+1)*antNumber0612] = (lcpAmp0612[freq,sunInd0:sunInd1,:].mean(axis=0))
    sefd0612[freq*antNumber0612:(freq+1)*antNumber0612] = (lcpAmp0612[freq,skyInd0:skyInd1,:].mean(axis=0))
    
full1224 = NP.zeros(freqs1224.shape[0] * antNumber1224)
sefd1224 = NP.zeros(freqs1224.shape[0] * antNumber1224)
for freq in range(freqs1224.shape[0]):
    full1224[freq*antNumber1224:(freq+1)*antNumber1224] = (lcpAmp1224[freq,sunInd0:sunInd1,:].mean(axis=0))
    sefd1224[freq*antNumber1224:(freq+1)*antNumber1224] = (lcpAmp1224[freq,skyInd0:skyInd1,:].mean(axis=0))

lcpAmp0306F1 = lcpAmp0306[1].mean(axis=1)
lcpAmp0306F14 = lcpAmp0306[14].mean(axis=1)
lcpAmp0612F1 = lcpAmp0612[1].mean(axis=1)
lcpAmp0612F14 = lcpAmp0612[14].mean(axis=1)
lcpAmp1224F1 = lcpAmp1224[1].mean(axis=1)
lcpAmp1224F14 = lcpAmp1224[14].mean(axis=1)

fig = PL.figure()
fig.suptitle('SRH sun->sky pointing 20220425')
pl = fig.subplots(nrows=3,ncols=2)
pl[0,0].plot(lcpAmp0306F1,color='red',label='%.1f GHz NF %.1f dB'%(freqs0306[1],10*NP.log10(lcpAmp0306F1[skyInd0] / Tref + 1)))
pl[0,1].plot(lcpAmp0306F14,color='blue',label='%.1f GHz NF %.1f dB'%(freqs0306[14],10*NP.log10(lcpAmp0306F14[skyInd0] / Tref + 1)))
pl[0,0].set_ylim(0,1.3e4)
pl[0,1].set_ylim(0,1.3e4)
pl[0,0].text(10,50,'sun')
pl[0,0].text(75,30,'sky')
pl[0,1].text(10,50,'sun')
pl[0,1].text(75,30,'sky')
pl[0,0].legend()
pl[0,1].legend()
pl[0,0].grid()
pl[0,1].grid()

pl[1,0].plot(lcpAmp0612F1,color='red',label='%.1f GHz NF %.1f dB'%(freqs0612[1],10*NP.log10(lcpAmp0612F1[skyInd0] / Tref + 1)))
pl[1,1].plot(lcpAmp0612F14,color='blue',label='%.1f GHz NF %.1f dB'%(freqs0612[14],10*NP.log10(lcpAmp0612F14[skyInd0] / Tref + 1)))
pl[1,0].set_ylim(0,1.3e4)
pl[1,1].set_ylim(0,1.3e4)
pl[1,0].text(10,200,'sun')
pl[1,0].text(75,150,'sky')
pl[1,1].text(10,500,'sun')
pl[1,1].text(75,150,'sky')
pl[1,0].legend()
pl[1,1].legend()
pl[1,0].grid()
pl[1,1].grid()

pl[2,0].plot(lcpAmp1224F1,color='red',label='%.1f GHz NF %.1f dB'%(freqs1224[1],10*NP.log10(lcpAmp1224F1[skyInd0] / Tref + 1)))
pl[2,1].plot(lcpAmp1224F14,color='blue',label='%.1f GHz NF %.1f dB'%(freqs1224[14],10*NP.log10(lcpAmp1224F14[skyInd0] / Tref + 1)))
pl[2,0].set_ylim(0,1.3e4)
pl[2,1].set_ylim(0,1.3e4)
pl[2,0].text(10,2000,'sun')
pl[2,0].text(75,1800,'sky')
pl[2,1].text(10,8000,'sun')
pl[2,1].text(75,7000,'sky')
pl[2,0].legend()
pl[2,1].legend()
pl[2,0].grid()
pl[2,1].grid()
