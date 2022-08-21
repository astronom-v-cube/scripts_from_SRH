#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 10:00:56 2021

@author: sergeyvlesovoi
"""

from astropy.io import fits
import numpy as NP
import pylab as PL
from ZirinTb import ZirinTb

def quiteSunSfuAsFrequncy(freq, freq0):
    return 80 + ((freq - freq0)*2E-6)**2

times = [13600,14200]
sky_times = [1000,1500]
sun_times = [2500,3000]
fluxAndCp = fits.open('srh_cp_20210704.fits')
frequencies = fluxAndCp[1].data['frequencies']

flareCp = fluxAndCp[2].data['I'][:,times[0]:times[1]].copy()
flareFlux = fluxAndCp[2].data['flux_I'][:,times[0]:times[1]].copy()
#flareFlux = fluxAndCp[2].data['flux_I'].copy()

meanSky = []
meanSun = []

ZT = ZirinTb()

for f in range(frequencies.shape[0]):
    meanSky.append(fluxAndCp[2].data['flux_I'][f,sky_times].mean())
    meanSun.append(fluxAndCp[2].data['flux_I'][f,sun_times].mean())

for f in range(frequencies.shape[0]):
    flareFlux[f,:] -= meanSky[f]
    flareFlux[f,:] /= (meanSun[f] - meanSky[f])
#    flareFlux[f,:] *= quiteSunSfuAsFrequncy(frequencies[f], frequencies[0]) 
    flareFlux[f,:] *= ZT.getSfuAtFrequncy(frequencies[f]*1e-6)

meanSky = []
meanSun = []

for f in range(frequencies.shape[0]):
    meanSky.append(fluxAndCp[2].data['I'][f,sky_times].mean())
    meanSun.append(fluxAndCp[2].data['I'][f,sun_times].mean())

for f in range(frequencies.shape[0]):
    flareCp[f,:] -= meanSky[f]
    flareCp[f,:] /= (meanSun[f] - meanSky[f])
    flareCp[f,:] *= quiteSunSfuAsFrequncy(frequencies[f], frequencies[0]) 

for f in range(frequencies.shape[0]):
    flareFlux[f,:] -= flareFlux[f,0]
    flareCp[f,:] -= flareCp[f,0]

PL.figure()
PL.imshow(flareCp,aspect=50)
PL.figure()
PL.imshow(flareFlux,aspect=50)

flareTb = flareCp + flareFlux
flareSize = flareFlux/flareTb
PL.figure()
#PL.ylim(0,1000)
for f in range(frequencies.shape[0]):
#    PL.plot(flareTb[f])
#    PL.plot(flareSize[f]*1000)
    PL.plot(flareFlux[f])
    #PL.plot(flareCp[f])
