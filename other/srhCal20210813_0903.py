#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 10:00:56 2021

@author: sergeyvlesovoi
"""

from astropy.io import fits
import numpy as NP
import pylab as PL

def quiteSunSfuAsFrequncy(freq, freq0):
    return 80 + ((freq - freq0)*2E-6)**2

#times = [13600,14200]
times = [0,20000]
sky_times = [24650,24670]
sun_times = [2500,3000]
fluxAndCp = fits.open('srh_cp_20210813.fits')
frequencies = fluxAndCp[1].data['frequencies']

flareCp = fluxAndCp[2].data['I']#[:,times[0]:times[1]].copy()
flareFlux = fluxAndCp[2].data['flux_I']#[:,times[0]:times[1]].copy()

meanSky = []
meanSun = []

for f in range(frequencies.shape[0]):
    meanSky.append(fluxAndCp[2].data['flux_I'][f,sky_times].mean())
    meanSun.append(fluxAndCp[2].data['flux_I'][f,sun_times].mean())

for f in range(frequencies.shape[0]):
    flareFlux[f,:] -= meanSky[f]
    flareFlux[f,:] /= (meanSun[f] - meanSky[f])
    flareFlux[f,:] *= quiteSunSfuAsFrequncy(frequencies[f], frequencies[0]) 

#meanSky = []
#meanSun = []
#
#for f in range(frequencies.shape[0]):
#    meanSky.append(fluxAndCp[2].data['I'][f,sky_times].mean())
#    meanSun.append(fluxAndCp[2].data['I'][f,sun_times].mean())

#for f in range(frequencies.shape[0]):
#    flareCp[f,:] -= meanSky[f]
#    flareCp[f,:] /= (meanSun[f] - meanSky[f])
#    flareCp[f,:] *= quiteSunSfuAsFrequncy(frequencies[f], frequencies[0]) 


PL.figure()
for f in range(frequencies.shape[0]):
    PL.plot(flareFlux[f])

PL.figure()
for f in range(frequencies.shape[0]):
    PL.plot(flareCp[f])
