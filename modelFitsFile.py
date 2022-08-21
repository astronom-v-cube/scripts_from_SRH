#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 09:29:16 2020

@author: mariagloba
"""

from astropy.io import fits
import numpy as NP
from scipy import ndimage
from astropy import coordinates
from astropy import constants
from BadaryRAO import BadaryRAO
import base2uvw as bl2uvw
import random
import datetime


srhFits = fits.open('/home/mariagloba/Work/fits/20200510/mf_20200510_030200.fit')

gainsL = NP.ones(48, dtype = 'complex')
for i in range(48):
    gainsL[i] =  random.uniform(1., 4.) *NP.exp(1j * random.uniform(-1., 1.))
gainsR = NP.ones(48, dtype = 'complex')
for i in range(48):
    gainsR[i] =  random.uniform(1., 4.) *NP.exp(1j * random.uniform(-1., 1.))
    
N = 256
arcsecPerPix = 10
radPerPix = NP.deg2rad(arcsecPerPix/3600.)
arcsecRadius = 1020
degRadius = NP.deg2rad(arcsecRadius/3600)
radius = int(arcsecRadius/arcsecPerPix +0.5)
modelL = NP.zeros((N, N))
modelR = NP.zeros((N, N))
for i in range(N):
    for j in range(N):
        x=i - N/2
        y=j - N/2
        if (NP.sqrt(x**2 + y**2) < radius):
            modelL[i, j] = 1.
            modelR[i, j] = 1.
modelL[128, 128] = 1000
modelR[128, 129] = 1000
FOV = N * radPerPix
sizeOfUv = 256

dateObs = datetime.datetime.fromisoformat(srhFits[0].header['DATE_OBS'])
timeObs = dateObs.time()
timeSec = int(datetime.timedelta(hours=timeObs.hour, minutes=timeObs.minute, seconds=timeObs.second).total_seconds())
RAO = BadaryRAO(dateObs.date().isoformat())
omegaEarth = coordinates.earth.OMEGA_EARTH.to_value()
hourAngle = omegaEarth * (timeSec - RAO.culmination)
declination = RAO.declination
freq = 4e9

visibilitiesL = NP.zeros(880, dtype = 'complex')
visibilitiesR = NP.zeros(880, dtype = 'complex')
uv = NP.zeros((sizeOfUv,sizeOfUv),dtype=complex);
x,y = NP.meshgrid(NP.linspace(-.5,.5,N), NP.linspace(-.5,.5,N))
#antennaA = []
#antennaB = []

for vis in range(512):
    i = vis // 32
    j = vis % 32 
    antA = j + 49
    antB = 192 - i
    uvw = bl2uvw.base2uvw(hourAngle,declination, antA, antB)
    cos_uv = NP.cos(2. * NP.pi * ((uvw[0]*freq/3e8)*x + (uvw[1]*freq/3e8)*y) * FOV)
    sin_uv = NP.sin(2. * NP.pi * ((uvw[0]*freq/3e8)*x + (uvw[1]*freq/3e8)*y) * FOV)
    real = NP.mean(cos_uv * modelL)
    imag = NP.mean(sin_uv * modelL)
    visibility = (real + imag * 1j)
    ant_based = gainsL[j] * NP.conj(gainsL[i+32])
#    baseline_based = NP.exp(1j * random.gauss(0, 0.05)) * random.gauss(1, 0.05)
#    noise = random.gauss(0, noiseLevel) * NP.exp(1j * random.gauss(0, phaseNoiseLevel))
    visibilitiesL[vis] = visibility * ant_based# + noise
    
    real = NP.mean(cos_uv * modelR)
    imag = NP.mean(sin_uv * modelR)
    visibility = (real + imag * 1j)
    ant_based = gainsR[i+32] * NP.conj(gainsR[j])
#    baseline_based = NP.exp(1j * random.gauss(0, 0.05)) * random.gauss(1, 0.05)
#    noise = random.gauss(0, noiseLevel) * NP.exp(1j * random.gauss(0, phaseNoiseLevel))
    visibilitiesR[vis] = visibility * ant_based# + noise
#    antennaA.append(antA)
#    antennaB.append(antB)
    
for vis in range(15):
    antA = 192 - vis - 1
    antB = 192 - vis
    uvw = bl2uvw.base2uvw(hourAngle,declination, antA, antB)
    cos_uv = NP.cos(2. * NP.pi * ((uvw[0]*freq/3e8)*x + (uvw[1]*freq/3e8)*y) * FOV)
    sin_uv = NP.sin(2. * NP.pi * ((uvw[0]*freq/3e8)*x + (uvw[1]*freq/3e8)*y) * FOV)
    real = NP.mean(cos_uv * modelL)
    imag = NP.mean(sin_uv * modelL)
    visibility = (real + imag * 1j)
    ant_based = gainsL[vis+32+1] * NP.conj(gainsL[vis+32])
#    baseline_based = NP.exp(1j * random.gauss(0, 0.05)) * random.gauss(1, 0.05)
#    noise = random.gauss(0, noiseLevel) * NP.exp(1j * random.gauss(0, phaseNoiseLevel))
    visibilitiesL[vis+512] = visibility * ant_based# + noise
    
    real = NP.mean(cos_uv * modelR)
    imag = NP.mean(sin_uv * modelR)
    visibility = (real + imag * 1j)
    ant_based = gainsR[vis+32] * NP.conj(gainsR[vis+32+1])
#    baseline_based = NP.exp(1j * random.gauss(0, 0.05)) * random.gauss(1, 0.05)
#    noise = random.gauss(0, noiseLevel) * NP.exp(1j * random.gauss(0, phaseNoiseLevel))
    visibilitiesR[vis+512] = visibility * ant_based# + noise
#    antennaA.append(antA)
#    antennaB.append(antB)
    
for vis in range(31):
    antA = vis + 49 + 1
    antB = vis + 49
    uvw = bl2uvw.base2uvw(hourAngle,declination, antA, antB)
    cos_uv = NP.cos(2. * NP.pi * ((uvw[0]*freq/3e8)*x + (uvw[1]*freq/3e8)*y) * FOV)
    sin_uv = NP.sin(2. * NP.pi * ((uvw[0]*freq/3e8)*x + (uvw[1]*freq/3e8)*y) * FOV)
    real = NP.mean(cos_uv * modelL)
    imag = NP.mean(sin_uv * modelL)
    visibility = (real + imag * 1j)
    ant_based = gainsL[vis+1] * NP.conj(gainsL[vis])
#    baseline_based = NP.exp(1j * random.gauss(0, 0.05)) * random.gauss(1, 0.05)
#    noise = random.gauss(0, noiseLevel) * NP.exp(1j * random.gauss(0, phaseNoiseLevel))
    visibilitiesL[vis+512+15] = visibility * ant_based# + noise
    
    real = NP.mean(cos_uv * modelR)
    imag = NP.mean(sin_uv * modelR)
    visibility = (real + imag * 1j)
    ant_based = gainsR[vis] * NP.conj(gainsR[vis+1])
#    baseline_based = NP.exp(1j * random.gauss(0, 0.05)) * random.gauss(1, 0.05)
#    noise = random.gauss(0, noiseLevel) * NP.exp(1j * random.gauss(0, phaseNoiseLevel))
    visibilitiesR[vis+512+15] = visibility * ant_based# + noise
#    antennaA.append(antA)
#    antennaB.append(antB)
    
for vis in range(14):
    antA = 192 - vis - 2
    antB = 192 - vis
    uvw = bl2uvw.base2uvw(hourAngle,declination, antA, antB)
    cos_uv = NP.cos(2. * NP.pi * ((uvw[0]*freq/3e8)*x + (uvw[1]*freq/3e8)*y) * FOV)
    sin_uv = NP.sin(2. * NP.pi * ((uvw[0]*freq/3e8)*x + (uvw[1]*freq/3e8)*y) * FOV)
    real = NP.mean(cos_uv * modelL)
    imag = NP.mean(sin_uv * modelL)
    visibility = (real + imag * 1j)
    ant_based = gainsL[vis+32+2] * NP.conj(gainsL[vis+32])
#    baseline_based = NP.exp(1j * random.gauss(0, 0.05)) * random.gauss(1, 0.05)
#    noise = random.gauss(0, noiseLevel) * NP.exp(1j * random.gauss(0, phaseNoiseLevel))
    visibilitiesL[vis+512+15+31] = visibility * ant_based# + noise
    
    real = NP.mean(cos_uv * modelR)
    imag = NP.mean(sin_uv * modelR)
    visibility = (real + imag * 1j)
    ant_based = gainsR[vis+32+2] * NP.conj(gainsR[vis+32])
#    baseline_based = NP.exp(1j * random.gauss(0, 0.05)) * random.gauss(1, 0.05)
#    noise = random.gauss(0, noiseLevel) * NP.exp(1j * random.gauss(0, phaseNoiseLevel))
    visibilitiesR[vis+512+15+31] = visibility * ant_based# + noise
#    antennaA.append(antA)
#    antennaB.append(antB)   
    
for vis in range(30):
    antA = vis + 49 + 2
    antB = vis + 49
    uvw = bl2uvw.base2uvw(hourAngle,declination, antA, antB)
    cos_uv = NP.cos(2. * NP.pi * ((uvw[0]*freq/3e8)*x + (uvw[1]*freq/3e8)*y) * FOV)
    sin_uv = NP.sin(2. * NP.pi * ((uvw[0]*freq/3e8)*x + (uvw[1]*freq/3e8)*y) * FOV)
    real = NP.mean(cos_uv * modelL)
    imag = NP.mean(sin_uv * modelL)
    visibility = (real + imag * 1j)
    ant_based = gainsL[vis+2] * NP.conj(gainsL[vis])
#    baseline_based = NP.exp(1j * random.gauss(0, 0.05)) * random.gauss(1, 0.05)
#    noise = random.gauss(0, noiseLevel) * NP.exp(1j * random.gauss(0, phaseNoiseLevel))
    visibilitiesL[vis+512+15+31+14] = visibility * ant_based# + noise
    
    real = NP.mean(cos_uv * modelR)
    imag = NP.mean(sin_uv * modelR)
    visibility = (real + imag * 1j)
    ant_based = gainsR[vis+2] * NP.conj(gainsR[vis])
#    baseline_based = NP.exp(1j * random.gauss(0, 0.05)) * random.gauss(1, 0.05)
#    noise = random.gauss(0, noiseLevel) * NP.exp(1j * random.gauss(0, phaseNoiseLevel))
    visibilitiesR[vis+512+15+31+14] = visibility * ant_based# + noise
#    antennaA.append(antA)
#    antennaB.append(antB)
    
#antennaA = NP.concatenate(NP.array(antennaA), NP.zeros(278))
#antennaB = NP.concatenate(NP.array(antennaB), NP.zeros(278))
#antenna = NP.array((192, 191, 190, 189, 188, 187, 186, 185, 184, 183, 182, 181, 180,
#        179, 178, 177,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,
#         59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,
#         72,  73,  74,  75,  76,  77,  78,  79,  80))

srhFits[1].data['VIS_LCP'][0, :880] = visibilitiesL
srhFits[1].data['VIS_RCP'][0, :880] = visibilitiesR
srhFits[1].data['AMP_LCP'][0, :48] = NP.ones(48)
srhFits[1].data['AMP_RCP'][0, :48] = NP.ones(48)

srhFits.writeto('model_test.fit')

