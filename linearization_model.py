#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 15:32:03 2021

@author: mariagloba
"""

import numpy as NP
import pylab as PL
import random
from BadaryRAO import BadaryRAO
from scipy.optimize import least_squares
import base2uvw_36
from astropy.time import Time, TimeDelta


baselinesNumber = 4
antNumberEW = 5

N = 2048
arcsecPerPix = 2
radPerPix = NP.deg2rad(arcsecPerPix/3600.)
arcsecRadius = 1020
frequency = 3e9
degRadius = NP.deg2rad(arcsecRadius/3600)
radius = int(arcsecRadius/arcsecPerPix +0.5)
model = NP.zeros((N, N))
for i in range(N):
    for j in range(N):
        x=i - N/2
        y=j - N/2
        if (NP.sqrt(x**2 + y**2) < radius):
            model[i, j] = 1.
#model[1024,1024]= 10000
RAO = BadaryRAO('2021-05-14')
gains = NP.ones(5, dtype = 'complex')
for i in range(5):
    gains[i] = random.uniform(1., 2) * NP.exp(1j * random.uniform(-NP.pi/2., NP.pi/2.))
            
ewGains = gains

declination = RAO.declination
noon = RAO.culmination


fitsDate ='2021-05-14T00:00:00';
scan = 0

scanDate = Time(fitsDate, format='isot',scale='utc');
#scanTime = noon
scanTime = 25200.
#scanTime = 28800.
scanDate += TimeDelta(scanTime,format='sec')
hourAngle = NP.deg2rad((scanTime - noon)*15./3600.)
O = 1024//2
FOV = N * radPerPix
x,y = NP.meshgrid(NP.linspace(-.5,.5,N), NP.linspace(-.5,.5,N))
ewSolVis = NP.zeros(baselinesNumber, dtype = 'complex')
ewRedundantVis = NP.array(())
for i in range(baselinesNumber):
    baseline = i+1
    uvw = base2uvw_36.base2uvw(hourAngle,declination, 1, 1+baseline)
    cos_uv = NP.cos(2. * NP.pi * ((uvw[0]*frequency/3e8)*x + (uvw[1]*frequency/3e8)*y) * FOV)
    sin_uv = NP.sin(2. * NP.pi * ((uvw[0]*frequency/3e8)*x + (uvw[1]*frequency/3e8)*y) * FOV)
    real = NP.sum(cos_uv * model)
    imag = NP.sum(sin_uv * model)
    ewSolVis[i] = (real + imag * 1j)/1e6
    ewRedundantVis = NP.append(ewRedundantVis, ewSolVis[i] * NP.conj(ewGains[:5-baseline]) * ewGains[baseline:])
    
y1_0 = 1 + 0j
y2_0 = 1 + 0j
y3_0 = 1 + 0j
y4_0 = 1 + 0j

amp1_0 = 1
amp2_0 = 1
amp3_0 = 1
amp4_0 = 1
amp5_0 = 1

pha1_0 = 0
pha2_0 = 0
pha3_0 = 0
pha4_0 = 0
pha5_0 = 0

e0 = NP.zeros(5, dtype = 'complex')
e0[0] = NP.exp(amp1_0 + 1j*pha1_0)
e0[1] = NP.exp(amp2_0 + 1j*pha2_0)
e0[2] = NP.exp(amp3_0 + 1j*pha3_0)
e0[3] = NP.exp(amp4_0 + 1j*pha4_0)
e0[4] = NP.exp(amp5_0 + 1j*pha5_0)

e0_coefs = NP.hstack((NP.conj(e0[:4]) * e0[1:5], NP.conj(e0[:3]) * e0[2:5], NP.conj(e0[:2]) * e0[3:5], NP.conj(e0[1]) * e0[4]))
y0_coefs = NP.hstack((NP.full(4, y1_0), NP.full(3, y2_0), NP.full(2, y3_0), y4_0))

delta = ewRedundantVis.copy()
delta -= y0_coefs*e0_coefs

matr = [[1, 0, 0, 0, 1, 1, 0, 0, 0, -1j, 1j, 0, 0, 0],\
          [1, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1j, 1j, 0, 0],\
          [1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1j, 1j, 0],\
          [1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1j, 1j],\
          [0, 1, 0, 0, 1, 0, 1, 0, 0, -1j, 0, 1j, 0, 0],\
          [0, 1, 0, 0, 0, 1, 0, 1, 0, 0, -1j, 0, 1j, 0],\
          [0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, -1j, 0, 1j],\
          [0, 0, 1, 0, 1, 0, 0, 1, 0, -1j, 0, 0, 1j, 0],\
          [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, -1j, 0, 0, 1j],\
          [0, 0, 0, 1, 1, 0, 0, 0, 1, -1j, 0, 0, 0, 1j]]
    
matr2 = [[1, 1j, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1j, 1j, 0, 0, 0],\
          [1, 1j, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1j, 1j, 0, 0],\
          [1, 1j, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1j, 1j, 0],\
          [1, 1j, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1j, 1j],\
          [0, 0, 1, 1j, 0, 0, 0, 0, 1, 0, 1, 0, 0, -1j, 0, 1j, 0, 0],\
          [0, 0, 1, 1j, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, -1j, 0, 1j, 0],\
          [0, 0, 1, 1j, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, -1j, 0, 1j],\
          [0, 0, 0, 0, 1, 1j, 0, 0, 1, 0, 0, 1, 0, -1j, 0, 0, 1j, 0],\
          [0, 0, 0, 0, 1, 1j, 0, 0, 0, 1, 0, 0, 1, 0, -1j, 0, 0, 1j],\
          [0, 0, 0, 0, 0, 0, 1, 1j, 1, 0, 0, 0, 1, -1j, 0, 0, 0, 1j]]
    
matrix = NP.array(matr2)
matrix *= e0_coefs[:, NP.newaxis]
matrix[:, 8:] *= y0_coefs[:, NP.newaxis]

res, c,d,e = NP.linalg.lstsq(matrix, delta, rcond=None)


