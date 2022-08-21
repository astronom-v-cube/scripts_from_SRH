#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 07:52:23 2018

@author: sergey
"""

from astropy.io import fits
import pylab as PL
import numpy as NP
import scipy.signal
import datetime as DT
import base2uvw as b2u
from BadaryRAO import BadaryRAO

def shift2D(arr):
    return NP.roll(arr, (arr.shape[0]//2, arr.shape[0]//2), axis=(0,1))
                   
def fftConvolution(arr1, arr2):
    size = arr1.shape[0]
    return NP.roll((NP.fft.ifft2(NP.fft.fft2(arr1) * NP.conjugate(NP.fft.fft2(arr2)))),(size//2,size//2),axis=(0,1))
    
def fftBeam(uvArr):
    return shift2D(NP.fft.fft2(shift2D(uvArr)))

N = 64
M = 16*N
L = 100

uvPlain0 = NP.zeros((M, M),dtype='complex')
uvPlain1 = NP.zeros((M, M),dtype='complex')

frequency = 8e9
lamb = 3e8 / frequency

ephDate = '2020-05-09'
RAO = BadaryRAO(ephDate)
declination = RAO.declination
declGrad = int(NP.rad2deg(declination) + .5)
antA = 64
antB = 65
mPerPix = 0.03
pixPerU = 1/2
pixPerM = int(1 / mPerPix + .5)
mRadius = 0.9
singleRadiusU = mRadius / lamb * pixPerU
areaA = NP.zeros((N, N))
areaB = NP.zeros((N, N))
pixPerRad = 16*N / pixPerU
arcsecRadius = 1080
radRadius = NP.deg2rad(arcsecRadius/3600)
qSunRadius = int(radRadius * pixPerRad +0.5)
FOV = NP.deg2rad(pixPerRad * N / 3600)

hourAngle = NP.linspace(-75,75,L)

u0 = M//2
v0 = M//2

beamPatterns = NP.zeros((N,N,L), dtype = 'float32')
qSun = NP.zeros((N, N), dtype = 'float32')

for i in range(N):
    for j in range(N):
        x=i - N/2
        y=j - N/2
        if (NP.sqrt(x**2 + y**2) < singleRadiusU):
            areaA[i, j] = 1.
            if (NP.sqrt((x - 3)**2 + y**2) > singleRadiusU):
                areaB[i, j] = 1.

for i in range(N):
    for j in range(N):
        x=i - N/2
        y=j - N/2
        if (NP.sqrt(x**2 + y**2) < qSunRadius):
            qSun[i, j] = 1.
            
for i in range(L):
    uvw = b2u.base2uvw(NP.deg2rad(hourAngle[i]), declination, antA, antB) / lamb * pixPerU
    u1 = M//2 +int(uvw[0] + .5)
    v1 = M//2 +int(uvw[1] + .5)
    uvPlain0[:,:] = 0
    uvPlain1[:,:] = 0
    uvPlain0[u0 - N//2:u0+N//2,v0 - N//2:v0+N//2] = areaA*complex(0,-1)
    uvPlain0[u1 - N//2:u1+N//2,v1 - N//2:v1+N//2] = areaB
    uvPlain1[u0 - N//2:u0+N//2,v0 - N//2:v0+N//2] = areaA
    uvPlain1[u1 - N//2:u1+N//2,v1 - N//2:v1+N//2] = areaB*complex(0,-1)
    beamPatterns[:,:,i] = (fftBeam(fftConvolution(uvPlain0, uvPlain1).real).real)[M//2-N//2:M//2+N//2,M//2-N//2:M//2+N//2]

#beamHeader = fits.Header()
#beamHeader['FREQ'] = ('%f'%frequency)
#beamHeader['DECL'] = ('%f'%declination)
#hdu = fits.PrimaryHDU(data=beamPatterns, header=beamHeader)
#hduList = fits.HDUList(hdu)
#hduList.writeto('srhPattern_%d%02d%02d_%03d_%03d_%04d.fit'%(year,month,day,antA,antB,int(frequency/1e6)))
