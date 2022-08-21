#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 01:20:56 2018

@author: sergey
"""
from astropy.io import fits
import pylab as PL
import numpy as NP
import scipy.signal
import datetime as DT
import ephem
import base2uvw as b2u
import scipy.special

N = 256
L = 3000
pix=20 # arcsec per pixel
arcsecRadius = 1080
degRadius = NP.deg2rad(arcsecRadius/3600)
radius = int(arcsecRadius/pix +0.5)
FOV = NP.deg2rad(pix * N / 3600)

qSun = NP.zeros((N, N))
qSunHalf = NP.zeros((N, N))

uv_arg = NP.deg2rad(pix * N / 3600)*NP.linspace(-1.,1.,N)
gx, gy = NP.meshgrid(uv_arg,uv_arg)
singleBeam = NP.exp(-(1000.*gx*gx + 1000.*gy*gy))#/np.deg2rad(2.*(1. - freq/64))**2)

for i in range(N):
    for j in range(N):
        x=i - N/2
        y=j - N/2
        if (NP.sqrt(x**2 + y**2) < radius and NP.sqrt(x**2 + y**2) > 0.95*radius and NP.abs(y) > NP.abs(2*x)):
            qSun[i, j] = 1.
        elif (NP.sqrt(x**2 + y**2) < radius):
            qSun[i, j] = 1.
            
for i in range(N):
    for j in range(N):
        x=i - N/2
        y=j - N/2
        if (NP.sqrt(x**2 + y**2) < radius and x > 0):
            qSunHalf[i, j] = 1.

qSun *= singleBeam
            
dL = 1*( 15//2) + 1
kern = NP.ones((dL,dL))
cqSun = scipy.signal.fftconvolve(qSun,kern) / dL**2
cqSun = cqSun[dL//2:dL//2+N,dL//2:dL//2+N]

qSunVisibility1D = NP.zeros(L)
qSunVisibility2D = NP.zeros(L)
qSunJ1 = NP.zeros(L)
qSunSinc = NP.zeros(L)
cqSunVisibility2D = NP.zeros(L)
for i in range(L):
    uv_arg = 2.*NP.pi*i*FOV*NP.linspace(-.5,.5,N)
    j_arg = 2*NP.pi*i*degRadius
    qSunVisibility1D[i] = (NP.cos(uv_arg)*qSun[N//2,:]).sum()
    qSunVisibility2D[i] = (NP.cos(uv_arg)*qSun).sum()
    cqSunVisibility2D[i] = (NP.cos(uv_arg)*cqSun).sum()
    qSunJ1[i] = NP.pi*scipy.special.jv(1,j_arg)/(j_arg)
    qSunSinc[i] = scipy.special.sinc(j_arg/NP.pi)

fig = PL.figure(figsize=(12,6));
sp0 = fig.add_subplot(3,2,1,title='Sharp Sun')
sp1 = fig.add_subplot(3,2,2,title='Smooth Sun')
sp2 = fig.add_subplot(3,2,3)
sp3 = fig.add_subplot(3,2,4)
sp4 = fig.add_subplot(3,2,5)
sp5 = fig.add_subplot(3,2,6)
sp2.set_xlim((0,3000))
sp3.set_xlim((0,3000))
sp2.set_ylim((-0.2,1.2))
sp3.set_ylim((-0.2,1.2))
#sp4.set_ylim((0.,1e5))
#sp5.set_ylim((0.,1e5))

sp0.imshow(qSun)
sp1.imshow(cqSun)
#PL.plot(qSunVisibility1D/qSunVisibility1D.max())
#PL.plot(scipy.special.sinc(NP.linspace(0,21,L)))
sp2.plot(qSunVisibility2D/qSunVisibility2D.max())
sp3.plot(cqSunVisibility2D/cqSunVisibility2D.max())
#PL.plot(qSunJ1/qSunJ1[1:].max())
#PL.plot(qSunSinc/qSunSinc[1:].max())
sp2.grid()
sp3.grid()
uv = NP.zeros(32,dtype=int)
hourAngle = NP.linspace(-1.5,1.5,160)
qResponse = NP.zeros((32,hourAngle.shape[0]))
cqResponse = NP.zeros((32,hourAngle.shape[0]))
decl = NP.deg2rad(-10)
for hh in range(hourAngle.shape[0]):
    for freq in range(5):
        lamb = 3e8/(4e9 + freq*(4e9/5))
        for ant2 in range(16):
            for ant1 in range(32):
#            for ant1 in [0,1,2,3,4,5,6,7,  24,25,26,27,28,29,30,31]:
#            for ant1 in [8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]:
                uvw = b2u.base2uvw(hourAngle[hh], decl, 49+ant1, 192-ant2)
                uv[ant1] = int(NP.sqrt(uvw[0]**2 + uvw[1]**2)/lamb + .5)
            if freq == 0 or freq == 4 and hh == 1:
                sp2.plot(uv,NP.ones(32)*.07*ant2,'.',markersize=1+0.5*freq/5)
                sp3.plot(uv,NP.ones(32)*.07*ant2,'.',markersize=1+0.5*freq/5)
            qResponse[freq,hh] += qSunVisibility2D[uv].sum()
            cqResponse[freq,hh] += cqSunVisibility2D[uv].sum()
    
#for i in range(1,40):
#    qResponse[:,i] /= qResponse[:,0]
#    cqResponse[:,i] /= cqResponse[:,0]

#sp4.imshow(qResponse[:,1:], origin='lower')
#sp5.imshow(cqResponse[:,1:], origin='lower')
    
for i in range(5):
    sp4.plot(qResponse[i,:])
    sp5.plot(cqResponse[i,:])

qResponse_m0 = qResponse
cqResponse_m0 = cqResponse
