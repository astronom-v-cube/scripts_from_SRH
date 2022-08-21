#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 01:20:56 2018

@author: sergey
"""
import pylab as PL
import numpy as NP
import scipy.signal
import base2uvw as b2u
import scipy.special
import datetime as DT
from BadaryRAO import BadaryRAO

N = 256
L = 3000
pix=20 # arcsec per pixel
arcsecRadius = 1020
degRadius = NP.deg2rad(arcsecRadius/3600)
radius = int(arcsecRadius/pix +0.5)
FOV = NP.deg2rad(pix * N / 3600)
dishDiameter = 1.8
freq = 6.8e9

qSun = NP.zeros((N, N))

uv_arg = NP.deg2rad(pix * N / 3600)*NP.linspace(-1.,1.,N)
gx, gy = NP.meshgrid(uv_arg,uv_arg)
#singleBeam = NP.exp(-(1000.*gx*gx + 1000.*gy*gy))#/np.deg2rad(2.*(1. - freq/64))**2)

for i in range(N):
    for j in range(N):
        x=i - N/2
        y=j - N/2
        if (NP.sqrt(x**2 + y**2) < radius):
            qSun[i, j] = 1.
            
#qSun *= singleBeam
            
dL = 2*( 1//2) + 1
kern = NP.ones((dL,dL))
arg_x = NP.linspace(-1.,1,dL)
arg_y = NP.linspace(-1.,1,dL)
xx, yy = NP.meshgrid(arg_x, arg_y)

gKern =   NP.exp(-0.5*(xx**2 + yy**2))
qSmoothSun = scipy.signal.fftconvolve(qSun,gKern) / dL**2
qSmoothSun = qSmoothSun[dL//2:dL//2+N,dL//2:dL//2+N]

qSmoothSunVisibility = NP.zeros(L)
qSunJ1Visibility = NP.zeros(L)
cos_uv_arg = NP.zeros((L,N))
for i in range(L):
    cos_uv_arg[i] = NP.cos(2.*NP.pi*i*FOV*NP.linspace(-.5,.5,N))
    qSmoothSunVisibility[i] = NP.abs((cos_uv_arg[i]*qSmoothSun).mean())
    j_arg = 2*NP.pi*i*degRadius
    qSunJ1Visibility[i] = NP.abs(NP.pi*scipy.special.jv(1,j_arg)/(j_arg))
qSmoothSunVisibility /= qSmoothSunVisibility.max()
qSunJ1Visibility[0] = qSunJ1Visibility[1] * qSmoothSunVisibility[0] / qSmoothSunVisibility[1]
qSunJ1Visibility /= qSunJ1Visibility.max()


deg0 = -45.
deg1 = 45.
decl_deg = -23
decl = NP.deg2rad(NP.linspace(-23,23,365))
hourAngle = NP.deg2rad(NP.linspace(deg0,deg1,200))
qJ1Response = NP.zeros((decl.shape[0], hourAngle.shape[0]))
qSmoothResponse = NP.zeros((decl.shape[0], hourAngle.shape[0]))

fig, ax = PL.subplots(nrows=1, figsize=(12,12))

ax.set_title('frequency %d MHz, radius %d"' % (freq/1e6, arcsecRadius))
ax.set_ylim(0, 0.03)
ax.set_ylabel('correlation')
ax.set_xlabel('hour angle')

uv = NP.zeros((32, 16),dtype=int)
lamb = 3e8/freq

date = DT.date(2016,7,1)
#for declIndex in range(decl.shape[0]):
#    rao = BadaryRAO(date)
#    sunRadius = rao.sunObject.radius
#    print(date,NP.rad2deg(sunRadius)*60)
#    
#    for i in range(L):
#        j_arg = 2*NP.pi*i*sunRadius
#        qSunJ1Visibility[i] = NP.abs(NP.pi*scipy.special.jv(1,j_arg)/(j_arg))
#    qSunJ1Visibility[0] = qSunJ1Visibility[1] * qSmoothSunVisibility[0] / qSmoothSunVisibility[1]
#    qSunJ1Visibility /= qSunJ1Visibility.max()
#    qSunJ1Visibility[0:100] = 0.
#    
#    for hAngleIndex in range(hourAngle.shape[0]):
#        for ant2 in range(16):
#            for ant1 in range(32):
#                uvw = b2u.base2uvw(hourAngle[hAngleIndex], rao.declination, 49+ant1, 192-ant2)
#                baseline = NP.sqrt(uvw[0]**2 + uvw[1]**2) / lamb
#                uv[ant1, ant2] = int(baseline + .5)
#        qJ1Response[declIndex, hAngleIndex] = qSunJ1Visibility[uv].mean()
#    date += DT.timedelta(days = 1)
#
for declIndex in range(decl.shape[0]):
    rao = BadaryRAO(date)
    sunRadius = rao.sunObject.radius
    print(date,NP.rad2deg(sunRadius)*60)
    
    _2pi_sunRadius = 2*NP.pi*sunRadius
    for hAngleIndex in range(hourAngle.shape[0]):
        baselineCount = 0
        for ant1 in range(32):
            uvw = b2u.base2uvwNext(hourAngle[hAngleIndex], rao.declination, 49+ant1)
            for ant2 in range(16):
                baseline = NP.sqrt(uvw[ant2][0]**2 + uvw[ant2][1]**2) / lamb
                if baseline > 100.:
                    j_arg = _2pi_sunRadius*baseline
                    qJ1Response[declIndex, hAngleIndex] += NP.abs(NP.pi*scipy.special.jv(1,j_arg)/(j_arg))
                    baselineCount += 1
        qJ1Response[declIndex, hAngleIndex] /= baselineCount
    date += DT.timedelta(days = 1)

for i in range(decl.shape[0]):
    ax.plot(hourAngle, qJ1Response[i], '.')
