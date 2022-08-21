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
import matplotlib.animation as animation

N = 1024
L = 12000
pix=5 # arcsec per pixel
arcsecRadius = 1020
degRadius = NP.deg2rad(arcsecRadius/3600)
radius = int(arcsecRadius/pix +0.5)
FOV = NP.deg2rad(pix * N / 3600)
dishDiameter = 1.8

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
            
dL = 2*( 5//2) + 1
kern = NP.ones((dL,dL))
arg_x = NP.linspace(-1.,1,dL)
arg_y = NP.linspace(-1.,1,dL)
xx, yy = NP.meshgrid(arg_x, arg_y)

gKern =   NP.exp(-0.5*(xx**2 + yy**2))
qSmoothSun = scipy.signal.fftconvolve(qSun,gKern) / dL**2
#qSmoothSun = qSmoothSun[dL//2:dL//2+N,dL//2:dL//2+N]
qSmoothSun = qSun

qSmoothSunVisibility = NP.zeros(L)
for i in range(L):
    uv_arg = 2.*NP.pi*i*FOV*NP.linspace(-.5,.5,N)
#    qSmoothSunVisibility[i] = (NP.cos(uv_arg)*qSmoothSun).sum()
    qSmoothSunVisibility[i] = (NP.cos(uv_arg)*qSun).sum()
qSmoothSunVisibility /= qSmoothSunVisibility.max()
qSmoothSunVisibility = NP.abs(qSmoothSunVisibility)
#7.5 
#qSmoothSunVisibility[0:10] = 0.
#qSmoothSunVisibility[350:] = 0.

#qSmoothSunVisibility[0:120] = 0.
#qSmoothSunVisibility[450:] = 0.

#qSmoothSunVisibility[0:850] = 0.
#qSmoothSunVisibility[2850:] = 0.

deg0 = -45.
deg1 = 0.
decl_deg = 0
decl = NP.deg2rad(decl_deg)
hourAngle = NP.deg2rad(NP.linspace(deg0,deg1,200))
qResponse = NP.zeros(hourAngle.shape[0])
EW_N = 128
S_N = 64
baselines = NP.zeros(EW_N*S_N)
freq = 7.5e9

fig, ax = PL.subplots(nrows=2, figsize=(12,12))
pl0, = ax[0].plot([], [], 'ro', color='red', animated=True)
pl1, = ax[1].plot([], [], '.', color='blue', animated=True)
ax[0].set_title('declination %d degrees, frequency %d MHz' % (decl_deg, freq/1e6))
ax[0].set_xlim(0,L)
ax[0].set_ylim(-15,1)
ax[0].set_ylabel('ln(visibility)')
ax[0].set_xlabel('spatial frequency (base/lambda)')
ax[0].plot(NP.log(qSmoothSunVisibility))
ax[1].set_xlim(deg0,deg1)
ax[1].set_ylim(0.00,.002)
ax[1].set_ylabel('correlation coefficient')
ax[1].set_xlabel('houe angle (degrees)')

def updateQSunVis(hAngleIndex):
    baselines[:] = 0
    qResponse[hAngleIndex] = 0.
    uv = NP.zeros((EW_N, S_N),dtype=int)
    uv_val = NP.zeros((EW_N, S_N))
    lamb = 3e8/freq
    for ant2 in range(S_N):
        for ant1 in range(EW_N):
#            uvw = b2u.base2uvw(hourAngle[hAngleIndex], decl, 49+ant1, 192-ant2)
            uvw = b2u.base2uvw(hourAngle[hAngleIndex], decl, 1+ant1, 192-ant2)
            baseline = NP.sqrt(uvw[0]**2 + uvw[1]**2) / lamb
            uv[ant1, ant2] = int(baseline + .5)
            baselines[ant2*EW_N + ant1] = baseline
    uv_val = qSmoothSunVisibility[uv]
    qResponse[hAngleIndex] = qSmoothSunVisibility[uv].mean()
      
    pl0.set_data(uv,NP.log(uv_val))
    pl1.set_data(NP.linspace(deg0,deg1,hourAngle.shape[0]), qResponse)
#    baselineHist = NP.histogram(baselines,bins=128)
#    pl1.set_data(baselineHist[1][1:],NP.log(baselineHist[0] * qSmoothSunVisibility[NP.array(baselineHist[1][1:], dtype='int')]))
#    print(hourAngle[hAngleIndex], qResponse[hAngleIndex])
    return pl0, pl1


ani = animation.FuncAnimation(fig, updateQSunVis, frames=hourAngle.shape[0], blit=True, repeat=False)
#ani.save('srh_cp_%d_%d.mpg' % (decl_deg, freq/1e6))

PL.show()
