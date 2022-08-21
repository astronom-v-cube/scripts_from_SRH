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

N = 256
L = 3000
pix=20 # arcsec per pixel
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
qSmoothSun = qSmoothSun[dL//2:dL//2+N,dL//2:dL//2+N]

qSmoothSunVisibility = NP.zeros(L)
for i in range(L):
    uv_arg = 2.*NP.pi*i*FOV*NP.linspace(-.5,.5,N)
#    qSmoothSunVisibility[i] = (NP.cos(uv_arg)*qSmoothSun).sum()
    qSmoothSunVisibility[i] = (NP.cos(uv_arg)*qSun).sum()
qSmoothSunVisibility /= qSmoothSunVisibility.max()
qSmoothSunVisibility = NP.abs(qSmoothSunVisibility)
#7.5 

dL = 101
dW = 5
L0 = 122
nL = 1

deg0 = -50.
deg1 = 0.
decl_deg = -13
decl = NP.deg2rad(decl_deg)
hourAngle = NP.deg2rad(NP.linspace(deg0,deg1,100))
qResponse = NP.zeros(hourAngle.shape[0])
baselines = NP.zeros(512)
visVal = NP.zeros(512)
uv_val_vs_time = NP.zeros((hourAngle.shape[0],512))
freq = 7.0e9

fig, ax = PL.subplots(nrows=2, figsize=(12,12))
pl1, = ax[0].plot([], [], '.', color='red', animated=True)
pl2, = ax[0].plot([], [], '.', color='blue', animated=True)
pl3, = ax[0].plot([], [], color='blue', animated=True)
pl4, = ax[1].plot([], [], '.', color='green', animated=True)
ax[0].set_title('declination %d degrees, frequency %d MHz, lobe %d' % (decl_deg, freq/1e6, nL))
ax[0].set_xlim(0,3000)
ax[0].set_ylim(5,15)
ax[0].set_ylabel('baseline number')
ax[0].set_xlabel('spatial frequency')
ax[0].plot(NP.log(qSmoothSunVisibility) + 15)
ax[0].grid()

ax[1].set_xlim(deg0,deg1)
ax[1].set_ylim(0.000,0.02)
ax[1].set_ylabel('hour angle')
ax[1].set_xlabel('correlation')

qSmoothSunVisibility[0:130] = 0.
for i in range(25):
    qSmoothSunVisibility[(L0 + i*dL) - dW:(L0 + i*dL) + dW] = 0.
    
def updateQSunVis(hAngleIndex):
    baselines[:] = 0
    qResponse[hAngleIndex] = 0.
    uv = NP.zeros((32, 16),dtype=int)
    uv_val = NP.zeros((32, 16))
    lamb = 3e8/freq
    for ant2 in range(16):
        for ant1 in range(32):
            uvw = b2u.base2uvw(hourAngle[hAngleIndex], decl, 49+ant1, 192-ant2)
            baseline = NP.sqrt(uvw[0]**2 + uvw[1]**2) / lamb
            uv[ant1, ant2] = int(baseline + .5)
            visVal[ant2*32 + ant1] = qSmoothSunVisibility[uv[ant1, ant2]]
            baselines[ant2*32 + ant1] = baseline*qSmoothSunVisibility[uv[ant1, ant2]]*100
    uv_val = qSmoothSunVisibility[uv]
    uv_val_vs_time[hAngleIndex] = visVal
    qResponse[hAngleIndex] = uv_val.mean()
      
    pl1.set_data(uv,NP.log(uv_val) + 15)
#    baselineHist = NP.histogram(baselines,bins=100)
#    pl2.set_data(baselineHist[1][1:],baselineHist[0]*.1)
#    pl3.set_data(baselineHist[1][1:],baselineHist[0]*.1)
    pl4.set_data(NP.linspace(deg0,deg1,hourAngle.shape[0]), qResponse)

    return pl1, pl2, pl3, pl4


ani = animation.FuncAnimation(fig, updateQSunVis, frames=hourAngle.shape[0], blit=True, repeat=False)
#ani.save('srh_cp_hist_%d_%d_%d.mp4' % (decl_deg, freq/1e6, arcsecRadius))

PL.show()
